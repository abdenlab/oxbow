use std::io;
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;
use noodles::core::region::Interval;

use crate::alignment::model::BatchBuilder;
use crate::alignment::model::Push as _;
use crate::util::query::BgzfChunkReader;

/// A record batch iterator yielding SAM or BAM records that intersect a genomic range.
pub struct BatchIterator<R> {
    reader: R,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
    header: noodles::sam::Header,
    reference_sequence_id: usize,
    interval: Interval,
}

impl<R> BatchIterator<R> {
    pub fn new(
        reader: R,
        header: noodles::sam::Header,
        reference_sequence_id: usize,
        interval: Interval,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        Self {
            reader,
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
            header,
            reference_sequence_id,
            interval,
        }
    }
}

impl<R> RecordBatchReader for BatchIterator<R>
where
    Self: Iterator<Item = Result<RecordBatch, ArrowError>>,
{
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

impl<R> Iterator for BatchIterator<noodles::sam::io::Reader<BgzfChunkReader<R>>>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = noodles::sam::Record::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => {
                    match intersects(
                        &self.header,
                        &record,
                        self.reference_sequence_id,
                        self.interval,
                    ) {
                        Ok(true) => match self.builder.push(&record) {
                            Ok(_) => {
                                self.count += 1;
                                count += 1;
                            }
                            Err(e) => return Some(Err(e.into())),
                        },
                        Ok(false) => {}
                        Err(e) => return Some(Err(e.into())),
                    };
                }
                Err(e) => return Some(Err(e.into())),
            }
        }

        if count == 0 {
            None
        } else {
            let batch = self.builder.finish();
            Some(batch)
        }
    }
}

impl<R> Iterator for BatchIterator<noodles::bam::io::Reader<BgzfChunkReader<R>>>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let header = self.builder.header();
        let mut record = noodles::bam::Record::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => {
                    match intersects(&header, &record, self.reference_sequence_id, self.interval) {
                        Ok(true) => {
                            match self.builder.push(&record) {
                                Ok(_) => {
                                    self.count += 1;
                                    count += 1;
                                }
                                Err(e) => return Some(Err(e.into())),
                            };
                        }
                        Ok(false) => {}
                        Err(e) => return Some(Err(e.into())),
                    };
                }
                Err(e) => return Some(Err(e.into())),
            }
        }

        if count == 0 {
            None
        } else {
            let batch = self.builder.finish();
            Some(batch)
        }
    }
}

pub fn intersects(
    header: &noodles::sam::Header,
    record: &impl noodles::sam::alignment::Record,
    region_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    match (
        record.reference_sequence_id(header).transpose()?,
        record.alignment_start().transpose()?,
        record.alignment_end().transpose()?,
    ) {
        (Some(record_id), Some(start), Some(end)) => {
            let record_interval = (start..=end).into();
            Ok(region_id == record_id && region_interval.intersects(record_interval))
        }
        _ => Ok(false),
    }
}

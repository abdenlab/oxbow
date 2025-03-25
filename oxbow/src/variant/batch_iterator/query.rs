use std::io;
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;
use noodles::core::region::Interval;

use crate::util::query::BgzfChunkReader;
use crate::variant::model::BatchBuilder;
use crate::variant::model::Push as _;

/// An iterator over records of an indexed file that intersects a given region.
pub struct BatchIterator<R> {
    reader: R,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
    header: noodles::vcf::Header,
    reference_sequence_name: String,
    interval: Interval,
}

impl<R> BatchIterator<R> {
    pub fn new(
        reader: R,
        header: noodles::vcf::Header,
        reference_sequence_name: String,
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
            reference_sequence_name,
            interval,
        }
    }
}

impl<R> RecordBatchReader for BatchIterator<R>
where
    BatchIterator<R>: Iterator<Item = Result<RecordBatch, ArrowError>>,
{
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

impl<R> Iterator for BatchIterator<noodles::vcf::io::Reader<BgzfChunkReader<R>>>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = noodles::vcf::Record::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => {
                    match intersects(
                        &self.header,
                        &record,
                        &self.reference_sequence_name,
                        self.interval,
                    ) {
                        Ok(true) => match self.builder.push(&record) {
                            Ok(()) => {
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

impl<R> Iterator for BatchIterator<noodles::bcf::io::Reader<BgzfChunkReader<R>>>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = noodles::bcf::Record::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => {
                    match intersects(
                        &self.header,
                        &record,
                        &self.reference_sequence_name,
                        self.interval,
                    ) {
                        Ok(true) => match self.builder.push(&record) {
                            Ok(()) => {
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

pub fn intersects(
    header: &noodles::vcf::Header,
    record: &impl noodles::vcf::variant::Record,
    region_chrom: &str,
    region_interval: Interval,
) -> io::Result<bool> {
    let record_chrom = record.reference_sequence_name(header)?;
    let Some(start) = record.variant_start().transpose()? else {
        return Ok(false);
    };
    let end = record.variant_end(header)?;
    let record_interval = Interval::from(start..=end);
    Ok(record_chrom == region_chrom && record_interval.intersects(region_interval))
}

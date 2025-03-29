use std::io;
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;

use noodles::core::region::Interval;
use noodles::csi::binning_index;

use crate::gxf::model::BatchBuilder;
use crate::gxf::model::Push as _;
use crate::util::query::BgzfChunkReader;

/// An iterator over records of an indexed file that intersects a given region.
pub struct BatchIterator<R> {
    reader: R,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
    header: binning_index::index::Header,
    reference_sequence_id: usize,
    interval: Interval,
}

impl<R> BatchIterator<R> {
    pub fn new(
        reader: R,
        header: binning_index::index::Header,
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
            reference_sequence_id,
            interval,
            header,
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

impl<R> Iterator for BatchIterator<noodles::gtf::io::Reader<BgzfChunkReader<R>>>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buf = String::new();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            buf.clear();
            match self.reader.read_line(&mut buf) {
                Ok(0) => break,
                Ok(_) => {
                    let line: noodles::gtf::Line = match buf.parse() {
                        Ok(line) => line,
                        Err(e) => return Some(Err(ArrowError::ExternalError(e.into()))),
                    };
                    match line {
                        noodles::gtf::Line::Comment(_) => continue,
                        noodles::gtf::Line::Record(record) => {
                            match intersects_gtf(
                                &self.header,
                                &record,
                                self.reference_sequence_id,
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
                            }
                        }
                    };
                }
                Err(e) => return Some(Err(e.into())),
            };
        }

        if count == 0 {
            None
        } else {
            let batch = self.builder.finish();
            Some(batch)
        }
    }
}

impl<R> Iterator for BatchIterator<noodles::gff::io::Reader<BgzfChunkReader<R>>>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = noodles::gff::Line::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.read_line(&mut line) {
                Ok(0) => break,
                Ok(_) => {
                    match line.as_record() {
                        Some(Ok(record)) => {
                            match intersects_gff(
                                &self.header,
                                &record,
                                self.reference_sequence_id,
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
                            }
                        }
                        _ => continue,
                    };
                }
                Err(e) => return Some(Err(e.into())),
            };
        }

        if count == 0 {
            None
        } else {
            let batch = self.builder.finish();
            Some(batch)
        }
    }
}

pub fn intersects_gtf(
    header: &noodles::csi::binning_index::index::Header,
    record: &noodles::gtf::Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    let rname = record.reference_sequence_name().as_bytes();
    match (
        header.reference_sequence_names().get_index_of(rname),
        record.start(),
        record.end(),
    ) {
        (Some(id), start, end) => {
            let alignment_interval = (start..=end).into();
            Ok(id == reference_sequence_id && region_interval.intersects(alignment_interval))
        }
        _ => Ok(false),
    }
}

pub fn intersects_gff(
    header: &noodles::csi::binning_index::index::Header,
    record: &noodles::gff::Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> io::Result<bool> {
    let rname = record.reference_sequence_name().as_bytes();
    match (
        header.reference_sequence_names().get_index_of(rname),
        record.start(),
        record.end(),
    ) {
        (Some(id), Ok(start), Ok(end)) => {
            let alignment_interval = (start..=end).into();
            Ok(id == reference_sequence_id && region_interval.intersects(alignment_interval))
        }
        _ => Ok(false),
    }
}

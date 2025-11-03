use std::io::{BufRead, Read};
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;

use crate::alignment::model::BatchBuilder;
use crate::alignment::model::Push as _;

/// A record batch iterator yielding SAM or BAM records from a readable stream.
pub struct BatchIterator<R> {
    reader: R,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
}

impl<R> BatchIterator<R> {
    pub fn new(reader: R, builder: BatchBuilder, batch_size: usize, limit: Option<usize>) -> Self {
        Self {
            reader,
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
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

impl<R> Iterator for BatchIterator<noodles::sam::io::Reader<R>>
where
    R: BufRead,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = noodles::sam::Record::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => match self.builder.push(&record) {
                    Ok(()) => {
                        self.count += 1;
                        count += 1;
                    }
                    Err(e) => return Some(Err(e.into())),
                },
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

impl<R> Iterator for BatchIterator<noodles::bam::io::Reader<R>>
where
    R: Read,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = noodles::bam::Record::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => match self.builder.push(&record) {
                    Ok(()) => {
                        self.count += 1;
                        count += 1;
                    }
                    Err(e) => return Some(Err(e.into())),
                },
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

impl<R> Iterator for BatchIterator<noodles::cram::io::Reader<R>>
where
    R: Read,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut count = 0;

        let header = self.builder.header();
        let mut records = self.reader.records(&header);

        while count < self.batch_size && self.count < self.limit {
            match records.next() {
                Some(Ok(record)) => match self.builder.push(&record) {
                    Ok(()) => {
                        self.count += 1;
                        count += 1;
                    }
                    Err(e) => return Some(Err(e.into())),
                },
                Some(Err(e)) => {
                    // noodles::cram::io::Records iterator doesn't terminate gracefully on EOF
                    if e.kind() == std::io::ErrorKind::UnexpectedEof {
                        return None;
                    }
                    return Some(Err(e.into()));
                }
                None => break,
            }
        }

        if count == 0 {
            None
        } else {
            match self.builder.finish() {
                Ok(batch) => Some(Ok(batch)),
                Err(e) => Some(Err(e)),
            }
        }
    }
}

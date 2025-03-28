use std::io::BufRead;
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;

use crate::gxf::model::BatchBuilder;
use crate::gxf::model::Push as _;

/// An iterator over all records a GTF/GFF file starting at the current position.
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

impl<R> Iterator for BatchIterator<noodles::gtf::io::Reader<R>>
where
    R: BufRead,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut count = 0;
        let mut records = self.reader.records();

        while count < self.batch_size && self.count < self.limit {
            match records.next() {
                Some(Ok(record)) => match self.builder.push(&record) {
                    Ok(()) => {
                        self.count += 1;
                        count += 1;
                    }
                    Err(e) => return Some(Err(e.into())),
                },
                Some(Err(e)) => return Some(Err(e.into())),
                None => break,
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

impl<R> Iterator for BatchIterator<noodles::gff::io::Reader<R>>
where
    R: BufRead,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut count = 0;
        let mut lines = self.reader.lines();

        while count < self.batch_size && self.count < self.limit {
            match lines.next() {
                Some(Ok(line)) => {
                    match line.as_record() {
                        Some(Ok(record)) => match self.builder.push(&record) {
                            Ok(()) => {
                                self.count += 1;
                                count += 1;
                            }
                            Err(e) => return Some(Err(e.into())),
                        },
                        _ => continue,
                    };
                }
                Some(Err(e)) => return Some(Err(e.into())),
                None => break,
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

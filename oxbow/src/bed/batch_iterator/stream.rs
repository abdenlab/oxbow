use std::io::BufRead;
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;

use crate::bed::model::BatchBuilder;
use crate::bed::model::Push as _;

/// A record batch iterator yielding BED records from a readable stream.
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

impl<R> Iterator for BatchIterator<noodles::bed::io::Reader<3, R>>
where
    R: BufRead,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = noodles::bed::Record::<3>::default();
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

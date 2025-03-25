use std::io::BufRead;
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;
use noodles::fasta::record::Definition;
use noodles::fasta::record::Sequence;

use crate::sequence::model::BatchBuilder;
use crate::sequence::model::Push as _;

/// A record batch iterator yielding FASTA or FASTQ records from a readable stream.
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

impl<R> Iterator for BatchIterator<noodles::fastq::io::Reader<R>>
where
    R: BufRead,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = noodles::fastq::Record::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.read_record(&mut record) {
                Ok(0) => break,
                Ok(_) => {
                    match self.builder.push(&record) {
                        Ok(()) => {
                            self.count += 1;
                            count += 1;
                        }
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

impl<R> Iterator for BatchIterator<noodles::fasta::io::Reader<R>>
where
    R: BufRead,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line_buf = String::new();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            line_buf.clear();
            let definition = match self.reader.read_definition(&mut line_buf) {
                Ok(0) => break,
                Ok(_) => match line_buf.parse::<Definition>() {
                    Ok(def) => def,
                    Err(e) => return Some(Err(ArrowError::ExternalError(Box::new(e)))),
                },
                Err(e) => return Some(Err(e.into())),
            };

            let mut sequence_buf = Vec::new();
            match self.reader.read_sequence(&mut sequence_buf) {
                Ok(_) => {
                    let record =
                        noodles::fasta::Record::new(definition, Sequence::from(sequence_buf));
                    match self.builder.push(&record) {
                        Ok(()) => {
                            self.count += 1;
                            count += 1;
                        }
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

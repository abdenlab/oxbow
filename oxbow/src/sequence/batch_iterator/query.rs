use std::io::{BufRead, Seek};
use std::vec::IntoIter;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;
use noodles::core::region::Region;

use crate::batch::{Push as _, RecordBatchBuilder as _};
use crate::sequence::model::BatchBuilder;

/// A record batch iterator that slices sequences from an indexed FASTA file.
pub struct BatchIterator<R> {
    reader: R,
    index: noodles::fasta::fai::Index,
    regions: IntoIter<Region>,
    builder: BatchBuilder,
    batch_size: usize,
}

impl<R> BatchIterator<noodles::fasta::io::Reader<R>>
where
    R: BufRead + Seek,
{
    pub fn new(
        reader: noodles::fasta::io::Reader<R>,
        index: noodles::fasta::fai::Index,
        regions: Vec<Region>,
        builder: BatchBuilder,
        batch_size: usize,
    ) -> Self {
        Self {
            reader,
            index,
            regions: regions.into_iter(),
            batch_size,
            builder,
        }
    }
}

impl<R> RecordBatchReader for BatchIterator<R>
where
    Self: Iterator<Item = Result<RecordBatch, ArrowError>>,
{
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        self.builder.schema()
    }
}

impl<R> Iterator for BatchIterator<noodles::fasta::io::Reader<R>>
where
    R: BufRead + Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut count = 0;

        while count < self.batch_size {
            let region = self.regions.next();
            match region {
                Some(region) => match self.reader.query(&self.index, &region) {
                    Ok(record) => {
                        match self.builder.push(&record) {
                            Ok(()) => {
                                count += 1;
                            }
                            Err(e) => return Some(Err(e.into())),
                        };
                    }
                    Err(e) => return Some(Err(e.into())),
                },
                None => break, // No more regions to read
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

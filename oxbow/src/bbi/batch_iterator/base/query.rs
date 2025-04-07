use std::io::{Read, Seek};
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;
use bigtools::{BBIReadError, BigBedRead, BigWigRead};
use noodles::core::Region;

use crate::bbi::model::base::Push as _;
use crate::bbi::model::base::{BatchBuilder, BigBedRecord, BigWigRecord};

/// A record batch iterator yielding BigWig records that intersect a genomic range.
pub struct BigWigBatchIterator {
    entries: Box<dyn Iterator<Item = Result<(String, bigtools::Value), BBIReadError>> + Send>,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
}

impl BigWigBatchIterator {
    pub fn new<R: Read + Seek + Send + 'static>(
        reader: BigWigRead<R>,
        region: Region,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        let chrom = region.name().to_string();
        let length = reader
            .chroms()
            .iter()
            .find(|item| item.name == chrom)
            .map(|item| item.length)
            .unwrap();
        let interval = region.interval();
        let start = interval
            .start()
            .map(|pos| usize::from(pos) as u32)
            .unwrap_or(1)
            - 1;
        let end = interval
            .end()
            .map(|pos| usize::from(pos) as u32)
            .unwrap_or_else(|| length);
        let iter = reader
            .get_interval_move(&chrom, start, end)
            .unwrap()
            .map(move |result| result.map(|value| (chrom.clone(), value)));
        Self {
            entries: Box::new(iter),
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
        }
    }
}

impl Iterator for BigWigBatchIterator {
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.entries.next() {
                Some(Ok((chrom, entry))) => {
                    let record = BigWigRecord::new(&chrom, &entry);
                    match self.builder.push(&record) {
                        Ok(()) => {
                            self.count += 1;
                            count += 1;
                        }
                        Err(e) => return Some(Err(e.into())),
                    }
                }
                Some(Err(e)) => return Some(Err(ArrowError::ExternalError(e.into()))),
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

impl RecordBatchReader for BigWigBatchIterator {
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

/// A record batch iterator yielding BigBed records that intersect a genomic range.
pub struct BigBedBatchIterator {
    entries: Box<dyn Iterator<Item = Result<(String, bigtools::BedEntry), BBIReadError>> + Send>,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
}

impl BigBedBatchIterator {
    pub fn new<R: Read + Seek + Send + 'static>(
        reader: BigBedRead<R>,
        region: Region,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        let chrom = region.name().to_string();
        let length = reader
            .chroms()
            .iter()
            .find(|item| item.name == chrom)
            .map(|item| item.length)
            .unwrap();
        let interval = region.interval();
        let start = interval
            .start()
            .map(|pos| usize::from(pos) as u32)
            .unwrap_or(1)
            - 1;
        let end = interval
            .end()
            .map(|pos| usize::from(pos) as u32)
            .unwrap_or_else(|| length);
        let iter = reader
            .get_interval_move(&chrom, start, end)
            .unwrap()
            .map(move |result| result.map(|value| (chrom.clone(), value)));
        Self {
            entries: Box::new(iter),
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
        }
    }
}

impl Iterator for BigBedBatchIterator {
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.entries.next() {
                Some(Ok((chrom, entry))) => {
                    let record = BigBedRecord::new(&chrom, &entry);
                    match self.builder.push(&record) {
                        Ok(()) => {
                            self.count += 1;
                            count += 1;
                        }
                        Err(e) => return Some(Err(e.into())),
                    }
                }
                Some(Err(e)) => return Some(Err(ArrowError::ExternalError(e.into()))),
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

impl RecordBatchReader for BigBedBatchIterator {
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

use std::io::{Read, Seek};
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;
use bigtools::{BigBedRead, BigWigRead, ZoomIntervalError, ZoomRecord};
use noodles::core::region::Region;

use crate::bbi::model::zoom::Push as _;
use crate::bbi::model::zoom::{BBIZoomRecord, BatchBuilder};

pub struct BatchIterator {
    builder: BatchBuilder,
    entries: Box<dyn Iterator<Item = Result<(String, ZoomRecord), ZoomIntervalError>> + Send>,
    batch_size: usize,
    limit: usize,
    count: usize,
}

impl BatchIterator {
    pub fn from_bigwig<R: Read + Seek + Send + 'static>(
        reader: BigWigRead<R>,
        region: Region,
        level: u32,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        let chrom = region.name().to_string();
        let chrominfo: Vec<_> = reader
            .chroms()
            .iter()
            .map(|item| (item.name.clone(), item.length))
            .collect();
        let length = chrominfo
            .iter()
            .find(|(name, _)| *name == chrom)
            .map(|(_, length)| length)
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
            .unwrap_or_else(|| *length);
        let iter = reader
            .get_zoom_interval_move(&chrom, start, end, level)
            .unwrap()
            .map(move |result| {
                result
                    .map(|record| (chrom.clone(), record))
                    .map_err(ZoomIntervalError::BBIReadError)
            });
        Self {
            builder,
            entries: Box::new(iter),
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
        }
    }

    pub fn from_bigbed<R: Read + Seek + Send + 'static>(
        reader: BigBedRead<R>,
        region: Region,
        level: u32,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        let chrom = region.name().to_string();
        let chrominfo: Vec<_> = reader
            .chroms()
            .iter()
            .map(|item| (item.name.clone(), item.length))
            .collect();
        let length = chrominfo
            .iter()
            .find(|(name, _)| *name == chrom)
            .map(|(_, length)| length)
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
            .unwrap_or_else(|| *length);
        let iter = reader
            .get_zoom_interval_move(&chrom, start, end, level)
            .unwrap()
            .map(move |result| {
                result
                    .map(|record| (chrom.clone(), record))
                    .map_err(ZoomIntervalError::BBIReadError)
            });
        Self {
            builder,
            entries: Box::new(iter),
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
        }
    }
}

impl Iterator for BatchIterator {
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.entries.next() {
                Some(Ok((chrom, entry))) => {
                    let record = BBIZoomRecord::new(&chrom, &entry);
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

impl RecordBatchReader for BatchIterator {
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

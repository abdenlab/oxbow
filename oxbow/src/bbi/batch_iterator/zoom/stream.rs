use std::io::{Read, Seek};
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;
use bigtools::{
    BigBedRead, BigWigRead, ChromInfo, ZoomIntervalError, ZoomIntervalIter, ZoomRecord,
};

use crate::bbi::model::zoom::Push as _;
use crate::bbi::model::zoom::{BBIZoomRecord, BatchBuilder};

enum Either<Left, Right> {
    Left(Left),
    Right(Right),
}

type ZoomIteratorState<B> = Either<B, (String, ZoomIntervalIter<B, B>)>;

/// Chains together zoom record iterators from all chromosomes of a BigWig or BigBed file.
struct BBIZoomRecords<B> {
    chroms: Vec<ChromInfo>,
    level: u32,
    state: Option<ZoomIteratorState<B>>,
}

impl<R> BBIZoomRecords<BigWigRead<R>>
where
    R: Read + Seek,
{
    pub fn from_bigwig(reader: BigWigRead<R>, level: u32) -> Self {
        let mut chroms = reader.chroms().to_vec();
        chroms.reverse();
        Self {
            chroms,
            level,
            state: Some(Either::Left(reader)),
        }
    }
}

impl<R> BBIZoomRecords<BigBedRead<R>>
where
    R: Read + Seek,
{
    pub fn from_bigbed(reader: BigBedRead<R>, level: u32) -> Self {
        let mut chroms = reader.chroms().to_vec();
        chroms.reverse();
        Self {
            chroms,
            level,
            state: Some(Either::Left(reader)),
        }
    }
}

impl<R> Iterator for BBIZoomRecords<BigWigRead<R>>
where
    R: Read + Seek,
{
    type Item = Result<(String, ZoomRecord), ZoomIntervalError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.state.take() {
                Some(Either::Left(bigwig)) => {
                    let info = self.chroms.pop()?;
                    let chrom = info.name;
                    let length = info.length;
                    let iter = match bigwig.get_zoom_interval_move(&chrom, 0, length, self.level) {
                        Ok(iter) => iter,
                        Err(e) => return Some(Err(e)),
                    };
                    self.state = Some(Either::Right((chrom, iter)));
                }
                Some(Either::Right((chrom, mut iter))) => match iter.next() {
                    Some(result) => {
                        let curr_chrom = chrom.clone();
                        self.state = Some(Either::Right((chrom, iter)));
                        let result = result
                            .map(|record| (curr_chrom, record))
                            .map_err(ZoomIntervalError::BBIReadError);
                        return Some(result);
                    }
                    None => {
                        let bigwig: BigWigRead<R> = iter.into();
                        self.state = Some(Either::Left(bigwig));
                    }
                },
                None => return None,
            }
        }
    }
}

impl<R> Iterator for BBIZoomRecords<BigBedRead<R>>
where
    R: Read + Seek,
{
    type Item = Result<(String, ZoomRecord), ZoomIntervalError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.state.take() {
                Some(Either::Left(bigbed)) => {
                    let info = self.chroms.pop()?;
                    let chrom = info.name;
                    let length = info.length;
                    let iter = match bigbed.get_zoom_interval_move(&chrom, 0, length, self.level) {
                        Ok(iter) => iter,
                        Err(e) => return Some(Err(e)),
                    };
                    self.state = Some(Either::Right((chrom, iter)));
                }
                Some(Either::Right((chrom, mut iter))) => match iter.next() {
                    Some(result) => {
                        let curr_chrom = chrom.clone();
                        self.state = Some(Either::Right((chrom, iter)));
                        let result = result
                            .map(|record| (curr_chrom, record))
                            .map_err(ZoomIntervalError::BBIReadError);
                        return Some(result);
                    }
                    None => {
                        let bigbed: BigBedRead<R> = iter.into();
                        self.state = Some(Either::Left(bigbed));
                    }
                },
                None => return None,
            }
        }
    }
}

/// A record batch iterator yielding batches from a [`BBIZoomRecords`] iterator.
pub struct BatchIterator<B> {
    entries: BBIZoomRecords<B>,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
}

impl<R> RecordBatchReader for BatchIterator<R>
where
    Self: Iterator<Item = Result<RecordBatch, ArrowError>>,
{
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

impl<R> BatchIterator<BigWigRead<R>>
where
    R: Read + Seek,
{
    pub fn from_bigwig(
        reader: BigWigRead<R>,
        level: u32,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        Self {
            entries: BBIZoomRecords::from_bigwig(reader, level),
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
        }
    }
}

impl<R> BatchIterator<BigBedRead<R>>
where
    R: Read + Seek,
{
    pub fn from_bigbed(
        reader: BigBedRead<R>,
        level: u32,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        Self {
            entries: BBIZoomRecords::from_bigbed(reader, level),
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
        }
    }
}

impl<R> Iterator for BatchIterator<BigWigRead<R>>
where
    R: Read + Seek,
{
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

impl<R> Iterator for BatchIterator<BigBedRead<R>>
where
    R: Read + Seek,
{
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

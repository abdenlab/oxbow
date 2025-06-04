use std::io::{Read, Seek};
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;
use bigtools::ChromInfo;
use bigtools::{BBIReadError, BigBedIntervalIter, BigBedRead, BigWigIntervalIter, BigWigRead};

use crate::bbi::model::base::Push as _;
use crate::bbi::model::base::{BatchBuilder, BigBedRecord, BigWigRecord};

enum Either<Left, Right> {
    Left(Left),
    Right(Right),
}

type BigWigIteratorState<R> = Either<BigWigRead<R>, (String, BigWigIntervalIter<R, BigWigRead<R>>)>;

/// Chains together record iterators from all chromosomes of a BigWig file.
struct BigWigRecords<R> {
    chroms: Vec<ChromInfo>,
    state: Option<BigWigIteratorState<R>>,
}

impl<R> BigWigRecords<R>
where
    R: Read + Seek,
{
    pub fn new(reader: BigWigRead<R>) -> Self {
        let mut chroms = reader.chroms().to_vec();
        chroms.reverse();
        Self {
            chroms,
            state: Some(Either::Left(reader)),
        }
    }
}

impl<R> Iterator for BigWigRecords<R>
where
    R: Read + Seek,
{
    type Item = Result<(String, bigtools::Value), BBIReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.state.take() {
                Some(Either::Left(bigwig)) => {
                    let info = self.chroms.pop()?;
                    let chrom = info.name;
                    let length = info.length;
                    let iter = match bigwig.get_interval_move(&chrom, 0, length) {
                        Ok(iter) => iter,
                        Err(e) => return Some(Err(e)),
                    };
                    self.state = Some(Either::Right((chrom, iter)));
                }
                Some(Either::Right((chrom, mut iter))) => match iter.next() {
                    Some(result) => {
                        let curr_chrom = chrom.clone();
                        self.state = Some(Either::Right((chrom, iter)));
                        return Some(result.map(|value| (curr_chrom, value)));
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

/// A record batch iterator yielding batches from a [`BigWigRecords`] iterator.
pub struct BigWigBatchIterator<R> {
    entries: BigWigRecords<R>,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
}

impl<R> BigWigBatchIterator<R>
where
    R: Read + Seek,
{
    pub fn new(
        reader: BigWigRead<R>,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        let iter = BigWigRecords::new(reader);
        Self {
            entries: iter,
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
        }
    }
}

impl<R> RecordBatchReader for BigWigBatchIterator<R>
where
    Self: Iterator<Item = Result<RecordBatch, ArrowError>>,
{
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

impl<R> Iterator for BigWigBatchIterator<R>
where
    R: Read + Seek,
{
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

type BigBedIteratorState<R> = Either<BigBedRead<R>, (String, BigBedIntervalIter<R, BigBedRead<R>>)>;

/// Chains together record iterators from all chromosomes of a BigBed file.
struct BigBedRecords<R> {
    chroms: Vec<ChromInfo>,
    state: Option<BigBedIteratorState<R>>,
}

impl<R> BigBedRecords<R>
where
    R: Read + Seek,
{
    pub fn new(reader: BigBedRead<R>) -> Self {
        let mut chroms = reader.chroms().to_vec();
        chroms.reverse();
        Self {
            chroms,
            state: Some(Either::Left(reader)),
        }
    }
}

impl<R> Iterator for BigBedRecords<R>
where
    R: Read + Seek,
{
    type Item = Result<(String, bigtools::BedEntry), BBIReadError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.state.take() {
                Some(Either::Left(bigbed)) => {
                    let info = self.chroms.pop()?;
                    let chrom = info.name;
                    let length = info.length;
                    let iter = match bigbed.get_interval_move(&chrom, 0, length) {
                        Ok(iter) => iter,
                        Err(e) => return Some(Err(e)),
                    };
                    self.state = Some(Either::Right((chrom, iter)));
                }
                Some(Either::Right((chrom, mut iter))) => match iter.next() {
                    Some(result) => {
                        let curr_chrom = chrom.clone();
                        self.state = Some(Either::Right((chrom, iter)));
                        return Some(result.map(|entry| (curr_chrom, entry)));
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

/// A record batch iterator yielding batches from a [`BigBedRecords`] iterator.
pub struct BigBedBatchIterator<R> {
    entries: BigBedRecords<R>,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
}

impl<R> BigBedBatchIterator<R>
where
    R: Read + Seek,
{
    pub fn new(
        reader: BigBedRead<R>,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
    ) -> Self {
        let iter = BigBedRecords::new(reader);
        Self {
            entries: iter,
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
        }
    }
}

impl<R> RecordBatchReader for BigBedBatchIterator<R>
where
    Self: Iterator<Item = Result<RecordBatch, ArrowError>>,
{
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

impl<R> Iterator for BigBedBatchIterator<R>
where
    R: Read + Seek,
{
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

use std::io::{BufRead, Read, Seek};
use std::sync::Arc;

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::record_batch::RecordBatchReader;

use crate::gxf::model::BatchBuilder;
use crate::gxf::model::Push as _;

/// A batch iterator over records until a given cursor position of the underlying reader is reached.
pub struct BatchIterator<R, I> {
    reader: R,
    builder: BatchBuilder,
    batch_size: usize,
    limit: usize,
    count: usize,
    stop: I,
}

impl<R, I> RecordBatchReader for BatchIterator<R, I>
where
    BatchIterator<R, I>: Iterator<Item = Result<RecordBatch, ArrowError>>,
{
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        Arc::new(self.builder.get_arrow_schema())
    }
}

/// Iterate until we reach the `stop` position of a seekable source.
impl<R, I> BatchIterator<R, I> {
    pub fn new(
        reader: R,
        builder: BatchBuilder,
        batch_size: usize,
        limit: Option<usize>,
        stop: I,
    ) -> Self {
        Self {
            reader,
            builder,
            batch_size,
            limit: limit.unwrap_or(usize::MAX),
            count: 0,
            stop,
        }
    }
}

impl<R> Iterator for BatchIterator<noodles::gtf::io::Reader<R>, u64>
where
    R: BufRead + Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buf = String::new();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.get_mut().stream_position() {
                Ok(pos) => {
                    if pos >= self.stop {
                        break;
                    }
                }
                Err(e) => return Some(Err(e.into())),
            };

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
                        noodles::gtf::Line::Record(record) => match self.builder.push(&record) {
                            Ok(()) => {
                                self.count += 1;
                                count += 1;
                            }
                            Err(e) => return Some(Err(e.into())),
                        },
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

impl<R> Iterator for BatchIterator<noodles::gff::io::Reader<R>, u64>
where
    R: BufRead + Seek,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = noodles::gff::Line::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            match self.reader.get_mut().stream_position() {
                Ok(pos) => {
                    if pos >= self.stop {
                        break;
                    }
                }
                Err(e) => return Some(Err(e.into())),
            };
            match self.reader.read_line(&mut line) {
                Ok(0) => break,
                Ok(_) => {
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

impl<R> Iterator
    for BatchIterator<
        noodles::gtf::io::Reader<noodles::bgzf::Reader<R>>,
        noodles::bgzf::VirtualPosition,
    >
where
    R: Read,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buf = String::new();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            if self.reader.get_ref().virtual_position() >= self.stop {
                break;
            }

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
                        noodles::gtf::Line::Record(record) => match self.builder.push(&record) {
                            Ok(()) => {
                                self.count += 1;
                                count += 1;
                            }
                            Err(e) => return Some(Err(e.into())),
                        },
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

impl<R> Iterator
    for BatchIterator<
        noodles::gff::io::Reader<noodles::bgzf::Reader<R>>,
        noodles::bgzf::VirtualPosition,
    >
where
    R: Read,
{
    type Item = Result<RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = noodles::gff::Line::default();
        let mut count = 0;

        while count < self.batch_size && self.count < self.limit {
            if self.reader.get_ref().virtual_position() >= self.stop {
                break;
            }

            match self.reader.read_line(&mut line) {
                Ok(0) => break,
                Ok(_) => {
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

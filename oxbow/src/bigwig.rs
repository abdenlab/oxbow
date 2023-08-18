use arrow::array::StringArray;
use arrow::array::{
    ArrayRef, Float32Array, Float32Builder, StringDictionaryBuilder, UInt32Array, UInt32Builder,
};
use arrow::datatypes::Int32Type;
use arrow::{error::ArrowError, record_batch::RecordBatch};
use bigtools::utils::reopen::ReopenableFile;
use bigtools::{BBIRead, BigWigRead};
use noodles::core::Region;
use std::io::{Read, Seek};
use std::sync::Arc;

use crate::batch_builder::{finish_batch, BatchBuilder};

/// A BigWig reader.
pub struct BigWigReader<R> {
    read: BigWigRead<R>,
}

pub struct BigWigRecord<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

impl BigWigReader<ReopenableFile> {
    /// Creates a BigWig reader from a given file path.
    pub fn new_from_path(path: &str) -> std::io::Result<Self> {
        let read = BigWigRead::open_file(path)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        Ok(Self { read })
    }
}

impl<R: Read + Seek> BigWigReader<R> {
    /// Creates a BigWig reader from a given file path.
    pub fn new(read: R) -> std::io::Result<Self> {
        let read = BigWigRead::open(read)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        Ok(Self { read })
    }

    /// Returns the records in the given region as Apache Arrow IPC.
    ///
    /// If the region is `None`, all records are returned.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use oxbow::bigwig::BigWigReader;
    ///
    /// let mut reader = BigWigReader::new_from_path("sample.bigWig").unwrap();
    /// let ipc = reader.records_to_ipc(Some("sq0:1-1000")).unwrap();
    /// ```
    pub fn records_to_ipc(&mut self, region: Option<&str>) -> Result<Vec<u8>, ArrowError> {
        let mut batch_builder = BigWigBatchBuilder::new(1024, &mut self.read)?;
        match region {
            Some(region) => {
                let region: Region = region.parse().unwrap();
                let chrom_name = region.name().to_owned();
                let (start, end) = match (region.interval().start(), region.interval().end()) {
                    (Some(start), Some(end)) => {
                        let start = start.get() as u32 - 1; // 1-based to 0-based
                        let end = end.get() as u32;
                        (start, end)
                    }
                    (Some(start), None) => {
                        let start = start.get() as u32 - 1; // 1-based to 0-based
                        let end = self
                            .read
                            .get_chroms()
                            .iter()
                            .find(|c| c.name == chrom_name)
                            .map(|c| c.length);
                        let end = end.ok_or_else(|| {
                            ArrowError::InvalidArgumentError("Invalid chromosome".to_string())
                        })?;
                        (start, end)
                    }
                    (None, Some(end)) => {
                        let start = 0;
                        let end = end.get() as u32;
                        (start, end)
                    }
                    (None, None) => {
                        let start = 0;
                        let end = self
                            .read
                            .get_chroms()
                            .iter()
                            .find(|c| c.name == chrom_name)
                            .map(|c| c.length);
                        let end = end.ok_or_else(|| {
                            ArrowError::InvalidArgumentError("Invalid chromosome".to_string())
                        })?;
                        (start, end)
                    }
                };
                let values = match self.read.get_interval(&chrom_name, start, end) {
                    Ok(v) => v,
                    Err(e) => {
                        return Err(ArrowError::ExternalError(Box::new(e)));
                    }
                };
                for value in values {
                    let v = value.unwrap();
                    let record = BigWigRecord {
                        chrom: &chrom_name,
                        start: v.start,
                        end: v.end,
                        value: v.value,
                    };
                    batch_builder.push(record);
                }
                finish_batch(batch_builder)
            }
            None => {
                // Can't use write_ipc, because we have separate iterators for each chrom
                let chroms = self.read.get_chroms().into_iter();
                for chrom in chroms {
                    let start = 0;
                    let end = chrom.length;
                    let values = self.read.get_interval(&chrom.name, start, end).unwrap();
                    for value in values {
                        let v = value.unwrap();
                        let record = BigWigRecord {
                            chrom: &chrom.name,
                            start: v.start,
                            end: v.end,
                            value: v.value,
                        };
                        batch_builder.push(record);
                    }
                }
                finish_batch(batch_builder)
            }
        }
    }
}

struct BigWigBatchBuilder {
    chrom: StringDictionaryBuilder<Int32Type>,
    start: UInt32Builder,
    end: UInt32Builder,
    value: Float32Builder,
}

impl BigWigBatchBuilder {
    pub fn new<R: Read + Seek>(
        capacity: usize,
        read: &mut BigWigRead<R>,
    ) -> Result<Self, ArrowError> {
        let chroms: Vec<String> = read.get_chroms().iter().map(|c| c.name.clone()).collect();
        let chroms: StringArray = StringArray::from(chroms);
        Ok(Self {
            chrom: StringDictionaryBuilder::<Int32Type>::new_with_dictionary(capacity, &chroms)?,
            start: UInt32Array::builder(capacity),
            end: UInt32Array::builder(capacity),
            value: Float32Array::builder(capacity),
        })
    }
}

impl BatchBuilder for BigWigBatchBuilder {
    type Record<'a> = BigWigRecord<'a>;

    fn push(&mut self, record: Self::Record<'_>) {
        self.chrom.append_value(record.chrom);
        self.start.append_value(record.start);
        self.end.append_value(record.end);
        self.value.append_value(record.value);
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        RecordBatch::try_from_iter(vec![
            ("chrom", Arc::new(self.chrom.finish()) as ArrayRef),
            ("start", Arc::new(self.start.finish()) as ArrayRef),
            ("end", Arc::new(self.end.finish()) as ArrayRef),
            ("value", Arc::new(self.value.finish()) as ArrayRef),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::ipc::reader::FileReader;
    use arrow::record_batch::RecordBatch;

    fn read_record_batch(region: Option<&str>) -> RecordBatch {
        let mut dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("../fixtures/valid.bigWig");
        let mut reader = BigWigReader::new_from_path(dir.to_str().unwrap()).unwrap();
        let ipc = reader.records_to_ipc(region).unwrap();
        let cursor = std::io::Cursor::new(ipc);
        let mut arrow_reader = FileReader::try_new(cursor, None).unwrap();
        // make sure we have one batch
        assert_eq!(arrow_reader.num_batches(), 1);
        arrow_reader.next().unwrap().unwrap()
    }

    #[test]
    fn test_read_all() {
        let record_batch = read_record_batch(None);
        assert_eq!(record_batch.num_rows(), 100000);
    }

    #[test]
    fn test_region_full() {
        let record_batch = read_record_batch(Some("chr17"));
        assert_eq!(record_batch.num_rows(), 100000);
    }

    #[test]
    fn rest_region_partial() {
        let record_batch = read_record_batch(Some("chr17:59000-60000"));
        assert_eq!(record_batch.num_rows(), 4);
    }
}

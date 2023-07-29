use arrow::array::{
    ArrayRef, UInt32Builder, Float32Builder, UInt32Array, Float32Array,
};
use arrow::ipc::writer::FileWriter;
use arrow::{
    error::ArrowError, record_batch::RecordBatch,
};
use bigtools::{BigWigRead, Value, BBIRead};
use bigtools::utils::reopen::ReopenableFile;
use noodles::core::Region;
use std::sync::Arc;

use crate::batch_builder::{write_ipc, BatchBuilder};

/// A BigWig reader.
pub struct BigWigReader {
    read: BigWigRead<ReopenableFile>,
}

impl BigWigReader {
    /// Creates a BigWig reader.
    pub fn new(path: &str) -> std::io::Result<Self> {
        let read = BigWigRead::open_file(path).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
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
    /// let mut reader = BigWigReader::new("sample.bigWig").unwrap();
    /// let ipc = reader.records_to_ipc(Some("sq0:1-1000")).unwrap();
    /// ```
    pub fn records_to_ipc(&mut self, region: Option<&str>) -> Result<Vec<u8>, ArrowError> {
        let mut batch_builder = BigWigBatchBuilder::new(1024)?;
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
                        let end = self.read.get_chroms().iter().find(|c| c.name == chrom_name).map(|c| c.length);
                        let end = end.ok_or_else(|| ArrowError::InvalidArgumentError("Invalid chromosome".to_string()))?;
                        (start, end)
                    }
                    (None, Some(end)) => {
                        let start = 0;
                        let end = end.get() as u32;
                        (start, end)
                    }
                    (None, None) => {
                        let start = 0;
                        let end = self.read.get_chroms().iter().find(|c| c.name == chrom_name).map(|c| c.length);
                        let end = end.ok_or_else(|| ArrowError::InvalidArgumentError("Invalid chromosome".to_string()))?;
                        (start, end)
                    }
                };
                let values = match self.read.get_interval(&chrom_name, start, end) {
                    Ok(v) => v.map(|v| v.unwrap()),
                    Err(e) => {
                        return Err(ArrowError::ExternalError(Box::new(e)));
                    }
                };
                write_ipc(values, batch_builder)
            }
            None => {
                // Can't use write_ipc, because we have separate iterators for each chrom
                let chroms = self.read.get_chroms().into_iter();
                for chrom in chroms {
                    let start = 0;
                    let end = chrom.length;
                    let values = self.read.get_interval(&chrom.name, start, end).unwrap();
                    values.for_each(|record| batch_builder.push(&record.unwrap()));
                }
                let batch = batch_builder.finish()?;
                let mut writer = FileWriter::try_new(Vec::new(), &batch.schema())?;
                writer.write(&batch)?;
                writer.finish()?;
                writer.into_inner()
            }
        }
    }
}

struct BigWigBatchBuilder {
    start: UInt32Builder,
    end: UInt32Builder,
    value: Float32Builder,
}

impl BigWigBatchBuilder {
    pub fn new(capacity: usize) -> Result<Self, ArrowError> {
        Ok(Self {
            start: UInt32Array::builder(capacity),
            end: UInt32Array::builder(capacity),
            value: Float32Array::builder(capacity),
        })
    }
}

impl<'a> BatchBuilder for BigWigBatchBuilder {
    type Record = Value;

    fn push(&mut self, record: &Self::Record) {
        self.start.append_value(record.start);
        self.end.append_value(record.end);
        self.value.append_value(record.value);
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        RecordBatch::try_from_iter(vec![
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
        let mut reader = BigWigReader::new(dir.to_str().unwrap()).unwrap();
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

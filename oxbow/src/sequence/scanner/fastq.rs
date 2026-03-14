use std::io::{BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;

use crate::sequence::model::batch_builder::BatchBuilder;
use crate::sequence::model::Model;
use crate::sequence::scanner::batch_iterator::BatchIterator;
use crate::util::query::{BgzfChunkReader, ByteRangeReader};
use noodles::bgzf::VirtualPosition;

/// A FASTQ scanner.
///
/// Schema parameters (fields) are declared at construction time. Scan methods
/// accept only column projection, batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::sequence::scanner::fastq::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.R1.fastq").map(BufReader::new).unwrap();
/// let fmt_reader = noodles::fastq::io::Reader::new(inner);
///
/// let scanner = Scanner::new(None).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {
    model: Model,
}

impl Scanner {
    /// Creates a FASTQ scanner from schema parameters.
    ///
    /// `fields`: field names. `None` → `["name", "description", "sequence", "quality"]`.
    pub fn new(fields: Option<Vec<String>>) -> crate::Result<Self> {
        let model = Model::new_fastq(fields)?;
        Ok(Self { model })
    }

    /// Creates a FASTQ scanner from a [`Model`].
    pub fn with_model(model: Model) -> Self {
        Self { model }
    }

    /// Returns a reference to the [`Model`].
    pub fn model(&self) -> &Model {
        &self.model
    }

    /// Returns the field names.
    pub fn field_names(&self) -> Vec<String> {
        self.model.field_names()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self) -> &Schema {
        self.model.schema()
    }

    /// Builds a BatchBuilder applying column projection.
    fn build_batch_builder(
        &self,
        columns: Option<Vec<String>>,
        capacity: usize,
    ) -> crate::Result<BatchBuilder> {
        match columns {
            None => BatchBuilder::from_model(&self.model, capacity),
            Some(cols) => {
                let projected = self.model.project(&cols)?;
                BatchBuilder::from_model(&projected, capacity)
            }
        }
    }
}

impl Scanner {
    /// Returns an iterator yielding record batches.
    pub fn scan<R: BufRead>(
        &self,
        fmt_reader: noodles::fastq::io::Reader<R>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches from specified byte ranges.
    pub fn scan_byte_ranges<R: BufRead + Seek>(
        &self,
        fmt_reader: noodles::fastq::io::Reader<R>,
        byte_ranges: Vec<(u64, u64)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let inner_reader = fmt_reader.into_inner();
        let range_reader = ByteRangeReader::new(inner_reader, byte_ranges);
        let fmt_reader = noodles::fastq::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches from specified virtual position ranges.
    pub fn scan_virtual_ranges<R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek>(
        &self,
        fmt_reader: noodles::fastq::io::Reader<R>,
        vpos_ranges: Vec<(VirtualPosition, VirtualPosition)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let bgzf_reader = fmt_reader.into_inner();
        let range_reader = BgzfChunkReader::new(bgzf_reader, vpos_ranges);
        let fmt_reader = noodles::fastq::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scanner_default() {
        let scanner = Scanner::new(None).unwrap();
        assert_eq!(
            scanner.field_names(),
            vec!["name", "description", "sequence", "quality"]
        );
    }

    #[test]
    fn test_scanner_schema() {
        let scanner = Scanner::new(None).unwrap();
        assert_eq!(scanner.schema().fields().len(), 4);
        let scanner = Scanner::new(Some(vec!["name".to_string(), "quality".to_string()])).unwrap();
        assert_eq!(scanner.schema().fields().len(), 2);
    }

    #[test]
    fn test_scanner_scan() {
        let data =
            b"@SEQ_ID\nGATTA\n+\n!!!!!\n@SEQ_ID2\nCATTAG\n+\n!!!!!!\n@SEQ_ID3\nTTAGGA\n+\n!!!!!!\n";
        let file = std::io::Cursor::new(data);
        let fmt_reader = noodles::fastq::io::Reader::new(file);

        let scanner = Scanner::new(None).unwrap();
        let mut batch_iter = scanner.scan(fmt_reader, None, Some(2), None).unwrap();

        let batch = batch_iter.next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);
        let batch = batch_iter.next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 1);
        assert!(batch_iter.next().is_none());
    }

    #[test]
    fn test_scanner_scan_with_limit() {
        let data =
            b"@SEQ_ID\nGATTA\n+\n!!!!!\n@SEQ_ID2\nCATTAG\n+\n!!!!!!\n@SEQ_ID3\nTTAGGA\n+\n!!!!!!\n";
        let file = std::io::Cursor::new(data);
        let fmt_reader = noodles::fastq::io::Reader::new(file);

        let scanner = Scanner::new(None).unwrap();
        let mut batch_iter = scanner.scan(fmt_reader, None, Some(3), Some(2)).unwrap();

        let batch = batch_iter.next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert!(batch_iter.next().is_none());
    }
}

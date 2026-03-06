use std::io::{BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::{Schema, SchemaRef};

use crate::batch::RecordBatchBuilder as _;
use crate::sequence::model::batch_builder::BatchBuilder;
use crate::sequence::model::field::FASTQ_DEFAULT_FIELD_NAMES;
use crate::sequence::scanner::batch_iterator::BatchIterator;
use crate::util::query::{BgzfChunkReader, ByteRangeReader};
use crate::OxbowError;
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
    fields: Option<Vec<String>>,
    schema: SchemaRef,
}

impl Scanner {
    /// Creates a FASTQ scanner from schema parameters.
    ///
    /// The schema is validated and cached at construction time.
    pub fn new(fields: Option<Vec<String>>) -> crate::Result<Self> {
        let batch_builder = BatchBuilder::new_fastq(fields.clone(), 0)?;
        let schema = batch_builder.schema();
        Ok(Self { fields, schema })
    }

    /// Returns the FASTQ field names.
    pub fn field_names(&self) -> Vec<String> {
        FASTQ_DEFAULT_FIELD_NAMES
            .iter()
            .map(|&s| s.to_string())
            .collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self) -> &Schema {
        &self.schema
    }

    /// Builds a BatchBuilder applying column projection.
    ///
    /// Returns an error if any requested column is not in the declared schema.
    fn build_batch_builder(
        &self,
        columns: Option<Vec<String>>,
        capacity: usize,
    ) -> crate::Result<BatchBuilder> {
        match columns {
            None => BatchBuilder::new_fastq(self.fields.clone(), capacity),
            Some(cols) => {
                let schema_names: Vec<&str> = self
                    .schema
                    .fields()
                    .iter()
                    .map(|f| f.name().as_str())
                    .collect();

                let unknown: Vec<&str> = cols
                    .iter()
                    .filter(|c| !schema_names.iter().any(|s| s.eq_ignore_ascii_case(c)))
                    .map(|c| c.as_str())
                    .collect();
                if !unknown.is_empty() {
                    return Err(OxbowError::invalid_input(format!(
                        "Unknown columns: {:?}. Available columns: {:?}",
                        unknown, schema_names
                    )));
                }

                let declared_field_names: Vec<String> = self.fields.clone().unwrap_or_else(|| {
                    FASTQ_DEFAULT_FIELD_NAMES
                        .iter()
                        .map(|&s| s.to_string())
                        .collect()
                });
                let projected_fields: Vec<String> = declared_field_names
                    .into_iter()
                    .filter(|name| cols.iter().any(|c| c.eq_ignore_ascii_case(name)))
                    .collect();

                BatchBuilder::new_fastq(Some(projected_fields), capacity)
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
            FASTQ_DEFAULT_FIELD_NAMES
                .iter()
                .map(|&s| s.to_string())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_scanner_schema() {
        let scanner = Scanner::new(None).unwrap();
        assert_eq!(
            scanner.schema().fields().len(),
            FASTQ_DEFAULT_FIELD_NAMES.len()
        );
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

use std::io::{self, BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use noodles::core::Region;

use crate::sequence::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::sequence::model::batch_builder::BatchBuilder;
use crate::sequence::model::field::FASTA_DEFAULT_FIELD_NAMES;

/// A FASTA scanner.
///
/// # Examples
///
/// ```no_run
/// use oxbow::sequence::format::fasta::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
/// use noodles::core::Region;
///
/// let inner = File::open("sample.fa").map(BufReader::new).unwrap();
/// let fmt_reader = noodles::fasta::io::Reader::new(inner);
/// let index = noodles::fasta::fai::read("sample.fa.fai").unwrap();
///
/// let scanner = Scanner::default();
/// let regions = vec!["chr1:1-1000", "chr1:1001-2000", "chr1:2001-3000", "chr1:3001-4000"];
/// let regions: Vec<Region> = regions.iter().map(|s| s.parse().unwrap()).collect();
/// let batches = scanner.scan_query(fmt_reader, regions, index, None, Some(2));
/// ```
pub struct Scanner {}

impl Default for Scanner {
    fn default() -> Self {
        Self::new()
    }
}

impl Scanner {
    // Creates a FASTA scanner.
    pub fn new() -> Self {
        Self {}
    }

    /// Returns the FASTA field names.
    pub fn field_names(&self) -> Vec<String> {
        FASTA_DEFAULT_FIELD_NAMES
            .iter()
            .map(|&s| s.to_string())
            .collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self, fields: Option<Vec<String>>) -> io::Result<Schema> {
        let batch_builder = BatchBuilder::new_fasta(fields, 0)?;
        Ok(batch_builder.get_arrow_schema())
    }
}

impl Scanner {
    /// Returns an iterator yielding record batches.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    ///
    /// # Note
    /// Since reference sequences are often large, the default batch size is
    /// set to 1.
    pub fn scan<R: BufRead>(
        &self,
        fmt_reader: noodles::fasta::io::Reader<R>,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1);
        let batch_builder = BatchBuilder::new_fasta(fields, batch_size)?;
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Fetches sequence slice records from a FASTA file by genomic range.
    ///
    /// To read from a BGZF-compressed FASTA file, use `R`: [`noodles::bgzf::IndexedReader`].
    pub fn scan_query<R: BufRead + Seek>(
        &self,
        fmt_reader: noodles::fasta::io::Reader<R>,
        regions: Vec<Region>,
        index: noodles::fasta::fai::Index,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new_fasta(fields, batch_size)?;
        let batch_iter =
            QueryBatchIterator::new(fmt_reader, index, regions, batch_builder, batch_size);
        Ok(batch_iter)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::fasta::fai;
    use std::io::BufReader;

    #[test]
    fn test_scanner_default() {
        let scanner = Scanner::default();
        assert_eq!(
            scanner.field_names(),
            FASTA_DEFAULT_FIELD_NAMES
                .iter()
                .map(|&s| s.to_string())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_scanner_schema() {
        let scanner = Scanner::new();
        let schema = scanner.schema(None).unwrap();
        assert_eq!(schema.fields().len(), FASTA_DEFAULT_FIELD_NAMES.len());
        let schema = scanner
            .schema(Some(vec!["name".to_string(), "sequence".to_string()]))
            .unwrap();
        assert_eq!(schema.fields().len(), 2);
    }

    #[test]
    fn test_scanner_scan() {
        let data = b">seq1\nACGTACGTACGT\n>seq2\nTGCATGCATGCA\n>seq3\nGATTACAGATTACA\n";
        let file = std::io::Cursor::new(data);
        let reader = BufReader::new(file);
        let fmt_reader = noodles::fasta::io::Reader::new(reader);

        let scanner = Scanner::new();
        let mut batch_iter = scanner.scan(fmt_reader, None, Some(2), Some(10)).unwrap();

        let batch = batch_iter.next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);
        let batch = batch_iter.next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 1);
        assert!(batch_iter.next().is_none());
    }

    #[test]
    fn test_scanner_scan_query() {
        let data = b">seq1\nACGTACGTACGT\n>seq2\nTGCATGCATGCA\n>seq3\nGATTACAGATTACA\n";
        let file = std::io::Cursor::new(data);
        let reader = BufReader::new(file);
        let fmt_reader = noodles::fasta::io::Reader::new(reader);
        let index = fai::Index::from(vec![
            fai::Record::new("seq1", 12, 0, 13, 13),
            fai::Record::new("seq2", 12, 12, 13, 13),
            fai::Record::new("seq3", 12, 24, 13, 13),
        ]);

        let scanner = Scanner::new();
        let regions = vec!["seq1:1-4", "seq2:1-4", "seq3:1-4"];
        let regions: Vec<Region> = regions.iter().map(|s| s.parse().unwrap()).collect();
        let mut batch_iter = scanner
            .scan_query(fmt_reader, regions, index, None, Some(2))
            .unwrap();

        let batch = batch_iter.next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);
        let batch = batch_iter.next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 1);
        assert!(batch_iter.next().is_none());
    }
}

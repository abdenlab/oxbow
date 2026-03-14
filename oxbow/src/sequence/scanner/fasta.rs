use std::io::{BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use noodles::core::Region;

use crate::sequence::model::batch_builder::BatchBuilder;
use crate::sequence::model::Model;
use crate::sequence::scanner::batch_iterator::{BatchIterator, QueryBatchIterator};

/// A FASTA scanner.
///
/// Schema parameters (fields) are declared at construction time. Scan methods
/// accept only column projection, batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::sequence::scanner::fasta::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
/// use noodles::core::Region;
///
/// let inner = File::open("sample.fa").map(BufReader::new).unwrap();
/// let fmt_reader = noodles::fasta::io::Reader::new(inner);
/// let index = noodles::fasta::fai::fs::read("sample.fa.fai").unwrap();
///
/// let scanner = Scanner::new(None).unwrap();
/// let regions = vec!["chr1:1-1000", "chr1:1001-2000"];
/// let regions: Vec<Region> = regions.iter().map(|s| s.parse().unwrap()).collect();
/// let batches = scanner.scan_query(fmt_reader, regions, index, None, Some(2));
/// ```
pub struct Scanner {
    model: Model,
}

impl Scanner {
    /// Creates a FASTA scanner from schema parameters.
    ///
    /// `fields`: field names. `None` → `["name", "description", "sequence"]`.
    pub fn new(fields: Option<Vec<String>>) -> crate::Result<Self> {
        let model = Model::new_fasta(fields)?;
        Ok(Self { model })
    }

    /// Creates a FASTA scanner from a [`Model`].
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
    ///
    /// # Note
    /// Since reference sequences are often large, the default batch size is
    /// set to 1.
    pub fn scan<R: BufRead>(
        &self,
        fmt_reader: noodles::fasta::io::Reader<R>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
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
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
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
        let scanner = Scanner::new(None).unwrap();
        assert_eq!(
            scanner.field_names(),
            vec!["name", "description", "sequence"]
        );
    }

    #[test]
    fn test_scanner_schema() {
        let scanner = Scanner::new(None).unwrap();
        assert_eq!(scanner.schema().fields().len(), 3);
        let scanner = Scanner::new(Some(vec!["name".to_string(), "sequence".to_string()])).unwrap();
        assert_eq!(scanner.schema().fields().len(), 2);
    }

    #[test]
    fn test_scanner_scan() {
        let data = b">seq1\nACGTACGTACGT\n>seq2\nTGCATGCATGCA\n>seq3\nGATTACAGATTACA\n";
        let file = std::io::Cursor::new(data);
        let reader = BufReader::new(file);
        let fmt_reader = noodles::fasta::io::Reader::new(reader);

        let scanner = Scanner::new(None).unwrap();
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

        let scanner = Scanner::new(None).unwrap();
        let regions = ["seq1:1-4", "seq2:1-4", "seq3:1-4"];
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

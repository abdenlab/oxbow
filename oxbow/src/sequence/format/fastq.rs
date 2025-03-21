use std::io::{self, BufRead};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;

use crate::sequence::batch_iterator::BatchIterator;
use crate::sequence::model::batch_builder::BatchBuilder;
use crate::sequence::model::field::FASTQ_DEFAULT_FIELD_NAMES;

/// A FASTQ scanner.
///
/// # Examples
///
/// ```no_run
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.R1.fastq").map(BufReader::new)?;
/// let fmt_reader = noodles::fastq::io::Reader::new(inner);
///
/// let scanner = Scanner::new();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {}

impl Scanner {
    // Creates a FASTQ scanner.
    pub fn new() -> Self {
        Self {}
    }

    /// Returns the FASTQ field names.
    pub fn field_names(&self) -> Vec<String> {
        FASTQ_DEFAULT_FIELD_NAMES
            .iter()
            .map(|&s| s.to_string())
            .collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self, fields: Option<Vec<String>>) -> io::Result<Schema> {
        let batch_builder = BatchBuilder::new_fastq(fields, 0)?;
        Ok(batch_builder.get_arrow_schema())
    }
}

impl Scanner {
    /// Returns an iterator yielding record batches.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn scan<R: BufRead>(
        &self,
        fmt_reader: noodles::fastq::io::Reader<R>,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new_fastq(fields, batch_size)?;
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

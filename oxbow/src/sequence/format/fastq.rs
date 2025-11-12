use std::io::{self, BufRead, Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;

use crate::sequence::batch_iterator::BatchIterator;
use crate::sequence::model::batch_builder::BatchBuilder;
use crate::sequence::model::field::FASTQ_DEFAULT_FIELD_NAMES;
use crate::util::query::{BgzfChunkReader, ByteRangeReader};
use noodles::bgzf::VirtualPosition;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;

/// A FASTQ scanner.
///
/// # Examples
///
/// ```no_run
/// use oxbow::sequence::format::fastq::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.R1.fastq").map(BufReader::new).unwrap();
/// let fmt_reader = noodles::fastq::io::Reader::new(inner);
///
/// let scanner = Scanner::default();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {}

impl Default for Scanner {
    fn default() -> Self {
        Self::new()
    }
}

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

    /// Returns an iterator yielding record batches from specified byte ranges.
    ///
    /// This operation requires a seekable (typically uncompressed) source.
    ///
    /// The scan will traverse the specified byte ranges without filtering by genomic coordinates.
    /// This is useful when you have pre-computed file offsets from a custom index. The byte ranges
    /// must align with record boundaries.
    pub fn scan_byte_ranges<R: BufRead + Seek>(
        &self,
        fmt_reader: noodles::fastq::io::Reader<R>,
        byte_ranges: Vec<(u64, u64)>,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new_fastq(fields, batch_size)?;

        let inner_reader = fmt_reader.into_inner();
        let range_reader = ByteRangeReader::new(inner_reader, byte_ranges);
        let fmt_reader = noodles::fastq::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches from specified virtual position ranges.
    ///
    /// This operation requires a BGZF-compressed source.
    ///
    /// The scan will traverse the specified virtual position ranges without filtering by genomic
    /// coordinates. This is useful when you have pre-computed virtual offsets from a custom index.
    pub fn scan_vpos_ranges<R: Read + Seek>(
        &self,
        fmt_reader: noodles::fastq::io::Reader<noodles::bgzf::Reader<R>>,
        vpos_ranges: Vec<(VirtualPosition, VirtualPosition)>,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new_fastq(fields, batch_size)?;

        // Convert virtual position tuples to Chunks
        let chunks: Vec<Chunk> = vpos_ranges
            .into_iter()
            .map(|(start, end)| Chunk::new(start, end))
            .collect();

        let bgzf_reader = fmt_reader.into_inner();
        let range_reader = BgzfChunkReader::new(bgzf_reader, chunks);
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
        let scanner = Scanner::default();
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
        let scanner = Scanner::new();
        let schema = scanner.schema(None).unwrap();
        assert_eq!(schema.fields().len(), FASTQ_DEFAULT_FIELD_NAMES.len());
        let schema = scanner
            .schema(Some(vec!["name".to_string(), "quality".to_string()]))
            .unwrap();
        assert_eq!(schema.fields().len(), 2);
    }

    #[test]
    fn test_scanner_scan() {
        let data =
            b"@SEQ_ID\nGATTA\n+\n!!!!!\n@SEQ_ID2\nCATTAG\n+\n!!!!!!\n@SEQ_ID3\nTTAGGA\n+\n!!!!!!\n";
        let file = std::io::Cursor::new(data);
        let fmt_reader = noodles::fastq::io::Reader::new(file);

        let scanner = Scanner::new();
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

        let scanner = Scanner::new();
        let mut batch_iter = scanner.scan(fmt_reader, None, Some(3), Some(2)).unwrap();

        let batch = batch_iter.next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert!(batch_iter.next().is_none());
    }
}

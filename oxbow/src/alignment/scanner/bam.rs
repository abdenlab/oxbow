use std::io::{BufRead, Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::{Schema, SchemaRef};
use noodles::bgzf::VirtualPosition;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::csi::BinningIndex;

use crate::alignment::model::field::DEFAULT_FIELD_NAMES;
use crate::alignment::model::tag::TagScanner;
use crate::alignment::model::BatchBuilder;
use crate::alignment::scanner::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::batch::RecordBatchBuilder as _;
use crate::util::query::{BgzfChunkReader, ByteRangeReader};
use crate::OxbowError;

/// A BAM scanner.
///
/// Schema parameters (fields, tag definitions) are declared at construction
/// time. Scan methods accept only column projection, batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::alignment::scanner::bam::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.bam").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::bam::io::Reader::new(inner);
/// let header = fmt_reader.read_header().unwrap();
///
/// let tag_defs = Scanner::tag_defs(&mut fmt_reader, Some(1000)).unwrap();
/// let scanner = Scanner::new(header, None, Some(tag_defs)).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {
    header: noodles::sam::Header,
    fields: Option<Vec<String>>,
    tag_defs: Option<Vec<(String, String)>>,
    schema: SchemaRef,
}

impl Scanner {
    /// Creates a BAM scanner from a SAM header and schema parameters.
    ///
    /// The schema is validated and cached at construction time. Scan methods
    /// will use this schema for projection.
    pub fn new(
        header: noodles::sam::Header,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
    ) -> crate::Result<Self> {
        let batch_builder = BatchBuilder::new(header.clone(), fields.clone(), tag_defs.clone(), 0)?;
        let schema = batch_builder.schema();
        Ok(Self {
            header,
            fields,
            tag_defs,
            schema,
        })
    }

    /// Returns a reference to the SAM header.
    pub fn header(&self) -> &noodles::sam::Header {
        &self.header
    }

    /// Returns the reference sequence names.
    pub fn chrom_names(&self) -> Vec<String> {
        self.header
            .reference_sequences()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect()
    }

    /// Returns the reference sequence names and lengths.
    pub fn chrom_sizes(&self) -> Vec<(String, u32)> {
        self.header
            .reference_sequences()
            .iter()
            .map(|(name, r)| (name.to_string(), r.length().get() as u32))
            .collect()
    }

    /// Returns the fixed field names.
    pub fn field_names(&self) -> Vec<String> {
        DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self) -> &Schema {
        &self.schema
    }

    /// Builds a BatchBuilder applying column projection.
    ///
    /// - `columns: None` → all declared columns
    /// - `columns: Some(cols)` → only the specified top-level columns;
    ///   fixed fields are intersected, "tags" includes all tag_defs if present
    ///
    /// Returns an error if any requested column is not in the declared schema.
    fn build_batch_builder(
        &self,
        columns: Option<Vec<String>>,
        capacity: usize,
    ) -> crate::Result<BatchBuilder> {
        match columns {
            None => BatchBuilder::new(
                self.header.clone(),
                self.fields.clone(),
                self.tag_defs.clone(),
                capacity,
            ),
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

                // Determine which declared field names to keep
                let declared_field_names: Vec<String> = self.fields.clone().unwrap_or_else(|| {
                    DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect()
                });
                let projected_fields: Vec<String> = declared_field_names
                    .into_iter()
                    .filter(|name| cols.iter().any(|c| c.eq_ignore_ascii_case(name)))
                    .collect();

                // Include tags only if "tags" is in the requested columns
                let tag_defs = if cols.iter().any(|c| c == "tags") {
                    self.tag_defs.clone()
                } else {
                    None
                };

                BatchBuilder::new(
                    self.header.clone(),
                    Some(projected_fields),
                    tag_defs,
                    capacity,
                )
            }
        }
    }
}

impl Scanner {
    /// Discovers tag definitions by scanning over records.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn tag_defs<R: Read>(
        fmt_reader: &mut noodles::bam::io::Reader<R>,
        scan_rows: Option<usize>,
    ) -> crate::Result<Vec<(String, String)>> {
        let records = fmt_reader.records();
        let mut tag_scanner = TagScanner::new();
        match scan_rows {
            None => {
                for result in records {
                    if let Ok(record) = result {
                        tag_scanner.push(&record);
                    } else {
                        eprintln!("Failed to read record");
                    }
                }
            }
            Some(n) => {
                for result in records.take(n) {
                    if let Ok(record) = result {
                        tag_scanner.push(&record);
                    } else {
                        eprintln!("Failed to read record");
                    }
                }
            }
        }
        Ok(tag_scanner.collect())
    }

    /// Returns an iterator yielding record batches.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn scan<R: Read>(
        &self,
        fmt_reader: noodles::bam::io::Reader<R>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches satisfying a genomic range query.
    ///
    /// This operation requires a BGZF source and an Index.
    ///
    /// The scan will traverse one or more virtual position ranges and filter
    /// for records that overlap the given region. The cursor will stop at the
    /// end of the last record scanned.
    pub fn scan_query<R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek>(
        &self,
        fmt_reader: noodles::bam::io::Reader<R>,
        region: noodles::core::Region,
        index: impl BinningIndex,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let interval = region.interval();

        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let reference_sequence_id = super::resolve_chrom_id(&self.header, region.name())?;
        let chunks = index.query(reference_sequence_id, interval)?;
        let bgzf_reader = fmt_reader.into_inner();
        let query_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::bam::io::Reader::from(query_reader);
        let batch_iter = QueryBatchIterator::new(
            fmt_reader,
            self.header.clone(),
            reference_sequence_id,
            interval,
            batch_builder,
            batch_size,
            limit,
        );
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches of unaligned reads.
    ///
    /// This operation requires a BGZF source and an Index.
    ///
    /// The scan will start at the offset where unmapped reads begin and
    /// continue until the source stream is exhausted.
    pub fn scan_unmapped<R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek>(
        &self,
        mut fmt_reader: noodles::bam::io::Reader<R>,
        index: impl BinningIndex,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        // This will make the reader seek to the beginning of the unmapped read records.
        let _ = fmt_reader.query_unmapped(&index)?;

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
        fmt_reader: noodles::bam::io::Reader<R>,
        byte_ranges: Vec<(u64, u64)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let inner_reader = fmt_reader.into_inner();
        let range_reader = ByteRangeReader::new(inner_reader, byte_ranges);
        let fmt_reader = noodles::bam::io::Reader::from(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches from specified virtual position ranges.
    ///
    /// This operation requires a BGZF-compressed source.
    ///
    /// The scan will traverse the specified virtual position ranges without filtering by genomic
    /// coordinates. This is useful when you have pre-computed virtual offsets from a custom index.
    pub fn scan_virtual_ranges<R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek>(
        &self,
        fmt_reader: noodles::bam::io::Reader<R>,
        vpos_ranges: Vec<(VirtualPosition, VirtualPosition)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        // Convert virtual position tuples to Chunks
        let chunks: Vec<Chunk> = vpos_ranges
            .into_iter()
            .map(|(start, end)| Chunk::new(start, end))
            .collect();

        let bgzf_reader = fmt_reader.into_inner();
        let range_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::bam::io::Reader::from(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::num::NonZero;

    fn mt_reader() -> (
        noodles::sam::Header,
        noodles::bam::io::Reader<noodles::bgzf::io::MultithreadedReader<std::fs::File>>,
    ) {
        let file = std::fs::File::open("../fixtures/sample.bam").unwrap();
        let mt_reader = noodles::bgzf::io::MultithreadedReader::with_worker_count(
            NonZero::new(2usize).unwrap(),
            file,
        );
        let mut fmt_reader = noodles::bam::io::Reader::from(mt_reader);
        let header = fmt_reader.read_header().unwrap();
        (header, fmt_reader)
    }

    #[test]
    fn test_scan_with_multithreaded_reader() {
        let (header, fmt_reader) = mt_reader();
        let scanner = Scanner::new(header, None, None).unwrap();
        let mut batches = scanner.scan(fmt_reader, None, None, Some(10)).unwrap();

        let batch = batches.next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 10);
        assert!(batch.num_columns() > 0);
    }

    #[test]
    fn test_scan_query_with_multithreaded_reader() {
        let (header, fmt_reader) = mt_reader();
        let scanner = Scanner::new(header, None, None).unwrap();

        let index = noodles::bam::bai::fs::read("../fixtures/sample.bam.bai").unwrap();

        let region = "chr1:1-100000".parse().unwrap();
        let mut batches = scanner
            .scan_query(fmt_reader, region, index, None, None, Some(10))
            .unwrap();

        let batch = batches.next().unwrap().unwrap();
        assert!(batch.num_rows() > 0);
        assert!(batch.num_columns() > 0);
    }
}

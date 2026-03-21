use std::io::{BufRead, Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use noodles::bgzf::VirtualPosition;
use noodles::csi::BinningIndex;

use crate::alignment::model::tag::TagScanner;
use crate::alignment::model::BatchBuilder;
use crate::alignment::scanner::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::alignment::AlignmentModel;
use crate::util::query::{BgzfChunkReader, ByteRangeReader};
use crate::{CoordSystem, Region, Select};

/// A BAM scanner.
///
/// Schema parameters (fields, tag definitions) are declared at construction
/// time. Scan methods accept only column projection, batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::alignment::scanner::bam::Scanner;
/// use oxbow::Select;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.bam").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::bam::io::Reader::new(inner);
/// let header = fmt_reader.read_header().unwrap();
///
/// let tag_defs = Scanner::tag_defs(&mut fmt_reader, Some(1000)).unwrap();
/// use oxbow::CoordSystem;
/// let scanner = Scanner::new(header, Select::All, Some(tag_defs), CoordSystem::OneClosed).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {
    header: noodles::sam::Header,
    model: AlignmentModel,
}

impl Scanner {
    /// Creates a BAM scanner from a SAM header and schema parameters.
    ///
    /// - `fields`: standard SAM field selection.
    /// - `tag_defs`: `None` → no tags column. `Some(vec![])` → empty struct.
    /// - `coord_system`: output coordinate system for position columns.
    pub fn new(
        header: noodles::sam::Header,
        fields: Select<String>,
        tag_defs: Option<Vec<(String, String)>>,
        coord_system: CoordSystem,
    ) -> crate::Result<Self> {
        let model = AlignmentModel::new(fields, tag_defs, coord_system)?;
        Ok(Self { header, model })
    }

    /// Creates a BAM scanner from an [`AlignmentModel`].
    pub fn with_model(header: noodles::sam::Header, model: AlignmentModel) -> Self {
        Self { header, model }
    }

    /// Returns a reference to the [`AlignmentModel`].
    pub fn model(&self) -> &AlignmentModel {
        &self.model
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

    /// Returns the field names declared in the model.
    pub fn field_names(&self) -> Vec<String> {
        self.model.field_names()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self) -> &Schema {
        self.model.schema()
    }

    /// Builds a BatchBuilder applying column projection.
    ///
    /// Returns an error if any requested column is not in the declared schema.
    fn build_batch_builder(
        &self,
        columns: Option<Vec<String>>,
        capacity: usize,
    ) -> crate::Result<BatchBuilder> {
        let model = match columns {
            None => self.model.clone(),
            Some(cols) => self.model.project(&cols)?,
        };
        BatchBuilder::from_model(&model, self.header.clone(), capacity)
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
        region: Region,
        index: impl BinningIndex,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let region = region.to_noodles()?;
        let interval = region.interval();

        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let reference_sequence_id = super::resolve_chrom_id(&self.header, region.name())?;
        let chunks = index.query(reference_sequence_id, interval)?;
        let chunks = chunks.into_iter().map(|c| (c.start(), c.end())).collect();
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

        let bgzf_reader = fmt_reader.into_inner();
        let range_reader = BgzfChunkReader::new(bgzf_reader, vpos_ranges);
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
        let scanner = Scanner::new(header, Select::All, None, CoordSystem::OneClosed).unwrap();
        let mut batches = scanner.scan(fmt_reader, None, None, Some(10)).unwrap();

        let batch = batches.next().unwrap().unwrap();
        assert_eq!(batch.num_rows(), 10);
        assert!(batch.num_columns() > 0);
    }

    #[test]
    fn test_scan_query_with_multithreaded_reader() {
        let (header, fmt_reader) = mt_reader();
        let scanner = Scanner::new(header, Select::All, None, CoordSystem::OneClosed).unwrap();

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

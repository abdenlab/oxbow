use std::io::{self, Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use noodles::csi::BinningIndex;

use crate::alignment::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::alignment::model::field::DEFAULT_FIELD_NAMES;
use crate::alignment::model::tag::TagScanner;
use crate::alignment::model::BatchBuilder;
use crate::util::query::BgzfChunkReader;

/// A BAM scanner.
///
/// # Examples
///
/// ```no_run
/// use oxbow::alignment::format::bam::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.bam").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::bam::io::Reader::new(inner);
/// let header = fmt_reader.read_header().unwrap();
///
/// let scanner = Scanner::new(header);
/// let tag_defs = scanner.tag_defs(&mut fmt_reader, Some(1000)).unwrap();
/// let batches = scanner.scan(fmt_reader, None, Some(tag_defs), None, Some(1000));
/// ```
pub struct Scanner {
    header: noodles::sam::Header,
}

impl Scanner {
    /// Creates a BAM scanner from a SAM header.
    pub fn new(header: noodles::sam::Header) -> Self {
        Self { header }
    }

    /// Returns the SAM header.
    pub fn header(&self) -> noodles::sam::Header {
        self.header.clone()
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

    /// Returns the standard field names.
    pub fn field_names(&self) -> Vec<String> {
        DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(
        &self,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
    ) -> io::Result<Schema> {
        let header = self.header();
        let batch_builder = BatchBuilder::new(header, fields, tag_defs, 0)?;
        Ok(batch_builder.get_arrow_schema())
    }
}

impl Scanner {
    /// Discovers tag definitions by scanning over records.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn tag_defs<R: Read>(
        &self,
        fmt_reader: &mut noodles::bam::io::Reader<R>,
        scan_rows: Option<usize>,
    ) -> io::Result<Vec<(String, String)>> {
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

    /// Returns an iterator yielding batches of records.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn scan<R: Read>(
        &self,
        fmt_reader: noodles::bam::io::Reader<R>,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let header = self.header();
        let batch_builder = BatchBuilder::new(header, fields, tag_defs, batch_size)?;
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
    #[allow(clippy::too_many_arguments)]
    pub fn scan_query<R: Read + Seek>(
        &self,
        fmt_reader: noodles::bam::io::Reader<noodles::bgzf::Reader<R>>,
        region: noodles::core::Region,
        index: impl BinningIndex,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let interval = region.interval();

        let batch_builder = BatchBuilder::new(self.header(), fields, tag_defs, batch_size)?;

        let reference_sequence_id = resolve_chrom_id(&self.header, region.name())?;
        let chunks = index.query(reference_sequence_id, interval).unwrap();
        let bgzf_reader = fmt_reader.into_inner();
        let query_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::bam::io::Reader::from(query_reader);
        let batch_iter = QueryBatchIterator::new(
            fmt_reader,
            self.header(),
            reference_sequence_id,
            interval,
            batch_builder,
            batch_size,
            limit,
        );
        Ok(batch_iter)
    }

    /// Returns an iterator yielding batches of records corresponding to unaligned reads.
    ///
    /// This operation requires a BGZF source and an Index.
    ///
    /// The scan will start at the offset where unmapped reads begin and
    /// continue until the source stream is exhausted.
    pub fn scan_unmapped<R: Read + Seek>(
        &self,
        mut fmt_reader: noodles::bam::io::Reader<noodles::bgzf::Reader<R>>,
        index: impl BinningIndex,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let header = self.header();
        let batch_builder = BatchBuilder::new(header, fields, tag_defs, batch_size)?;

        // This will make the reader seek to the beginning of the unmapped read records.
        let _ = fmt_reader.query_unmapped(&index)?;

        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

fn resolve_chrom_id(
    header: &noodles::sam::Header,
    reference_sequence_name: &[u8],
) -> io::Result<usize> {
    let Some(id) = header
        .reference_sequences()
        .get_index_of(reference_sequence_name)
    else {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "Reference sequence {:?} not found in index header.",
                reference_sequence_name
            ),
        ));
    };
    Ok(id)
}

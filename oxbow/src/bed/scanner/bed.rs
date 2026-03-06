use std::io::{self, BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::{Schema, SchemaRef};
use noodles::bgzf::VirtualPosition;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::csi::BinningIndex;

use crate::batch::RecordBatchBuilder as _;
use crate::bed::model::BatchBuilder;
use crate::bed::model::BedSchema;
use crate::bed::scanner::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::util::query::{BgzfChunkReader, ByteRangeReader};

/// A BED scanner.
///
/// Schema parameters (fields) are declared at construction time. Scan methods
/// accept only column projection, batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::bed::scanner::bed::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.bed").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::bed::io::Reader::new(inner);
///
/// let bed_schema = "bed6+3".parse().unwrap();
/// let scanner = Scanner::new(bed_schema, None).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000)).unwrap();
/// ```
pub struct Scanner {
    bed_schema: BedSchema,
    fields: Option<Vec<String>>,
    schema: SchemaRef,
}

impl Scanner {
    /// Creates a BED scanner from a BED schema and optional field names.
    ///
    /// The schema is validated and cached at construction time.
    pub fn new(bed_schema: BedSchema, fields: Option<Vec<String>>) -> io::Result<Self> {
        let batch_builder = BatchBuilder::new(fields.clone(), &bed_schema, 0)?;
        let schema = batch_builder.schema();
        Ok(Self {
            bed_schema,
            fields,
            schema,
        })
    }

    /// Returns the Arrow schema.
    pub fn schema(&self) -> &Schema {
        &self.schema
    }

    /// Returns the BED field names from the schema.
    pub fn field_names(&self) -> Vec<String> {
        self.bed_schema.field_names()
    }

    /// Builds a BatchBuilder applying column projection.
    ///
    /// Returns an error if any requested column is not in the declared schema.
    fn build_batch_builder(
        &self,
        columns: Option<Vec<String>>,
        capacity: usize,
    ) -> io::Result<BatchBuilder> {
        match columns {
            None => BatchBuilder::new(self.fields.clone(), &self.bed_schema, capacity),
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
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Unknown columns: {:?}. Available columns: {:?}",
                            unknown, schema_names
                        ),
                    ));
                }

                BatchBuilder::new(Some(cols), &self.bed_schema, capacity)
            }
        }
    }
}

impl Scanner {
    /// Returns an iterator yielding record batches.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn scan<R: BufRead>(
        &self,
        fmt_reader: noodles::bed::io::Reader<3, R>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches satisfying a genomic range query.
    ///
    /// This operation requires a BGZF source and an Index.
    ///
    /// The scan will consume contiguous "chunks" of BGZF blocks and filter for
    /// records that overlap the given region. The cursor will move to the end
    /// of the last record scanned.
    pub fn scan_query<R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek>(
        &self,
        fmt_reader: noodles::bed::io::Reader<3, R>,
        region: noodles::core::Region,
        index: impl BinningIndex,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let reference_sequence_name = region.name().to_string();
        let interval = region.interval();

        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let Some(header) = index.header() else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "Index header not found.",
            ));
        };
        let reference_sequence_id = resolve_chrom_id(header, &reference_sequence_name)?;
        let chunks = index.query(reference_sequence_id, interval)?;
        let bgzf_reader = fmt_reader.into_inner();
        let query_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::bed::io::Reader::new(query_reader);
        let batch_iter = QueryBatchIterator::new(
            fmt_reader,
            header.clone(),
            reference_sequence_id,
            interval,
            batch_builder,
            batch_size,
            limit,
        );
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
        fmt_reader: noodles::bed::io::Reader<3, R>,
        byte_ranges: Vec<(u64, u64)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let inner_reader = fmt_reader.into_inner();
        let range_reader = ByteRangeReader::new(inner_reader, byte_ranges);
        let fmt_reader = noodles::bed::io::Reader::new(range_reader);
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
        fmt_reader: noodles::bed::io::Reader<3, R>,
        vpos_ranges: Vec<(VirtualPosition, VirtualPosition)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        // Convert virtual position tuples to Chunks
        let chunks: Vec<Chunk> = vpos_ranges
            .into_iter()
            .map(|(start, end)| Chunk::new(start, end))
            .collect();

        let bgzf_reader = fmt_reader.into_inner();
        let range_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::bed::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

fn resolve_chrom_id(
    header: &noodles::csi::binning_index::index::Header,
    reference_sequence_name: &str,
) -> io::Result<usize> {
    let Some(id) = header
        .reference_sequence_names()
        .get_index_of(reference_sequence_name.as_bytes())
    else {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "Reference sequence {} not found in index header.",
                reference_sequence_name
            ),
        ));
    };
    Ok(id)
}

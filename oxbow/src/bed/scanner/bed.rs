use std::io::{BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use noodles::bgzf::VirtualPosition;
use noodles::csi::BinningIndex;

use crate::bed::model::BatchBuilder;
use crate::bed::model::BedSchema;
use crate::bed::model::Model;
use crate::bed::scanner::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::util::query::{BgzfChunkReader, ByteRangeReader};
use crate::OxbowError;

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
    model: Model,
}

impl Scanner {
    /// Creates a BED scanner from a BED schema and optional field projection.
    ///
    /// - `bed_schema`: the parsing interpretation.
    /// - `fields`: column names to project. `None` → all fields from the schema.
    pub fn new(bed_schema: BedSchema, fields: Option<Vec<String>>) -> crate::Result<Self> {
        let model = Model::new(bed_schema, fields)?;
        Ok(Self { model })
    }

    /// Creates a BED scanner from a [`Model`].
    pub fn with_model(model: Model) -> Self {
        Self { model }
    }

    /// Returns a reference to the [`Model`].
    pub fn model(&self) -> &Model {
        &self.model
    }

    /// Returns the Arrow schema.
    pub fn schema(&self) -> &Schema {
        self.model.schema()
    }

    /// Returns the field names.
    pub fn field_names(&self) -> Vec<String> {
        self.model.field_names()
    }

    /// Builds a BatchBuilder applying column projection.
    fn build_batch_builder(
        &self,
        columns: Option<Vec<String>>,
        capacity: usize,
    ) -> crate::Result<BatchBuilder> {
        match columns {
            None => BatchBuilder::new(self.model.bed_schema(), None, capacity),
            Some(cols) => {
                let projected = self.model.project(&cols)?;
                BatchBuilder::new(
                    projected.bed_schema(),
                    Some(projected.field_names()),
                    capacity,
                )
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
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let reference_sequence_name = region.name().to_string();
        let interval = region.interval();

        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let Some(header) = index.header() else {
            return Err(OxbowError::not_found("Index header not found."));
        };
        let reference_sequence_id = resolve_chrom_id(header, &reference_sequence_name)?;
        let chunks = index.query(reference_sequence_id, interval)?;
        let chunks = chunks.into_iter().map(|c| (c.start(), c.end())).collect();
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
    ) -> crate::Result<impl RecordBatchReader> {
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
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let bgzf_reader = fmt_reader.into_inner();
        let range_reader = BgzfChunkReader::new(bgzf_reader, vpos_ranges);
        let fmt_reader = noodles::bed::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

fn resolve_chrom_id(
    header: &noodles::csi::binning_index::index::Header,
    reference_sequence_name: &str,
) -> crate::Result<usize> {
    let Some(id) = header
        .reference_sequence_names()
        .get_index_of(reference_sequence_name.as_bytes())
    else {
        return Err(OxbowError::not_found(format!(
            "Reference sequence {} not found in index header.",
            reference_sequence_name
        )));
    };
    Ok(id)
}

use std::io::{BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use noodles::bgzf::VirtualPosition;
use noodles::csi::binning_index;
use noodles::csi::BinningIndex;

use crate::gxf::model::attribute::AttributeScanner;
use crate::gxf::model::attribute::Push as _;
use crate::gxf::model::BatchBuilder;
use crate::gxf::model::Model;
use crate::gxf::scanner::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::util::query::{BgzfChunkReader, ByteRangeReader};
use crate::{CoordSystem, OxbowError, Region, Select};

/// A GTF scanner.
///
/// Schema parameters (fields, attribute definitions) are declared at
/// construction time. Scan methods accept only column projection,
/// batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::gxf::scanner::gtf::Scanner;
/// use oxbow::Select;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.gtf").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::gtf::io::Reader::new(inner);
///
/// let attr_defs = Scanner::attribute_defs(&mut fmt_reader, Some(1000)).unwrap();
/// use oxbow::CoordSystem;
/// let scanner = Scanner::new(None, Select::All, Some(attr_defs), CoordSystem::OneClosed).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {
    header: Option<binning_index::index::Header>,
    model: Model,
}

impl Scanner {
    /// Creates a GTF scanner from schema parameters.
    ///
    /// - `fields`: standard GXF field selection. `All` → all 8 standard fields.
    /// - `attr_defs`: `None` → no attributes column. `Some(vec![])` → empty struct.
    /// - `coord_system`: output coordinate system for position columns.
    pub fn new(
        header: Option<binning_index::index::Header>,
        fields: Select<String>,
        attr_defs: Option<Vec<(String, String)>>,
        coord_system: CoordSystem,
    ) -> crate::Result<Self> {
        let model = Model::new(fields, attr_defs, coord_system)?;
        Ok(Self { header, model })
    }

    /// Creates a GTF scanner from a [`Model`].
    pub fn with_model(header: Option<binning_index::index::Header>, model: Model) -> Self {
        Self { header, model }
    }

    /// Returns a reference to the [`Model`].
    pub fn model(&self) -> &Model {
        &self.model
    }

    /// Returns the index header if one was provided.
    pub fn header(&self) -> Option<&binning_index::index::Header> {
        self.header.as_ref()
    }

    /// Returns the reference sequence names if an index header was provided.
    pub fn chrom_names(&self) -> crate::Result<Vec<String>> {
        if let Some(header) = &self.header {
            Ok(header
                .reference_sequence_names()
                .iter()
                .map(|name| name.to_string())
                .collect())
        } else {
            Err(OxbowError::not_found("Index header not found."))
        }
    }

    /// Returns the fixed field names.
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
    /// Discovers attribute definitions by scanning over records.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    pub fn attribute_defs<R: BufRead>(
        fmt_reader: &mut noodles::gtf::io::Reader<R>,
        scan_rows: Option<usize>,
    ) -> crate::Result<Vec<(String, String)>> {
        let lines = fmt_reader.lines();
        let mut attr_scanner = AttributeScanner::new();
        match scan_rows {
            None => {
                for line in lines {
                    match line {
                        Ok(line) => match line.as_record() {
                            Some(result) => attr_scanner.push(result?),
                            None => continue,
                        },
                        Err(e) => eprintln!("Failed to read line: {}", e),
                    }
                }
            }
            Some(n) => {
                for line in lines.take(n) {
                    match line {
                        Ok(line) => match line.as_record() {
                            Some(result) => attr_scanner.push(result?),
                            None => continue,
                        },
                        Err(e) => eprintln!("Failed to read line: {}", e),
                    }
                }
            }
        }
        Ok(attr_scanner.collect())
    }

    /// Returns an iterator yielding record batches.
    pub fn scan<R: BufRead>(
        &self,
        fmt_reader: noodles::gtf::io::Reader<R>,
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
    pub fn scan_query<R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek>(
        &self,
        fmt_reader: noodles::gtf::io::Reader<R>,
        region: Region,
        index: impl BinningIndex,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let region = region.to_noodles()?;
        let reference_sequence_name = region.name().to_string();
        let interval = region.interval();

        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let Some(header) = index.header() else {
            return Err(OxbowError::not_found("Index header not found."));
        };
        let reference_sequence_id = super::resolve_chrom_id(header, &reference_sequence_name)?;
        let chunks = index.query(reference_sequence_id, interval)?;
        let chunks = chunks.into_iter().map(|c| (c.start(), c.end())).collect();
        let bgzf_reader = fmt_reader.into_inner();
        let query_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::gff::io::Reader::new(query_reader);
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
    pub fn scan_byte_ranges<R: BufRead + Seek>(
        &self,
        fmt_reader: noodles::gtf::io::Reader<R>,
        byte_ranges: Vec<(u64, u64)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let inner_reader = fmt_reader.into_inner();
        let range_reader = ByteRangeReader::new(inner_reader, byte_ranges);
        let fmt_reader = noodles::gtf::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches from specified virtual position ranges.
    pub fn scan_virtual_ranges<R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek>(
        &self,
        fmt_reader: noodles::gtf::io::Reader<R>,
        vpos_ranges: Vec<(VirtualPosition, VirtualPosition)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let bgzf_reader = fmt_reader.into_inner();
        let range_reader = BgzfChunkReader::new(bgzf_reader, vpos_ranges);
        let fmt_reader = noodles::gtf::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

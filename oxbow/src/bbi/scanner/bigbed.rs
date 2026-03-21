use std::io::{Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema as ArrowSchema;
use bigtools::BigBedRead;

use crate::bbi::model::base::BatchBuilder;
use crate::bbi::model::base::BedSchema;
use crate::bbi::model::base::Model;
use crate::bbi::scanner::batch_iterator::base::{BigBedBatchIterator, BigBedQueryBatchIterator};
use crate::{CoordSystem, Region, Select};

/// A BigBed scanner.
///
/// Schema parameters (fields) are declared at construction time. Scan methods
/// accept only column projection, batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::bbi::scanner::bigbed::Scanner;
///
/// let mut fmt_reader = bigtools::BigBedRead::open_file("sample.bigBed").unwrap();
/// let info = fmt_reader.info();
///
/// use oxbow::Select;
/// use oxbow::CoordSystem;
/// let scanner = Scanner::new("bed12".parse().unwrap(), info.clone(), Select::All, CoordSystem::ZeroHalfOpen).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
pub struct Scanner {
    model: Model,
    info: bigtools::BBIFileInfo,
}

impl Scanner {
    /// Creates a BigBed scanner from a BED schema, BBI file info, and optional field names.
    ///
    /// - `bed_schema`: the parsing interpretation.
    /// - `info`: the BBI file info.
    /// - `fields`: column names to project.
    /// - `coord_system`: output coordinate system for position columns.
    pub fn new(
        bed_schema: BedSchema,
        info: bigtools::BBIFileInfo,
        fields: Select<String>,
        coord_system: CoordSystem,
    ) -> crate::Result<Self> {
        let model = Model::new(bed_schema, fields, coord_system)?;
        Ok(Self { model, info })
    }

    /// Creates a BigBed scanner from a [`Model`].
    pub fn with_model(model: Model, info: bigtools::BBIFileInfo) -> Self {
        Self { model, info }
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
    pub fn schema(&self) -> &ArrowSchema {
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
    /// Returns the BBI file info.
    pub fn info(&self) -> &bigtools::BBIFileInfo {
        &self.info
    }

    /// Returns the reference sequence names.
    pub fn chrom_names(&self) -> Vec<String> {
        let chroms = &self.info.chrom_info;
        chroms.iter().map(|info| info.name.clone()).collect()
    }

    /// Returns the reference sequence names and lengths.
    pub fn chrom_sizes(&self) -> Vec<(String, u32)> {
        let chroms = &self.info.chrom_info;
        chroms
            .iter()
            .map(|info| (info.name.clone(), info.length))
            .collect()
    }

    /// Returns the zoom/reduction level resolutions in the BigBed file.
    pub fn zoom_levels(&self) -> Vec<u32> {
        let zoom_levels: Vec<u32> = self
            .info
            .zoom_headers
            .iter()
            .map(|header| header.reduction_level)
            .collect();
        zoom_levels
    }

    /// Returns an iterator yielding record batches.
    pub fn scan<R: Read + Seek>(
        &self,
        fmt_reader: BigBedRead<R>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        let batch_iter = BigBedBatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches satisfying a genomic range query.
    pub fn scan_query<R: Read + Seek + Send + 'static>(
        &self,
        fmt_reader: BigBedRead<R>,
        region: Region,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let region = region.to_noodles()?;
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        let batch_iter =
            BigBedQueryBatchIterator::new(fmt_reader, region, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

use std::io::{Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema as ArrowSchema;
use bigtools::BigWigRead;

use crate::bbi::model::base::BatchBuilder;
use crate::bbi::model::base::BedSchema;
use crate::bbi::scanner::batch_iterator::base::{BigWigBatchIterator, BigWigQueryBatchIterator};
use crate::bed::model::Model;

/// A BigWig scanner.
///
/// Schema parameters (fields) are declared at construction time. Scan methods
/// accept only column projection, batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::bbi::scanner::bigwig::Scanner;
/// let mut fmt_reader = bigtools::BigWigRead::open_file("sample.bigWig").unwrap();
/// let info = fmt_reader.info();
///
/// let scanner = Scanner::new(info.clone(), None).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {
    model: Model,
    info: bigtools::BBIFileInfo,
}

impl Scanner {
    /// Creates a BigWig scanner from BBI file info and optional field names.
    pub fn new(info: bigtools::BBIFileInfo, fields: Option<Vec<String>>) -> crate::Result<Self> {
        let bed_schema: BedSchema = "bedGraph".parse().unwrap();
        let model = Model::new(bed_schema, fields)?;
        Ok(Self { model, info })
    }

    /// Creates a BigWig scanner from a [`Model`].
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
            None => BatchBuilder::new(
                self.model.bed_schema().clone(),
                Some(self.model.field_names()),
                capacity,
            ),
            Some(cols) => {
                let projected = self.model.project(&cols)?;
                BatchBuilder::new(
                    projected.bed_schema().clone(),
                    Some(projected.field_names()),
                    capacity,
                )
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

    /// Returns the zoom/reduction level resolutions in the BigWig file.
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
        fmt_reader: BigWigRead<R>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        let batch_iter = BigWigBatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches satisfying a genomic range query.
    pub fn scan_query<R: Read + Seek + Send + 'static>(
        &self,
        fmt_reader: BigWigRead<R>,
        region: noodles::core::Region,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        let batch_iter =
            BigWigQueryBatchIterator::new(fmt_reader, region, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

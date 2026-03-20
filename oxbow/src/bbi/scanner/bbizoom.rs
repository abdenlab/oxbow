use std::io::{Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema as ArrowSchema;

pub use super::BBIReader;
use crate::bbi::model::zoom::BatchBuilder;
use crate::bbi::model::zoom::Model;
use crate::bbi::scanner::batch_iterator::zoom::{BBIZoomBatchIterator, BBIZoomQueryBatchIterator};
use crate::{CoordSystem, Select};

/// A scanner for the summary statistics from BBI file zoom level.
///
/// Schema parameters (fields) are declared at construction time. Scan methods
/// accept only column projection, batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::bbi::scanner::bbizoom::Scanner;
/// use oxbow::bbi::BBIReader;
///
/// let mut fmt_reader = bigtools::BigWigRead::open_file("sample.bigWig").unwrap();
/// let info = fmt_reader.info();
/// let ref_names = info.chrom_info.iter().map(|c| c.name.clone()).collect();
/// let zoom_levels: Vec<u32> = info.zoom_headers.iter().map(|h| h.reduction_level).collect();
/// use oxbow::{CoordSystem, Select};
/// let scanner = Scanner::new(ref_names, zoom_levels[0], Select::All, CoordSystem::ZeroHalfOpen).unwrap();
/// let batches = scanner.scan(BBIReader::BigWig(fmt_reader), None, None, Some(1000));
pub struct Scanner {
    ref_names: Vec<String>,
    zoom_level: u32,
    model: Model,
}

impl Scanner {
    /// Creates a BBI zoom level scanner.
    pub fn new(
        ref_names: Vec<String>,
        zoom_level: u32,
        fields: Select<String>,
        coord_system: CoordSystem,
    ) -> crate::Result<Self> {
        let model = Model::new(fields, coord_system)?;
        Ok(Self {
            ref_names,
            zoom_level,
            model,
        })
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
        let cs = self.model.coord_system();
        match columns {
            None => BatchBuilder::new(
                &self.ref_names,
                Some(self.model.field_names()),
                cs,
                capacity,
            ),
            Some(cols) => {
                let projected = self.model.project(&cols)?;
                BatchBuilder::new(&self.ref_names, Some(projected.field_names()), cs, capacity)
            }
        }
    }
}

impl Scanner {
    /// Returns an iterator over record batches.
    pub fn scan<R: Read + Seek + Send + 'static>(
        &self,
        reader: BBIReader<R>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<Box<dyn RecordBatchReader + Send>> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        match reader {
            BBIReader::BigWig(reader) => {
                let batch_iter = BBIZoomBatchIterator::from_bigwig(
                    reader,
                    self.zoom_level,
                    batch_builder,
                    batch_size,
                    limit,
                );
                Ok(Box::new(batch_iter))
            }
            BBIReader::BigBed(reader) => {
                let batch_iter = BBIZoomBatchIterator::from_bigbed(
                    reader,
                    self.zoom_level,
                    batch_builder,
                    batch_size,
                    limit,
                );
                Ok(Box::new(batch_iter))
            }
        }
    }

    /// Returns an iterator over record batches satisfying a genomic range query.
    pub fn scan_query<R: Read + Seek + Send + 'static>(
        &self,
        reader: BBIReader<R>,
        region: noodles::core::Region,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        match reader {
            BBIReader::BigWig(reader) => {
                let batch_iter = BBIZoomQueryBatchIterator::from_bigwig(
                    reader,
                    region,
                    self.zoom_level,
                    batch_builder,
                    batch_size,
                    limit,
                );
                Ok(batch_iter)
            }
            BBIReader::BigBed(reader) => {
                let batch_iter = BBIZoomQueryBatchIterator::from_bigbed(
                    reader,
                    region,
                    self.zoom_level,
                    batch_builder,
                    batch_size,
                    limit,
                );
                Ok(batch_iter)
            }
        }
    }
}

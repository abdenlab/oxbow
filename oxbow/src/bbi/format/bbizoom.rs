use std::io::{self, Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema as ArrowSchema;

pub use super::BBIReader;
use crate::bbi::batch_iterator::zoom::{BBIZoomBatchIterator, BBIZoomQueryBatchIterator};
use crate::bbi::model::zoom::field::DEFAULT_FIELD_NAMES;
use crate::bbi::model::zoom::BatchBuilder;

/// A scanner for the summary statistics from BBI file zoom level.
///
/// # Examples
///
/// ```no_run
/// use oxbow::bbi::format::bbizoom::Scanner;
/// use oxbow::bbi::BBIReader;
///
/// let mut fmt_reader = bigtools::BigWigRead::open_file("sample.bigWig").unwrap();
/// let info = fmt_reader.info();
/// let ref_names = info.chrom_info.iter().map(|c| c.name.clone()).collect();
/// let zoom_levels: Vec<u32> = info.zoom_headers.iter().map(|h| h.reduction_level).collect();
/// let scanner = Scanner::new(ref_names, zoom_levels[0]);
/// let batches = scanner.scan(BBIReader::BigWig(fmt_reader), None, None, Some(1000));
pub struct Scanner {
    ref_names: Vec<String>,
    zoom_level: u32,
}

impl Scanner {
    /// Creates a BBI zoom level scanner.
    pub fn new(ref_names: Vec<String>, zoom_level: u32) -> Self {
        Self {
            ref_names,
            zoom_level,
        }
    }

    /// Returns the reference sequence names.
    pub fn field_names(&self) -> Vec<String> {
        let fields = DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect();
        fields
    }

    /// Returns the Arrow schema.
    pub fn schema(&self, fields: Option<Vec<String>>) -> io::Result<ArrowSchema> {
        let batch_builder = BatchBuilder::new(&self.ref_names, fields, 0)?;
        Ok(batch_builder.get_arrow_schema())
    }
}

impl Scanner {
    /// Returns an iterator over record batches.
    pub fn scan<R: Read + Seek + Send + 'static>(
        &self,
        reader: BBIReader<R>,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<Box<dyn RecordBatchReader + Send>> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new(&self.ref_names, fields, batch_size)?;
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
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new(&self.ref_names, fields, batch_size)?;
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

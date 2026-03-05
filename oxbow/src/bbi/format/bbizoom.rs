use std::io::{self, Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::{Schema as ArrowSchema, SchemaRef};

pub use super::BBIReader;
use crate::batch::RecordBatchBuilder as _;
use crate::bbi::batch_iterator::zoom::{BBIZoomBatchIterator, BBIZoomQueryBatchIterator};
use crate::bbi::model::zoom::field::DEFAULT_FIELD_NAMES;
use crate::bbi::model::zoom::BatchBuilder;

/// A scanner for the summary statistics from BBI file zoom level.
///
/// Schema parameters (fields) are declared at construction time. Scan methods
/// accept only column projection, batch size, and limit.
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
/// let scanner = Scanner::new(ref_names, zoom_levels[0], None).unwrap();
/// let batches = scanner.scan(BBIReader::BigWig(fmt_reader), None, None, Some(1000));
pub struct Scanner {
    ref_names: Vec<String>,
    zoom_level: u32,
    fields: Option<Vec<String>>,
    schema: SchemaRef,
}

impl Scanner {
    /// Creates a BBI zoom level scanner.
    ///
    /// The schema is validated and cached at construction time.
    pub fn new(
        ref_names: Vec<String>,
        zoom_level: u32,
        fields: Option<Vec<String>>,
    ) -> io::Result<Self> {
        let batch_builder = BatchBuilder::new(&ref_names, fields.clone(), 0)?;
        let schema = batch_builder.schema();
        Ok(Self {
            ref_names,
            zoom_level,
            fields,
            schema,
        })
    }

    /// Returns the field names.
    pub fn field_names(&self) -> Vec<String> {
        DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self) -> &ArrowSchema {
        &self.schema
    }

    /// Builds a BatchBuilder applying column projection.
    fn build_batch_builder(
        &self,
        columns: Option<Vec<String>>,
        capacity: usize,
    ) -> io::Result<BatchBuilder> {
        match columns {
            None => BatchBuilder::new(&self.ref_names, self.fields.clone(), capacity),
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

                BatchBuilder::new(&self.ref_names, Some(cols), capacity)
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
    ) -> io::Result<Box<dyn RecordBatchReader + Send>> {
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
    ) -> io::Result<impl RecordBatchReader> {
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

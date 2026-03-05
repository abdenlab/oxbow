use std::io::{self, Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::{Schema as ArrowSchema, SchemaRef};
use bigtools::BigBedRead;

use crate::batch::RecordBatchBuilder as _;
use crate::bbi::batch_iterator::base::{BigBedBatchIterator, BigBedQueryBatchIterator};
use crate::bbi::model::base::BedSchema;
use crate::bbi::model::base::BatchBuilder;

/// A BigBed scanner.
///
/// Schema parameters (fields) are declared at construction time. Scan methods
/// accept only column projection, batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::bbi::format::bigbed::Scanner;
///
/// let mut fmt_reader = bigtools::BigBedRead::open_file("sample.bigBed").unwrap();
/// let info = fmt_reader.info();
///
/// let scanner = Scanner::new("bed12".parse().unwrap(), info.clone(), None).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
pub struct Scanner {
    bed_schema: BedSchema,
    info: bigtools::BBIFileInfo,
    fields: Option<Vec<String>>,
    schema: SchemaRef,
}

impl Scanner {
    /// Creates a BigBed scanner from a BED schema, BBI file info, and optional field names.
    ///
    /// The schema is validated and cached at construction time.
    pub fn new(
        bed_schema: BedSchema,
        info: bigtools::BBIFileInfo,
        fields: Option<Vec<String>>,
    ) -> io::Result<Self> {
        let batch_builder = BatchBuilder::new(bed_schema.clone(), fields.clone(), 0)?;
        let schema = batch_builder.schema();
        Ok(Self {
            bed_schema,
            info,
            fields,
            schema,
        })
    }

    /// Returns the field names.
    pub fn field_names(&self) -> Vec<String> {
        let fields = self.bed_schema.fields();
        fields.iter().map(|def| def.name.clone()).collect()
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
            None => BatchBuilder::new(self.bed_schema.clone(), self.fields.clone(), capacity),
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

                BatchBuilder::new(self.bed_schema.clone(), Some(cols), capacity)
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
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        let batch_iter = BigBedBatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches satisfying a genomic range query.
    pub fn scan_query<R: Read + Seek + Send + 'static>(
        &self,
        fmt_reader: BigBedRead<R>,
        region: noodles::core::Region,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;
        let batch_iter =
            BigBedQueryBatchIterator::new(fmt_reader, region, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

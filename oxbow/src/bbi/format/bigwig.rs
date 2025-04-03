use std::io::{self, Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema as ArrowSchema;
use bigtools::BigWigRead;

use crate::bbi::batch_iterator::base::{BigWigBatchIterator, BigWigQueryBatchIterator};
use crate::bbi::model::base::schema::BedSchema;
use crate::bbi::model::base::BatchBuilder;

/// A BigWig scanner.
///
/// # Examples
///
/// ```no_run
/// use oxbow::bbi::format::bigwig::Scanner;
/// let mut fmt_reader = bigtools::BigWigRead::open_file("sample.bigWig").unwrap();
/// let info = fmt_reader.info();
///
/// let scanner = Scanner::new(info.clone());
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {
    bed_schema: BedSchema,
    info: bigtools::BBIFileInfo,
}

impl Scanner {
    /// Creates a BigBed scanner from a BED schema and BBI file info.
    pub fn new(info: bigtools::BBIFileInfo) -> Self {
        let bed_schema = "bedGraph".parse().unwrap();
        Self { bed_schema, info }
    }

    /// Returns the field names.
    pub fn field_names(&self) -> Vec<String> {
        let fields = self.bed_schema.fields();
        fields.iter().map(|def| def.name.clone()).collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self, fields: Option<Vec<String>>) -> io::Result<ArrowSchema> {
        let batch_builder = BatchBuilder::new(self.bed_schema.clone(), fields, 0)?;
        Ok(batch_builder.get_arrow_schema())
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

    /// Returns an iterator over record batches.
    pub fn scan<R: Read + Seek>(
        &self,
        fmt_reader: BigWigRead<R>,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new(self.bed_schema.clone(), fields, batch_size)?;
        let batch_iter = BigWigBatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator over record batches satisfying a genomic range query.
    pub fn scan_query<R: Read + Seek + Send + 'static>(
        &self,
        fmt_reader: BigWigRead<R>,
        region: noodles::core::Region,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = BatchBuilder::new(self.bed_schema.clone(), fields, batch_size)?;
        let batch_iter =
            BigWigQueryBatchIterator::new(fmt_reader, region, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

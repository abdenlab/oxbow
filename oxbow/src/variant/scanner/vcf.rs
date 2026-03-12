use std::io::{BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::{Schema, SchemaRef};
use noodles::bgzf::VirtualPosition;
use noodles::csi::BinningIndex;

use crate::batch::RecordBatchBuilder as _;
use crate::util::query::{BgzfChunkReader, ByteRangeReader};
use crate::variant::model::field::DEFAULT_FIELD_NAMES;
use crate::variant::model::{BatchBuilder, GenotypeBy};
use crate::variant::scanner::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::OxbowError;

/// A VCF scanner.
///
/// Schema parameters (fields, info fields, genotype fields, samples, genotype_by)
/// are declared at construction time. Scan methods accept only column projection,
/// batch size, and limit.
///
/// # Examples
///
/// ```no_run
/// use oxbow::variant::scanner::vcf::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.vcf").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::vcf::io::Reader::new(inner);
/// let header = fmt_reader.read_header().unwrap();
///
/// let scanner = Scanner::new(header, None, None, None, None, None, None).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {
    header: noodles::vcf::Header,
    fields: Option<Vec<String>>,
    info_fields: Option<Vec<String>>,
    genotype_fields: Option<Vec<String>>,
    samples: Option<Vec<String>>,
    genotype_by: GenotypeBy,
    sample_prefix: Option<String>,
    schema: SchemaRef,
}

impl Scanner {
    /// Creates a VCF scanner from a VCF header and schema parameters.
    ///
    /// The schema is validated and cached at construction time.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        header: noodles::vcf::Header,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<GenotypeBy>,
        sample_prefix: Option<String>,
    ) -> crate::Result<Self> {
        let genotype_by = genotype_by.unwrap_or(GenotypeBy::Sample);
        let batch_builder = BatchBuilder::new(
            header.clone(),
            fields.clone(),
            info_fields.clone(),
            genotype_fields.clone(),
            samples.clone(),
            genotype_by.clone(),
            sample_prefix.clone(),
            0,
        )?;
        let schema = batch_builder.schema();
        Ok(Self {
            header,
            fields,
            info_fields,
            genotype_fields,
            samples,
            genotype_by,
            sample_prefix,
            schema,
        })
    }

    /// Returns a reference to the VCF header.
    pub fn header(&self) -> &noodles::vcf::Header {
        &self.header
    }

    /// Returns the reference sequence names.
    pub fn chrom_names(&self) -> Vec<String> {
        self.header
            .contigs()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect()
    }

    /// Returns the reference sequence names and lengths.
    pub fn chrom_sizes(&self) -> Vec<(String, u32)> {
        self.header
            .contigs()
            .iter()
            .map(|(name, r)| (name.to_string(), r.length().unwrap_or(0) as u32))
            .collect()
    }

    /// Returns the fixed field names.
    pub fn field_names(&self) -> Vec<String> {
        DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect()
    }

    /// Returns the INFO field names.
    pub fn info_field_names(&self) -> Vec<String> {
        self.header.infos().iter().map(|(k, _)| k.clone()).collect()
    }

    /// Returns the INFO field definitions.
    pub fn info_field_defs(&self) -> Vec<(String, String, String)> {
        use noodles::vcf::header::record::value::map::info::Number;
        self.header
            .infos()
            .iter()
            .map(|(k, def)| {
                let n = match def.number() {
                    Number::Count(n) => n.to_string(),
                    Number::AlternateBases => "A".to_string(),
                    Number::ReferenceAlternateBases => "R".to_string(),
                    Number::Samples => "G".to_string(),
                    Number::Unknown => ".".to_string(),
                };
                let ty = def.ty().to_string();
                (k.clone(), n, ty)
            })
            .collect()
    }

    /// Returns the FORMAT field names.
    pub fn genotype_field_names(&self) -> Vec<String> {
        self.header
            .formats()
            .iter()
            .map(|(k, _)| k.clone())
            .collect()
    }

    /// Returns the FORMAT field definitions.
    pub fn genotype_field_defs(&self) -> Vec<(String, String, String)> {
        use noodles::vcf::header::record::value::map::format::Number;
        self.header
            .formats()
            .iter()
            .map(|(k, def)| {
                let n = match def.number() {
                    Number::Count(n) => n.to_string(),
                    Number::AlternateBases => "A".to_string(),
                    Number::ReferenceAlternateBases => "R".to_string(),
                    Number::Samples => "G".to_string(),
                    Number::LocalAlternateBases => "LA".to_string(),
                    Number::LocalReferenceAlternateBases => "LR".to_string(),
                    Number::LocalSamples => "LG".to_string(),
                    Number::Ploidy => "P".to_string(),
                    Number::BaseModifications => "M".to_string(),
                    Number::Unknown => ".".to_string(),
                };
                let ty = def.ty().to_string();
                (k.clone(), n, ty)
            })
            .collect()
    }

    /// Returns the sample names.
    pub fn sample_names(&self) -> Vec<String> {
        self.header.sample_names().iter().cloned().collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self) -> &Schema {
        &self.schema
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
            None => BatchBuilder::new(
                self.header.clone(),
                self.fields.clone(),
                self.info_fields.clone(),
                self.genotype_fields.clone(),
                self.samples.clone(),
                self.genotype_by.clone(),
                self.sample_prefix.clone(),
                capacity,
            ),
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
                    return Err(OxbowError::invalid_input(format!(
                        "Unknown columns: {:?}. Available columns: {:?}",
                        unknown, schema_names
                    )));
                }

                // Fixed fields: intersect with declared
                let declared_field_names: Vec<String> = self.fields.clone().unwrap_or_else(|| {
                    DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect()
                });
                let projected_fields: Vec<String> = declared_field_names
                    .into_iter()
                    .filter(|name| cols.iter().any(|c| c.eq_ignore_ascii_case(name)))
                    .collect();

                // Include info only if "info" is in the requested columns
                let info_fields = if cols.iter().any(|c| c == "info") {
                    self.info_fields.clone()
                } else {
                    Some(vec![])
                };

                // Genotype/sample columns: filter based on genotype_by mode
                let (genotype_fields, samples) = match self.genotype_by {
                    GenotypeBy::Sample => {
                        // Non-fixed, non-"info" columns are sample names
                        let declared_samples: Vec<String> =
                            self.samples.clone().unwrap_or_else(|| self.sample_names());
                        let projected_samples: Vec<String> = declared_samples
                            .into_iter()
                            .filter(|name| {
                                let col_name = match &self.sample_prefix {
                                    Some(prefix) => format!("{}{}", prefix, name),
                                    None => name.clone(),
                                };
                                cols.iter().any(|c| c == &col_name)
                            })
                            .collect();
                        if projected_samples.is_empty() {
                            (Some(vec![]), Some(vec![]))
                        } else {
                            (self.genotype_fields.clone(), Some(projected_samples))
                        }
                    }
                    GenotypeBy::Field => {
                        // Non-fixed, non-"info" columns are genotype field names
                        let declared_gt_fields: Vec<String> = self
                            .genotype_fields
                            .clone()
                            .unwrap_or_else(|| self.genotype_field_names());
                        let projected_gt: Vec<String> = declared_gt_fields
                            .into_iter()
                            .filter(|name| cols.iter().any(|c| c == name))
                            .collect();
                        if projected_gt.is_empty() {
                            (Some(vec![]), Some(vec![]))
                        } else {
                            (Some(projected_gt), self.samples.clone())
                        }
                    }
                };

                BatchBuilder::new(
                    self.header.clone(),
                    Some(projected_fields),
                    info_fields,
                    genotype_fields,
                    samples,
                    self.genotype_by.clone(),
                    self.sample_prefix.clone(),
                    capacity,
                )
            }
        }
    }
}

impl Scanner {
    /// Returns an iterator yielding record batches.
    pub fn scan<R: BufRead>(
        &self,
        fmt_reader: noodles::vcf::io::Reader<R>,
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
        fmt_reader: noodles::vcf::io::Reader<R>,
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

        let reference_sequence_id =
            resolve_chrom_id(&self.header, &index, &reference_sequence_name)?;
        let chunks = index.query(reference_sequence_id, interval)?;
        let chunks = chunks.into_iter().map(|c| (c.start(), c.end())).collect();
        let bgzf_reader = fmt_reader.into_inner();
        let query_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::vcf::io::Reader::new(query_reader);
        let batch_iter = QueryBatchIterator::new(
            fmt_reader,
            self.header.clone(),
            reference_sequence_name,
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
        fmt_reader: noodles::vcf::io::Reader<R>,
        byte_ranges: Vec<(u64, u64)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let inner_reader = fmt_reader.into_inner();
        let range_reader = ByteRangeReader::new(inner_reader, byte_ranges);
        let fmt_reader = noodles::vcf::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches from specified virtual position ranges.
    pub fn scan_virtual_ranges<R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek>(
        &self,
        fmt_reader: noodles::vcf::io::Reader<R>,
        vpos_ranges: Vec<(VirtualPosition, VirtualPosition)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let batch_builder = self.build_batch_builder(columns, batch_size)?;

        let bgzf_reader = fmt_reader.into_inner();
        let range_reader = BgzfChunkReader::new(bgzf_reader, vpos_ranges);
        let fmt_reader = noodles::vcf::io::Reader::new(range_reader);
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }
}

fn resolve_chrom_id(
    header: &noodles::vcf::Header,
    index: &impl BinningIndex,
    chrom: &str,
) -> crate::Result<usize> {
    // For VCF, first try the index file's header, then try the source file's header.
    let id = index
        .header()
        .and_then(|index_header| {
            index_header
                .reference_sequence_names()
                .get_index_of(chrom.as_bytes())
        })
        .or_else(|| header.contigs().get_index_of(chrom));

    id.ok_or_else(|| {
        OxbowError::invalid_input(format!(
            "Reference sequence '{}' not found in Index header or VCF header.",
            chrom
        ))
    })
}

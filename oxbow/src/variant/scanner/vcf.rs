use std::io::{BufRead, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use noodles::bgzf::VirtualPosition;
use noodles::csi::BinningIndex;

use crate::util::query::{BgzfChunkReader, ByteRangeReader};
use crate::variant::model::{BatchBuilder, GenotypeBy, Model};
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
/// let scanner = Scanner::new(header, None, None, None, None, None).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {
    header: noodles::vcf::Header,
    model: Model,
}

impl Scanner {
    /// Creates a VCF scanner from a VCF header and schema parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        header: noodles::vcf::Header,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<GenotypeBy>,
    ) -> crate::Result<Self> {
        let model = Model::from_header(
            &header,
            fields,
            info_fields,
            genotype_fields,
            samples,
            genotype_by,
        )?;
        Ok(Self { header, model })
    }

    /// Creates a VCF scanner from a [`Model`].
    pub fn with_model(header: noodles::vcf::Header, model: Model) -> Self {
        Self { header, model }
    }

    /// Returns a reference to the [`Model`].
    pub fn model(&self) -> &Model {
        &self.model
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
        self.model.field_names()
    }

    /// Returns the INFO field names from the header.
    pub fn info_field_names(&self) -> Vec<String> {
        self.header.infos().iter().map(|(k, _)| k.clone()).collect()
    }

    /// Returns the INFO field definitions from the header.
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

    /// Returns the FORMAT field names from the header.
    pub fn genotype_field_names(&self) -> Vec<String> {
        self.header
            .formats()
            .iter()
            .map(|(k, _)| k.clone())
            .collect()
    }

    /// Returns the FORMAT field definitions from the header.
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

    /// Returns the sample names from the header.
    pub fn sample_names(&self) -> Vec<String> {
        self.header.sample_names().iter().cloned().collect()
    }

    /// Returns the Arrow schema.
    pub fn schema(&self) -> &Schema {
        self.model.schema()
    }

    /// Builds a BatchBuilder applying column projection.
    fn build_batch_builder(
        &self,
        columns: Option<Vec<String>>,
        capacity: usize,
    ) -> crate::Result<BatchBuilder> {
        match columns {
            None => BatchBuilder::from_model(&self.model, self.header.clone(), capacity),
            Some(cols) => {
                let projected = self.model.project(&cols)?;
                BatchBuilder::from_model(&projected, self.header.clone(), capacity)
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

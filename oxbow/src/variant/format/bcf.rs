use std::io::{self, Read, Seek};

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use noodles::csi::BinningIndex;

use crate::util::query::BgzfChunkReader;
use crate::variant::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::variant::model::field::DEFAULT_FIELD_NAMES;
use crate::variant::model::{BatchBuilder, GenotypeBy};

/// A BCF scanner.
///
/// # Examples
///
/// ```no_run
/// use oxbow::variant::format::bcf::Scanner;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.bcf").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::bcf::io::Reader::new(inner);
/// let header = fmt_reader.read_header().unwrap();
///
/// let scanner = Scanner::new(header);
/// let batches = scanner.scan(fmt_reader, None, None, None, None, None, None, Some(1000));
/// ```
pub struct Scanner {
    header: noodles::vcf::Header,
}

impl Scanner {
    /// Creates a BCF scanner from a VCF header.
    pub fn new(header: noodles::vcf::Header) -> Self {
        Self { header }
    }

    /// Returns the VCF header.
    pub fn header(&self) -> noodles::vcf::Header {
        self.header.clone()
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
    pub fn schema(
        &self,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<GenotypeBy>,
    ) -> io::Result<Schema> {
        let header = self.header();
        let genotype_by = genotype_by.unwrap_or(GenotypeBy::Sample);
        let batch_builder = BatchBuilder::new(
            header,
            fields,
            info_fields,
            genotype_fields,
            samples,
            genotype_by,
            0,
        )?;
        Ok(batch_builder.get_arrow_schema())
    }
}

impl Scanner {
    /// Returns an iterator yielding batches of records.
    ///
    /// The scan will begin at the current position of the reader and will
    /// move the cursor to the end of the last record scanned.
    #[allow(clippy::too_many_arguments)]
    pub fn scan<R: Read>(
        &self,
        fmt_reader: noodles::bcf::io::Reader<R>,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<GenotypeBy>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let genotype_by = genotype_by.unwrap_or(GenotypeBy::Sample);
        let header = self.header();
        let batch_builder = BatchBuilder::new(
            header,
            fields,
            info_fields,
            genotype_fields,
            samples,
            genotype_by,
            batch_size,
        )?;
        let batch_iter = BatchIterator::new(fmt_reader, batch_builder, batch_size, limit);
        Ok(batch_iter)
    }

    /// Returns an iterator yielding record batches satisfying a genomic range query.
    ///
    /// This operation requires a BGZF source and an Index.
    ///
    /// The scan will consume contiguous "chunks" of BGZF blocks and filter for
    /// records that overlap the given region. The cursor will move to the end
    /// of the last record scanned.
    #[allow(clippy::too_many_arguments)]
    pub fn scan_query<R: Read + Seek>(
        &self,
        fmt_reader: noodles::bcf::io::Reader<noodles::bgzf::Reader<R>>,
        region: noodles::core::Region,
        index: impl BinningIndex,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<GenotypeBy>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> io::Result<impl RecordBatchReader> {
        let genotype_by = genotype_by.unwrap_or(GenotypeBy::Sample);
        let batch_size = batch_size.unwrap_or(1024);
        let reference_sequence_name = region.name().to_string();
        let interval = region.interval();

        let batch_builder = BatchBuilder::new(
            self.header(),
            fields,
            info_fields,
            genotype_fields,
            samples,
            genotype_by,
            batch_size,
        )?;

        let reference_sequence_id =
            resolve_chrom_id(&self.header, &index, &reference_sequence_name)?;
        let chunks = index.query(reference_sequence_id, interval)?;
        let bgzf_reader = fmt_reader.into_inner();
        let query_reader = BgzfChunkReader::new(bgzf_reader, chunks);
        let fmt_reader = noodles::bcf::io::Reader::from(query_reader);
        let batch_iter = QueryBatchIterator::new(
            fmt_reader,
            self.header(),
            reference_sequence_name,
            interval,
            batch_builder,
            batch_size,
            limit,
        );
        Ok(batch_iter)
    }
}

fn resolve_chrom_id(
    header: &noodles::vcf::Header,
    index: &impl BinningIndex,
    chrom: &str,
) -> io::Result<usize> {
    // For BCF, first try the source file's header, then try the index file's header.
    let id = header.contigs().get_index_of(chrom).or_else(|| {
        index.header().and_then(|index_header| {
            index_header
                .reference_sequence_names()
                .get_index_of(chrom.as_bytes())
        })
    });

    id.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Reference sequence '{}' not found in VCF header or Index header.",
                chrom
            ),
        )
    })
}

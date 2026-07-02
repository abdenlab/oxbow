use std::io::{BufRead, Seek};
use std::pin::Pin;

use arrow::array::RecordBatchReader;
use arrow::datatypes::Schema;
use futures::{stream, Stream, StreamExt};
use noodles::bgzf::VirtualPosition;
use noodles::csi::BinningIndex;
use tokio::io::AsyncBufRead;

use crate::async_scanner::AsyncScanner;
use crate::batch::{Push, RecordBatchBuilder};
use crate::util::query::{BgzfChunkReader, ByteRangeReader};
use crate::variant::model::{BatchBuilder, GenotypeBy, Model};
use crate::variant::scanner::batch_iterator::{BatchIterator, QueryBatchIterator};
use crate::{CoordSystem, OxbowError, Region, Result, Select};

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
/// use oxbow::Select;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.vcf").map(BufReader::new).unwrap();
/// let mut fmt_reader = noodles::vcf::io::Reader::new(inner);
/// let header = fmt_reader.read_header().unwrap();
///
/// use oxbow::CoordSystem;
/// let scanner = Scanner::new(header, Select::All, Select::All, Select::All, None, Select::All, None, CoordSystem::OneClosed).unwrap();
/// let batches = scanner.scan(fmt_reader, None, None, Some(1000));
/// ```
pub struct Scanner {
    header: noodles::vcf::Header,
    model: Model,
}

impl Scanner {
    /// Creates a VCF scanner from a VCF header and schema parameters.
    ///
    /// - `header`: the VCF header, used for schema inference and validation.
    /// - `fields`: standard SAM field selection.
    /// - `info_fields`: INFO field selection.
    /// - `genotype_fields`: FORMAT field selection.
    /// - `genotype_by`: how to group genotype fields and samples.
    /// - `samples`: sample selection for genotype fields.
    /// - `samples_nested`: whether to nest sample-genotype columns under a single samples column.
    /// - `coord_system`: output coordinate system for position columns.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        header: noodles::vcf::Header,
        fields: Select,
        info_fields: Select,
        genotype_fields: Select,
        genotype_by: Option<GenotypeBy>,
        samples: Select,
        samples_nested: Option<bool>,
        coord_system: CoordSystem,
    ) -> crate::Result<Self> {
        let model = Model::from_header(
            &header,
            fields,
            info_fields,
            genotype_fields,
            genotype_by,
            samples,
            samples_nested,
            coord_system,
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
        region: Region,
        index: impl BinningIndex,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> crate::Result<impl RecordBatchReader> {
        let batch_size = batch_size.unwrap_or(1024);
        let region = region.to_noodles()?;
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

impl AsyncScanner for Scanner {
    fn scan<R: AsyncBufRead + Unpin + Send + 'static>(
        &self,
        reader: R,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> Pin<Box<dyn Stream<Item = Result<arrow_array::RecordBatch>>>> {
        let batch_size = batch_size.unwrap_or(1024);
        let vcf_reader = noodles::vcf::r#async::io::Reader::new(reader);
        let total: usize = 0;
        let batch_builder = BatchBuilder::from_model(&self.model, self.header.clone(), batch_size)
            .expect("BatchBuilder construction should not fail");

        let stream = stream::try_unfold(
            (vcf_reader, total, batch_builder, batch_size, limit),
            |state| async move {
                let (mut vcf_reader, mut total, mut batch_builder, batch_size, limit) = state;
                let mut count = 0;
                let mut record = noodles::vcf::Record::default();

                while count < batch_size && total < limit.unwrap_or(usize::MAX) {
                    match vcf_reader.read_record(&mut record).await {
                        Ok(0) => break,
                        Ok(_) => {
                            batch_builder.push(&record)?;
                            count += 1;
                            total += 1;
                        }
                        Err(e) => return Err(e.into()),
                    }
                }

                if count == 0 {
                    Ok(None)
                } else {
                    let batch = batch_builder.finish()?;
                    Ok(Some((
                        batch,
                        (vcf_reader, total, batch_builder, batch_size, limit),
                    )))
                }
            },
        )
        .boxed();

        stream
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::model::ModelBuilder;
    use futures::StreamExt;
    use noodles::vcf::header::record::value::map::{Contig, Info, Map};
    use noodles::vcf::Header;
    use tokio::io::BufReader as AsyncBufReader;

    fn create_test_scanner() -> Scanner {
        let header = Header::builder()
            .add_contig("sq0", Map::<Contig>::new())
            .add_info("DP", Map::<Info>::from("DP"))
            .add_sample_name("sample1")
            .build();
        let model = ModelBuilder::new(header.clone())
            .genotype_fields(Select::Omit)
            .samples(Select::Omit)
            .build()
            .unwrap();
        Scanner::with_model(header, model)
    }

    #[tokio::test]
    async fn test_async_scan_streams_batches() {
        let record_data = "\
sq0\t1\t.\tA\tT\t.\tPASS\tDP=30\t0/1
sq0\t2\t.\tC\tG\t.\tPASS\tDP=50\t1/1
sq0\t3\t.\tG\tA\t.\tPASS\tDP=10\t0/0
";
        let reader = AsyncBufReader::new(record_data.as_bytes());
        let scanner = create_test_scanner();
        let mut stream = <Scanner as AsyncScanner>::scan(&scanner, reader, Some(2), None);

        let batch = stream.next().await.expect("expected first batch").unwrap();
        assert_eq!(batch.num_rows(), 2);

        let batch = stream.next().await.expect("expected second batch").unwrap();
        assert_eq!(batch.num_rows(), 1);

        assert!(stream.next().await.is_none());
    }
}

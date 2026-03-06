use std::collections::HashMap;
use std::sync::Arc;

use crate::OxbowError;

use arrow::array::{ArrayRef, StructArray};
use arrow::datatypes::{DataType, Field as ArrowField, FieldRef, SchemaRef};
use arrow::error::ArrowError;
use arrow::record_batch::{RecordBatch, RecordBatchOptions};
use indexmap::IndexMap;
use noodles::vcf::variant::record::info::field::Value as InfoFieldValue;
use noodles::vcf::variant::record::samples::series::value::Value as SampleFieldValue;
use noodles::vcf::variant::record::samples::Sample;
use noodles::vcf::variant::record::samples::Series;

use crate::batch::{Push, RecordBatchBuilder};

use super::field::Push as _;
use super::field::{Field, FieldBuilder, DEFAULT_FIELD_NAMES};
use super::genotype::{GenotypeDef, SampleStructBuilder, SeriesStructBuilder};
use super::info::{InfoBuilder, InfoDef};

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum GenotypeBy {
    Sample,
    Field,
}

pub enum GenotypeDataBuilder {
    BySample(SampleStructBuilder),
    ByField(SeriesStructBuilder),
}

/// A builder for an Arrow record batch of variant calls.
pub struct BatchBuilder {
    schema: SchemaRef,
    row_count: usize,
    header: noodles::vcf::Header,
    info_defs: Vec<InfoDef>,
    genotype_defs: Vec<GenotypeDef>,
    sample_names: Vec<String>,
    genotype_by: GenotypeBy,
    field_builders: IndexMap<Field, FieldBuilder>,
    info_builders: IndexMap<InfoDef, InfoBuilder>,
    genotype_builders: IndexMap<String, GenotypeDataBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` for VCF/BCF records.
    pub fn new(
        header: noodles::vcf::Header,
        field_names: Option<Vec<String>>,
        info_field_names: Option<Vec<String>>,
        genotype_field_names: Option<Vec<String>>,
        sample_names: Option<Vec<String>>,
        genotype_by: GenotypeBy,
        capacity: usize,
    ) -> crate::Result<Self> {
        let ref_names = header
            .contigs()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect::<Vec<String>>();

        let default_field_names: Vec<String> = DEFAULT_FIELD_NAMES
            .iter()
            .map(|name| name.to_string())
            .collect();
        let fields: Vec<Field> = field_names
            .unwrap_or(default_field_names)
            .into_iter()
            .map(|name| name.parse())
            .collect::<Result<Vec<_>, _>>()?;
        let mut field_builders = IndexMap::new();
        for field in &fields {
            let builder = match field {
                Field::Chrom => FieldBuilder::with_refs(field.clone(), capacity, &ref_names)
                    .map_err(|e| crate::OxbowError::invalid_data(e.to_string()))?,
                _ => FieldBuilder::new(field.clone(), capacity),
            };
            field_builders.insert(field.clone(), builder);
        }

        let default_info_names: Vec<String> = header
            .infos()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect();
        let info_defs = info_field_names
            .unwrap_or(default_info_names)
            .into_iter()
            .filter_map(|name| {
                let info = header.infos().get(&name)?;
                Some(InfoDef::new(name, &info.number(), &info.ty()))
            })
            .collect::<Vec<InfoDef>>();
        let mut info_builders = IndexMap::new();
        for def in &info_defs {
            let builder = InfoBuilder::new(&def.ty);
            info_builders.insert(def.clone(), builder);
        }

        let default_sample_names = header
            .sample_names()
            .iter()
            .cloned()
            .collect::<Vec<String>>();
        let sample_names = sample_names.unwrap_or(default_sample_names);

        let default_genotype_names: Vec<String> = header
            .formats()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect();
        let genotype_defs = genotype_field_names
            .unwrap_or(default_genotype_names)
            .into_iter()
            .filter_map(|name| {
                let format = header.formats().get(&name)?;
                Some(GenotypeDef::new(name, &format.number(), &format.ty()))
            })
            .collect::<Vec<GenotypeDef>>();
        let genotype_builders = match genotype_by {
            GenotypeBy::Sample => {
                let sample_builders = sample_names
                    .iter()
                    .map(|sample_name| {
                        let builder = SampleStructBuilder::new(genotype_defs.clone());
                        (
                            sample_name.to_string(),
                            GenotypeDataBuilder::BySample(builder),
                        )
                    })
                    .collect::<IndexMap<String, GenotypeDataBuilder>>();
                sample_builders
            }
            GenotypeBy::Field => {
                let series_builders = genotype_defs
                    .iter()
                    .map(|def| {
                        let builder = SeriesStructBuilder::new(def.clone(), sample_names.clone());
                        (def.name.clone(), GenotypeDataBuilder::ByField(builder))
                    })
                    .collect::<IndexMap<String, GenotypeDataBuilder>>();
                series_builders
            }
        };

        // Build schema once
        let mut arrow_fields: Vec<ArrowField> =
            fields.iter().map(|f| f.get_arrow_field()).collect();
        if !info_defs.is_empty() {
            let nested: Vec<ArrowField> =
                info_defs.iter().map(|def| def.get_arrow_field()).collect();
            arrow_fields.push(ArrowField::new(
                "info",
                DataType::Struct(arrow::datatypes::Fields::from(nested)),
                true,
            ));
        }
        if !sample_names.is_empty() && !genotype_defs.is_empty() {
            match genotype_by {
                GenotypeBy::Sample => {
                    for (name, builder) in &genotype_builders {
                        let nested = match builder {
                            GenotypeDataBuilder::BySample(b) => b.get_arrow_fields(),
                            _ => panic!("Invalid builder type for sample: {:?}", name),
                        };
                        arrow_fields.push(ArrowField::new(
                            name,
                            DataType::Struct(arrow::datatypes::Fields::from(nested)),
                            true,
                        ));
                    }
                }
                GenotypeBy::Field => {
                    for (name, builder) in &genotype_builders {
                        let nested = match builder {
                            GenotypeDataBuilder::ByField(b) => b.get_arrow_fields(),
                            _ => panic!("Invalid builder type for field: {:?}", name),
                        };
                        arrow_fields.push(ArrowField::new(
                            name,
                            DataType::Struct(arrow::datatypes::Fields::from(nested)),
                            true,
                        ));
                    }
                }
            }
        }
        let schema = Arc::new(arrow::datatypes::Schema::new(arrow_fields));

        Ok(Self {
            schema,
            row_count: 0,
            header,
            info_defs,
            genotype_defs,
            sample_names,
            genotype_by,
            field_builders,
            info_builders,
            genotype_builders,
        })
    }

    pub fn header(&self) -> noodles::vcf::Header {
        self.header.clone()
    }
}

impl RecordBatchBuilder for BatchBuilder {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        let mut columns: Vec<ArrayRef> = self
            .field_builders
            .iter_mut()
            .map(|(_, builder)| builder.finish())
            .collect();

        // info (optional)
        if !self.info_defs.is_empty() {
            let info_arrays: Vec<(FieldRef, ArrayRef)> = self
                .info_builders
                .iter_mut()
                .map(|(def, builder)| {
                    let arrow_field = def.get_arrow_field();
                    let array_ref = builder.finish();
                    (Arc::new(arrow_field), array_ref)
                })
                .collect();
            let info = StructArray::from(info_arrays);
            columns.push(Arc::new(info));
        }

        // genotype data (optional)
        if !self.sample_names.is_empty() && !self.genotype_defs.is_empty() {
            match self.genotype_by {
                GenotypeBy::Sample => {
                    for (_, builder) in &mut self.genotype_builders {
                        let builder = match builder {
                            GenotypeDataBuilder::BySample(b) => b,
                            _ => unreachable!(),
                        };
                        columns.push(Arc::new(builder.finish()));
                    }
                }
                GenotypeBy::Field => {
                    for (_, builder) in &mut self.genotype_builders {
                        let builder = match builder {
                            GenotypeDataBuilder::ByField(b) => b,
                            _ => unreachable!(),
                        };
                        columns.push(Arc::new(builder.finish()));
                    }
                }
            }
        }

        let batch = if columns.is_empty() {
            RecordBatch::try_new_with_options(
                self.schema.clone(),
                columns,
                &RecordBatchOptions::new().with_row_count(Some(self.row_count)),
            )
        } else {
            RecordBatch::try_new(self.schema.clone(), columns)
        };
        self.row_count = 0;
        batch
    }
}

/// Iterate through an INFO field once and collect all successfully parsed values.
///
/// This is more resilient than calling `info.get()` per field, because `get()` scans
/// from the beginning each time and stops at the first tokenization error (e.g. `;;`
/// in malformed VCF files). By using `iter()`, we recover all parseable fields on
/// both sides of malformed tokens.
fn collect_info_fields<'a>(
    info: &'a dyn noodles::vcf::variant::record::Info,
    header: &'a noodles::vcf::Header,
) -> HashMap<&'a str, Option<InfoFieldValue<'a>>> {
    info.iter(header).filter_map(|result| result.ok()).collect()
}

/// Append a VCF record to the batch.
impl Push<&noodles::vcf::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::vcf::Record) -> crate::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record, &self.header)?;
        }

        // info (optional)
        if !self.info_defs.is_empty() {
            let info = record.info();
            let parsed = collect_info_fields(&info, &self.header);
            for (def, builder) in self.info_builders.iter_mut() {
                match parsed.get(def.name.as_str()) {
                    Some(Some(value)) => {
                        // info field was parsed successfully
                        builder.append_value(value)?;
                    }
                    Some(None) => {
                        // info field value is empty (missing, ".")
                        builder.append_null();
                    }
                    None => {
                        // info field is present in this record or parsing failed
                        builder.append_null();
                    }
                }
            }
        }

        // genotype data (optional)
        if !self.sample_names.is_empty() && !self.genotype_defs.is_empty() {
            let record_samples = record.samples();
            let keys: Vec<String> = self
                .genotype_defs
                .iter()
                .map(|def| def.name.clone())
                .collect();

            match self.genotype_by {
                GenotypeBy::Sample => {
                    for (sample_name, builder) in self.genotype_builders.iter_mut() {
                        let builder = match builder {
                            GenotypeDataBuilder::BySample(b) => b,
                            _ => {
                                return Err(OxbowError::invalid_data(format!(
                                    "Invalid builder type for sample: {:?}",
                                    sample_name
                                )));
                            }
                        };

                        // Get the sample named `sample_name` in this record.
                        let sample = match record_samples.get(&self.header, sample_name) {
                            Some(sample) => sample,
                            None => {
                                return Err(OxbowError::not_found(format!(
                                    "Sample not found: {}",
                                    sample_name
                                )))
                            }
                        };

                        // Collect the values for each desired field from the given sample in this record.
                        let data = keys
                            .iter()
                            .map(|key| {
                                let value = match sample.get(&self.header, key) {
                                    Some(Ok(v)) => v,
                                    _ => None,
                                };
                                (key.clone(), value)
                            })
                            .collect::<IndexMap<String, Option<SampleFieldValue>>>();

                        builder.push(data)?;
                    }
                }
                GenotypeBy::Field => {
                    for (key, builder) in self.genotype_builders.iter_mut() {
                        let builder = match builder {
                            GenotypeDataBuilder::ByField(b) => b,
                            _ => {
                                return Err(OxbowError::invalid_data(format!(
                                    "Invalid builder type for field: {:?}",
                                    key
                                )));
                            }
                        };

                        // Get the "series" for field `key` in this record.
                        let series = match record_samples.select(key) {
                            Some(result) => result,
                            None => {
                                // A series for this field could not be generated from this record.
                                // Presumably, because the field does not appear in any sample in
                                // this particular record. Store nulls.
                                let data = self
                                    .sample_names
                                    .iter()
                                    .map(|name| (name.clone(), None))
                                    .collect::<IndexMap<String, Option<SampleFieldValue>>>();

                                builder.push(data)?;
                                continue;
                            }
                        };

                        // Collect the values for field `key` from each desired sample in this record.
                        let data = self
                            .sample_names
                            .iter()
                            .map(|sample_name| {
                                let i = self
                                    .header
                                    .sample_names()
                                    .get_index_of(sample_name)
                                    .unwrap();
                                let maybe_result = series.get(&self.header, i).flatten();
                                let option = match maybe_result {
                                    Some(Ok(value)) => Some(value),
                                    _ => None,
                                };
                                (sample_name.clone(), option)
                            })
                            .collect::<IndexMap<String, Option<SampleFieldValue>>>();

                        builder.push(data)?;
                    }
                }
            }
        }
        self.row_count += 1;
        Ok(())
    }
}

/// Append a BCF record to the batch.
impl Push<&noodles::bcf::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::bcf::Record) -> crate::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record, &self.header)?;
        }

        // info (optional)
        if !self.info_defs.is_empty() {
            let info = record.info();
            let parsed = collect_info_fields(&info, &self.header);
            for (def, builder) in self.info_builders.iter_mut() {
                match parsed.get(def.name.as_str()) {
                    Some(Some(value)) => {
                        builder.append_value(value)?;
                    }
                    Some(None) => {
                        builder.append_null();
                    }
                    None => {
                        builder.append_null();
                    }
                }
            }
        }

        // genotype data (optional)
        if !self.sample_names.is_empty() && !self.genotype_defs.is_empty() {
            let record_samples = record.samples()?;
            let keys: Vec<String> = self
                .genotype_defs
                .iter()
                .map(|def| def.name.clone())
                .collect();

            match self.genotype_by {
                GenotypeBy::Sample => {
                    for (sample_name, builder) in self.genotype_builders.iter_mut() {
                        let builder = match builder {
                            GenotypeDataBuilder::BySample(b) => b,
                            _ => {
                                return Err(OxbowError::invalid_data(format!(
                                    "Invalid builder type for sample: {:?}",
                                    sample_name
                                )));
                            }
                        };

                        // Get the sample named `sample_name` in this record.
                        let sample = match record_samples.get(&self.header, sample_name) {
                            Some(sample) => sample,
                            None => {
                                return Err(OxbowError::not_found(format!(
                                    "Sample not found: {}",
                                    sample_name
                                )))
                            }
                        };

                        // Collect the values for each desired field from the given sample in this record.
                        let data = keys
                            .iter()
                            .map(|key| {
                                let value = match sample.get(&self.header, key) {
                                    Some(Ok(value)) => value,
                                    _ => None,
                                };
                                (key.to_string(), value)
                            })
                            .collect::<IndexMap<String, Option<SampleFieldValue>>>();

                        builder.push(data)?;
                    }
                }
                GenotypeBy::Field => {
                    for (key, builder) in self.genotype_builders.iter_mut() {
                        let builder = match builder {
                            GenotypeDataBuilder::ByField(b) => b,
                            _ => {
                                return Err(OxbowError::invalid_data(format!(
                                    "Invalid builder type for field: {:?}",
                                    key
                                )));
                            }
                        };

                        // Get the "series" for field `key` in this record.
                        let result = match record_samples.select(&self.header, key) {
                            Some(result) => result,
                            None => {
                                // A series for this field could not be generated from this record.
                                // Presumably, because the field does not appear in any sample in
                                // this particular record. Store nulls.
                                let data = self
                                    .sample_names
                                    .iter()
                                    .map(|name| (name.clone(), None))
                                    .collect::<IndexMap<String, Option<SampleFieldValue>>>();
                                builder.push(data)?;
                                continue;
                            }
                        };
                        let series = match result {
                            Ok(series) => series,
                            Err(e) => {
                                // The field is present in at least one sample in this record, but
                                // something went wrong when generating the series.
                                eprintln!("Error parsing series: {:?}", e);
                                let data = self
                                    .sample_names
                                    .iter()
                                    .map(|name| (name.clone(), None))
                                    .collect::<IndexMap<String, Option<SampleFieldValue>>>();
                                builder.push(data)?;
                                continue;
                            }
                        };

                        // Collect the values for field `key` from each desired sample in this record.
                        let data = self
                            .sample_names
                            .iter()
                            .map(|sample_name| {
                                let i = self
                                    .header
                                    .sample_names()
                                    .get_index_of(sample_name)
                                    .unwrap();
                                let maybe_result = series.get(&self.header, i).unwrap();
                                let option = match maybe_result {
                                    Some(Ok(value)) => Some(value),
                                    _ => None,
                                };
                                (sample_name.clone(), option)
                            })
                            .collect::<IndexMap<String, Option<SampleFieldValue>>>();

                        builder.push(data)?;
                    }
                }
            }
        }
        self.row_count += 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::vcf::header::record::value::map::{Contig, Format, Info, Map};
    use noodles::vcf::variant::io::Write;
    use noodles::vcf::variant::record_buf::samples::sample::Value;
    use noodles::vcf::variant::record_buf::samples::Keys;
    use noodles::vcf::variant::record_buf::Samples;
    use noodles::vcf::{Header, Record as VcfRecord};

    fn create_test_header() -> Header {
        let contig = Map::<Contig>::new();
        let info = Map::<Info>::from("DP");
        let format = Map::<Format>::from("GT");

        Header::builder()
            .add_contig("sq0", contig.clone())
            .add_info("DP", info)
            .add_format("GT", format)
            .add_sample_name("sample1")
            .add_sample_name("sample2")
            .build()
    }

    fn create_test_vcfrecord(header: &Header) -> VcfRecord {
        let keys: Keys = vec!["GT".to_string()].into_iter().collect();
        let samples = Samples::new(
            keys,
            vec![
                vec![Some(Value::from("0|0"))],
                vec![Some(Value::from("1/1"))],
            ],
        );
        let record_buf = noodles::vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name("sq0")
            .set_variant_start(Position::MIN)
            .set_reference_bases("A")
            .set_samples(samples)
            .build();
        let mut writer = noodles::vcf::io::Writer::new(Vec::new());
        writer.write_variant_record(header, &record_buf).unwrap();
        let buf = writer.into_inner();
        let record = noodles::vcf::Record::try_from(buf.as_slice()).unwrap();
        record
    }

    #[test]
    fn test_batch_builder_new() {
        let header = create_test_header();
        let batch_builder = BatchBuilder::new(
            header.clone(),
            None,
            None,
            None,
            None,
            GenotypeBy::Sample,
            10,
        )
        .unwrap();

        assert_eq!(batch_builder.header(), header);
        assert_eq!(batch_builder.sample_names.len(), 2);
        assert_eq!(batch_builder.genotype_by, GenotypeBy::Sample);
    }

    #[test]
    fn test_schema() {
        let header = create_test_header();
        let batch_builder =
            BatchBuilder::new(header, None, None, None, None, GenotypeBy::Sample, 10).unwrap();

        let schema = batch_builder.schema();
        assert!(schema.fields().len() > 0);
    }

    #[test]
    fn test_push_vcf_record() {
        let header = create_test_header();
        let mut batch_builder = BatchBuilder::new(
            header.clone(),
            None,
            None,
            None,
            None,
            GenotypeBy::Sample,
            10,
        )
        .unwrap();
        let record = create_test_vcfrecord(&header);

        assert!(batch_builder.push(&record).is_ok());
    }

    #[test]
    fn test_finish() {
        let header = create_test_header();
        let mut batch_builder = BatchBuilder::new(
            header.clone(),
            None,
            None,
            None,
            None,
            GenotypeBy::Sample,
            10,
        )
        .unwrap();
        let record = create_test_vcfrecord(&header);
        batch_builder.push(&record).unwrap();
        let record_batch = batch_builder.finish().unwrap();

        assert_eq!(record_batch.num_rows(), 1);
        assert_eq!(
            record_batch.num_columns(),
            batch_builder.schema().fields().len()
        );
    }

    #[test]
    fn test_push_vcf_record_with_double_semicolons_in_info() {
        use noodles::vcf::header::record::value::map::info::{Number, Type};

        // Build a header with multiple INFO fields
        let mut builder = Header::builder().add_contig("sq0", Map::<Contig>::new());
        for (id, ty, number) in [
            ("E_gnomAD", Type::Flag, Number::Count(0)),
            ("MA", Type::String, Number::Count(1)),
            ("MAF", Type::Float, Number::Count(1)),
        ] {
            let info = Map::<Info>::new(number, ty, "");
            builder = builder.add_info(id, info);
        }
        let header = builder.build();

        // Create a record with ;; in the INFO string (like Ensembl VCFs)
        let line = b"sq0\t1\t.\tA\t.\t.\t.\tE_gnomAD;;MA=C;MAF=0.123\n";
        let record = VcfRecord::try_from(&line[..]).unwrap();

        let mut batch_builder = BatchBuilder::new(
            header,
            None,
            None,
            None,
            Some(vec![]),
            GenotypeBy::Sample,
            10,
        )
        .unwrap();

        // Should not error - fields on both sides of ;; should be recovered
        batch_builder.push(&record).unwrap();
        let batch = batch_builder.finish().unwrap();

        assert_eq!(batch.num_rows(), 1);

        // INFO fields are nested under an "info" struct column
        let info_col = batch
            .column_by_name("info")
            .unwrap()
            .as_any()
            .downcast_ref::<StructArray>()
            .unwrap();

        // E_gnomAD is before ;; → should be parsed
        assert!(
            !info_col.column_by_name("E_gnomAD").unwrap().is_null(0),
            "E_gnomAD should not be null"
        );
        // MA is after ;; → should also be parsed thanks to iter() recovery
        assert!(
            !info_col.column_by_name("MA").unwrap().is_null(0),
            "MA should not be null"
        );
        // MAF is after ;; → should also be parsed
        assert!(
            !info_col.column_by_name("MAF").unwrap().is_null(0),
            "MAF should not be null"
        );
    }
}

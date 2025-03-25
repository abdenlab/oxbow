use std::io;
use std::sync::Arc;

use arrow::array::{ArrayRef, StructArray};
use arrow::datatypes::{DataType, Field as ArrowField, FieldRef, Schema};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use indexmap::IndexMap;
use noodles::vcf::variant::record::samples::series::value::Value as SampleFieldValue;
use noodles::vcf::variant::record::samples::Sample;
use noodles::vcf::variant::record::samples::Series;

use super::field::Push as _;
use super::field::{Field, FieldBuilder, DEFAULT_FIELD_NAMES};
use super::genotype::{GenotypeDef, SampleStructBuilder, SeriesStructBuilder};
use super::info::{InfoBuilder, InfoDef};

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
    header: noodles::vcf::Header,
    fields: Vec<Field>,
    info_defs: Vec<InfoDef>,
    genotype_defs: Vec<GenotypeDef>,
    sample_names: Vec<String>,
    genotype_by: GenotypeBy,
    field_builders: IndexMap<Field, FieldBuilder>,
    info_builders: IndexMap<InfoDef, InfoBuilder>,
    genotype_builders: IndexMap<String, GenotypeDataBuilder>,
}

impl BatchBuilder {
    pub fn new(
        header: noodles::vcf::Header,
        field_names: Option<Vec<String>>,
        info_field_names: Option<Vec<String>>,
        genotype_field_names: Option<Vec<String>>,
        sample_names: Option<Vec<String>>,
        genotype_by: GenotypeBy,
        capacity: usize,
    ) -> io::Result<Self> {
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
                    .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?,
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

        Ok(Self {
            header,
            fields,
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

    pub fn get_arrow_fields(&self) -> Vec<ArrowField> {
        // fixed fields
        let mut fields: Vec<ArrowField> = self
            .fields
            .iter()
            .map(|field| field.get_arrow_field())
            .collect();

        // info (optional)
        if !self.info_defs.is_empty() {
            let nested_fields: Vec<ArrowField> = self
                .info_builders
                .iter()
                .map(|(def, b)| b.get_arrow_field(&def.name))
                .collect();
            let info_field = ArrowField::new(
                "info",
                DataType::Struct(arrow::datatypes::Fields::from(nested_fields)),
                true,
            );
            fields.push(info_field);
        }

        // genotype data (optional)
        if !self.sample_names.is_empty() && !self.genotype_defs.is_empty() {
            match self.genotype_by {
                GenotypeBy::Sample => {
                    for (sample_name, builder) in &self.genotype_builders {
                        let nested_fields = match builder {
                            GenotypeDataBuilder::BySample(builder) => builder.get_arrow_fields(),
                            _ => panic!("Invalid builder type for sample: {:?}", sample_name),
                        };
                        let sample_field = ArrowField::new(
                            sample_name,
                            DataType::Struct(arrow::datatypes::Fields::from(nested_fields)),
                            true,
                        );
                        fields.push(sample_field);
                    }
                }
                GenotypeBy::Field => {
                    for (field_name, builder) in &self.genotype_builders {
                        let nested_fields = match builder {
                            GenotypeDataBuilder::ByField(builder) => builder.get_arrow_fields(),
                            _ => panic!("Invalid builder type for field: {:?}", field_name),
                        };
                        let field_field = ArrowField::new(
                            field_name,
                            DataType::Struct(arrow::datatypes::Fields::from(nested_fields)),
                            true,
                        );
                        fields.push(field_field);
                    }
                }
            }
        }

        fields
    }

    pub fn get_arrow_schema(&self) -> Schema {
        Schema::new(self.get_arrow_fields())
    }

    pub fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        let mut name_to_array: Vec<(&str, ArrayRef)> = self
            .field_builders
            .iter_mut()
            .map(|(field, builder)| {
                let name = field.name();
                (name, builder.finish())
            })
            .collect();

        // info (optional)
        if !self.info_defs.is_empty() {
            let info_arrays: Vec<(FieldRef, ArrayRef)> = self
                .info_builders
                .iter_mut()
                .map(|(def, builder)| {
                    let arrow_field = builder.get_arrow_field(&def.name);
                    let array_ref = builder.finish();
                    (Arc::new(arrow_field), array_ref)
                })
                .collect();
            let info = StructArray::from(info_arrays);
            name_to_array.push(("info", Arc::new(info)));
        }

        // genotype data (optional)
        if !self.sample_names.is_empty() && !self.genotype_defs.is_empty() {
            match self.genotype_by {
                GenotypeBy::Sample => {
                    for (sample_name, builder) in &mut self.genotype_builders {
                        let builder = match builder {
                            GenotypeDataBuilder::BySample(b) => b,
                            _ => panic!("Invalid builder type for sample: {:?}", sample_name),
                        };
                        let sample = builder.finish();
                        name_to_array.push((sample_name, Arc::new(sample)));
                    }
                }
                GenotypeBy::Field => {
                    for (field_name, builder) in &mut self.genotype_builders {
                        let builder = match builder {
                            GenotypeDataBuilder::ByField(b) => b,
                            _ => panic!("Invalid builder type for field: {:?}", field_name),
                        };
                        let series = builder.finish();
                        name_to_array.push((field_name, Arc::new(series)));
                    }
                }
            }
        }

        RecordBatch::try_from_iter(name_to_array)
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Append a VCF record to the batch.
impl Push<&noodles::vcf::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::vcf::Record) -> io::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record, &self.header)?;
        }

        // info (optional)
        if !self.info_defs.is_empty() {
            let info = record.info();
            for (def, builder) in self.info_builders.iter_mut() {
                match info.get(&self.header, &def.name) {
                    // info field is present in this record
                    Some(result) => {
                        match result {
                            // info field parsed successfully
                            Ok(Some(value)) => {
                                builder.append_value(&value)?;
                            }
                            // info field value is empty
                            Ok(None) => {
                                builder.append_null();
                            }
                            // info field could not be parsed
                            Err(_) => {
                                eprintln!("Error parsing INFO field: {:?}", def.name);
                                builder.append_null();
                            }
                        }
                    }
                    // info field is not present in this record
                    None => {
                        builder.append_null();
                    }
                };
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
                                return Err(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("Invalid builder type for sample: {:?}", sample_name),
                                ));
                            }
                        };

                        // Get the sample named `sample_name` in this record.
                        let sample = match record_samples.get(&self.header, sample_name) {
                            Some(sample) => sample,
                            None => {
                                return Err(io::Error::new(
                                    io::ErrorKind::NotFound,
                                    format!("Sample not found: {}", sample_name),
                                ))
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
                                return Err(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("Invalid builder type for field: {:?}", key),
                                ));
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
        Ok(())
    }
}

/// Append a BCF record to the batch.
impl Push<&noodles::bcf::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::bcf::Record) -> io::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record, &self.header)?;
        }

        // info (optional)
        if !self.info_defs.is_empty() {
            let info = record.info();
            for (def, builder) in self.info_builders.iter_mut() {
                match info.get(&self.header, &def.name) {
                    // info field is present in this record
                    Some(result) => {
                        match result {
                            // info field parsed successfully
                            Ok(Some(value)) => {
                                builder.append_value(&value)?;
                            }
                            // info field value is empty
                            Ok(None) => {
                                builder.append_null();
                            }
                            // info field could not be parsed
                            Err(_) => {
                                eprintln!("Error parsing INFO field: {:?}", def.name);
                                builder.append_null();
                            }
                        }
                    }
                    // info field is not present in this record
                    None => {
                        builder.append_null();
                    }
                };
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
                                return Err(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("Invalid builder type for sample: {:?}", sample_name),
                                ));
                            }
                        };

                        // Get the sample named `sample_name` in this record.
                        let sample = match record_samples.get(&self.header, sample_name) {
                            Some(sample) => sample,
                            None => {
                                return Err(io::Error::new(
                                    io::ErrorKind::NotFound,
                                    format!("Sample not found: {}", sample_name),
                                ))
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
                                return Err(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("Invalid builder type for field: {:?}", key),
                                ));
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
        Ok(())
    }
}

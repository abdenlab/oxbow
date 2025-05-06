use std::io;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, Float32Builder, GenericStringBuilder, Int32Builder, ListBuilder, StringArray,
    StringDictionaryBuilder,
};
use arrow::datatypes::{DataType, Field as ArrowField, Int32Type};
use arrow::error::ArrowError;

use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::Filters;
use noodles::vcf::variant::record::Ids;

use crate::util::reset_dictarray_builder;

pub const DEFAULT_FIELD_NAMES: [&str; 7] = ["chrom", "pos", "id", "ref", "alt", "qual", "filter"];

/// A variant (VCF/BCF) standard field.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum Field {
    Chrom,
    Pos,
    Id,
    Ref,
    Alt,
    Qual,
    Filter,
}

impl Field {
    pub fn name(&self) -> &str {
        match self {
            Self::Chrom => "chrom",
            Self::Pos => "pos",
            Self::Id => "id",
            Self::Ref => "ref",
            Self::Alt => "alt",
            Self::Qual => "qual",
            Self::Filter => "filter",
        }
    }

    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::Chrom => {
                DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8))
            }
            Self::Pos => DataType::Int32,
            Self::Id => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                DataType::List(Arc::new(item))
            }
            Self::Ref => DataType::Utf8,
            Self::Alt => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                DataType::List(Arc::new(item))
            }
            Self::Qual => DataType::Float32,
            Self::Filter => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                DataType::List(Arc::new(item))
            }
        }
    }

    pub fn get_arrow_field(&self) -> ArrowField {
        ArrowField::new(self.name(), self.arrow_type(), true)
    }
}

impl FromStr for Field {
    type Err = std::io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "chrom" => Ok(Self::Chrom),
            "pos" => Ok(Self::Pos),
            "id" => Ok(Self::Id),
            "ref" => Ok(Self::Ref),
            "alt" => Ok(Self::Alt),
            "qual" => Ok(Self::Qual),
            "filter" => Ok(Self::Filter),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Invalid field name: {}", s),
            )),
        }
    }
}

/// Builds an Arrow array (column) corresponding to a variant standard field.
pub enum FieldBuilder {
    Chrom(StringDictionaryBuilder<Int32Type>),
    Pos(Int32Builder),
    Id(ListBuilder<GenericStringBuilder<i32>>),
    Ref(GenericStringBuilder<i32>),
    Alt(ListBuilder<GenericStringBuilder<i32>>),
    Qual(Float32Builder),
    Filter(ListBuilder<GenericStringBuilder<i32>>),
}

impl FieldBuilder {
    pub fn new(field: Field, capacity: usize) -> Self {
        match field {
            Field::Chrom => Self::Chrom(StringDictionaryBuilder::<Int32Type>::new()),
            Field::Pos => Self::Pos(Int32Builder::with_capacity(capacity)),
            Field::Id => Self::Id(ListBuilder::with_capacity(
                GenericStringBuilder::new(),
                capacity,
            )),
            Field::Ref => Self::Ref(GenericStringBuilder::with_capacity(capacity, 1024)),
            Field::Alt => Self::Alt(ListBuilder::with_capacity(
                GenericStringBuilder::new(),
                capacity,
            )),
            Field::Qual => Self::Qual(Float32Builder::with_capacity(capacity)),
            Field::Filter => Self::Filter(ListBuilder::with_capacity(
                GenericStringBuilder::new(),
                capacity,
            )),
        }
    }

    pub fn with_refs(
        field: Field,
        capacity: usize,
        ref_names: &[String],
    ) -> Result<Self, ArrowError> {
        let field = match field {
            Field::Chrom => {
                let refs = StringArray::from(ref_names.to_owned());
                Self::Chrom(StringDictionaryBuilder::<Int32Type>::new_with_dictionary(
                    capacity, &refs,
                )?)
            }
            _ => {
                return Err(ArrowError::InvalidArgumentError(format!(
                    "Field {:?} does not require reference names",
                    field
                )))
            }
        };
        Ok(field)
    }

    pub fn finish(&mut self) -> ArrayRef {
        match self {
            Self::Chrom(builder) => {
                let array = reset_dictarray_builder(builder);
                Arc::new(array)
            }
            Self::Pos(builder) => Arc::new(builder.finish()),
            Self::Id(builder) => Arc::new(builder.finish()),
            Self::Ref(builder) => Arc::new(builder.finish()),
            Self::Alt(builder) => Arc::new(builder.finish()),
            Self::Qual(builder) => Arc::new(builder.finish()),
            Self::Filter(builder) => Arc::new(builder.finish()),
        }
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T, header: &noodles::vcf::Header) -> io::Result<()>;
}

/// Append a field value from a VCF record to the column.
impl Push<&noodles::vcf::Record> for FieldBuilder {
    fn push(
        &mut self,
        record: &noodles::vcf::Record,
        header: &noodles::vcf::Header,
    ) -> io::Result<()> {
        match self {
            Self::Chrom(builder) => {
                let rname = record.reference_sequence_name();
                builder.append_value(rname);
            }
            Self::Pos(builder) => {
                builder.append_option(
                    record
                        .variant_start()
                        .transpose()
                        .unwrap_or(None)
                        .map(|x| x.get() as i32),
                );
            }
            Self::Id(builder) => {
                let ids = record
                    .ids()
                    .iter()
                    .map(|id| id.to_string())
                    .collect::<Vec<String>>();
                for id in ids {
                    builder.values().append_value(id);
                }
                builder.append(true);
            }
            Self::Ref(builder) => {
                builder.append_value(record.reference_bases());
            }
            Self::Alt(builder) => {
                let alt_bases = record
                    .alternate_bases()
                    .iter()
                    .map(|alt| alt.map(|a| a.to_string()))
                    .collect::<Result<Vec<String>, _>>()?;
                for alt in alt_bases {
                    builder.values().append_value(alt);
                }
                builder.append(true);
            }
            Self::Qual(builder) => {
                builder.append_option(record.quality_score().transpose().unwrap_or(None));
            }
            Self::Filter(builder) => {
                let filters = record
                    .filters()
                    .iter(header)
                    .map(|filter| filter.map(|f| f.to_string()))
                    .collect::<Result<Vec<String>, _>>()?;

                if filters.is_empty() {
                    // Something went wrong. There should be at least one value.
                    // Make the entry invalid.
                    builder.append(false);
                } else {
                    match filters[0].as_str() {
                        "PASS" => {
                            // No filters failed. We store an empty list.
                            // The entry is valid.
                            builder.append(true);
                        }
                        "." => {
                            // Unknown filter state.
                            // Make the entry invalid.
                            builder.append(false);
                        }
                        _ => {
                            // Some known filters failed.
                            // The entry is valid.
                            for filter in filters {
                                builder.values().append_value(filter);
                            }
                            builder.append(true);
                        }
                    }
                };
            }
        }
        Ok(())
    }
}

/// Append a field value from a BCF record to the column.
impl Push<&noodles::bcf::Record> for FieldBuilder {
    fn push(
        &mut self,
        record: &noodles::bcf::Record,
        header: &noodles::vcf::Header,
    ) -> io::Result<()> {
        match self {
            Self::Chrom(builder) => {
                let id = record.reference_sequence_id()?;
                let rname = header
                    .contigs()
                    .get_index(id)
                    .map(|(name, _)| name.to_string());
                builder.append_option(rname);
            }
            Self::Pos(builder) => {
                builder.append_option(
                    record
                        .variant_start()
                        .transpose()
                        .unwrap_or(None)
                        .map(|x| x.get() as i32),
                );
            }
            Self::Id(builder) => {
                let ids = record
                    .ids()
                    .iter()
                    .map(|id| id.to_string())
                    .collect::<Vec<String>>();
                for id in ids {
                    builder.values().append_value(id);
                }
                builder.append(true);
            }
            Self::Ref(builder) => {
                let reference_bases = record.reference_bases();
                let ref_bytes = reference_bases.as_ref();
                let ref_str = String::from_utf8(ref_bytes.to_vec())
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                builder.append_value(ref_str);
            }
            Self::Alt(builder) => {
                let alt_bases = record
                    .alternate_bases()
                    .iter()
                    .map(|alt| alt.map(|a| a.to_string()))
                    .collect::<Result<Vec<String>, _>>()?;
                for alt in alt_bases {
                    builder.values().append_value(alt);
                }
                builder.append(true);
            }
            Self::Qual(builder) => {
                builder.append_option(record.quality_score()?);
            }
            Self::Filter(builder) => {
                let filters = record
                    .filters()
                    .iter(header)
                    .map(|filter| filter.map(|f| f.to_string()))
                    .collect::<Result<Vec<String>, _>>()?;

                if filters.is_empty() {
                    // Something went wrong. There should be at least one value.
                    // Make the entry invalid.
                    builder.append(false);
                } else {
                    match filters[0].as_str() {
                        "PASS" => {
                            // No filters failed. We store an empty list.
                            // The entry is valid.
                            builder.append(true);
                        }
                        "." => {
                            // Unknown filter state.
                            // Make the entry invalid.
                            builder.append(false);
                        }
                        _ => {
                            // Some known filters failed.
                            // The entry is valid.
                            for filter in filters {
                                builder.values().append_value(filter);
                            }
                            builder.append(true);
                        }
                    }
                };
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field_arrow_type() {
        for field in [
            Field::Chrom,
            Field::Pos,
            Field::Id,
            Field::Ref,
            Field::Alt,
            Field::Qual,
            Field::Filter,
        ] {
            let mut builder = FieldBuilder::new(field.clone(), 10);
            let data_type = builder.finish().data_type().clone();
            assert_eq!(field.arrow_type(), data_type);
        }
    }

    #[test]
    fn test_field_from_str() {
        assert_eq!(Field::from_str("chrom").unwrap(), Field::Chrom);
        assert_eq!(Field::from_str("pos").unwrap(), Field::Pos);
        assert_eq!(Field::from_str("id").unwrap(), Field::Id);
        assert_eq!(Field::from_str("ref").unwrap(), Field::Ref);
        assert_eq!(Field::from_str("alt").unwrap(), Field::Alt);
        assert_eq!(Field::from_str("qual").unwrap(), Field::Qual);
        assert_eq!(Field::from_str("filter").unwrap(), Field::Filter);
        assert!(Field::from_str("invalid").is_err());
    }

    #[test]
    fn test_field_builder_with_refs() {
        let ref_names = vec!["chr1".to_string(), "chr2".to_string()];
        let builder = FieldBuilder::with_refs(Field::Chrom, 10, &ref_names).unwrap();
        if let FieldBuilder::Chrom(mut b) = builder {
            let arr = b.finish();
            for ref_name in &ref_names {
                match arr.lookup_key(ref_name) {
                    Some(_) => (),
                    None => panic!("Reference name '{}' not found in dictionary", ref_name),
                }
            }
        }
        let result = FieldBuilder::with_refs(Field::Pos, 10, &ref_names);
        assert!(result.is_err());
    }

    #[test]
    fn test_field_builder_push() {
        for field in [
            Field::Chrom,
            Field::Pos,
            Field::Id,
            Field::Ref,
            Field::Alt,
            Field::Qual,
            Field::Filter,
        ] {
            let mut builder = FieldBuilder::new(field, 10);
            let record = noodles::vcf::Record::default();
            let header = noodles::vcf::Header::default();
            assert!(builder.push(&record, &header).is_ok());
        }
    }
}

use std::io;
use std::sync::Arc;

use arrow::array::{ArrayRef, StructArray};
use arrow::datatypes::{DataType, Field as ArrowField, FieldRef, Schema};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use indexmap::IndexMap;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::Data;

use super::field::Push as _;
use super::field::{Field, FieldBuilder, DEFAULT_FIELD_NAMES};
use super::tag::{TagBuilder, TagDef};

/// A builder for an Arrow record batch of alignments.
pub struct BatchBuilder {
    header: noodles::sam::Header,
    fields: Vec<Field>,
    tag_defs: Vec<TagDef>,
    field_builders: IndexMap<Field, FieldBuilder>,
    tag_builders: IndexMap<TagDef, TagBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` for SAM/BAM records.
    pub fn new(
        header: noodles::sam::Header,
        field_names: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        capacity: usize,
    ) -> io::Result<Self> {
        let ref_names = header
            .reference_sequences()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect::<Vec<String>>();

        let default_field_names: Vec<String> = DEFAULT_FIELD_NAMES
            .into_iter()
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
                Field::Rname | Field::Rnext => {
                    FieldBuilder::with_refs(field.clone(), capacity, &ref_names)
                        .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?
                }
                _ => FieldBuilder::new(field.clone(), capacity),
            };
            field_builders.insert(field.clone(), builder);
        }

        let tag_defs: Vec<TagDef> = tag_defs
            .unwrap_or_default()
            .into_iter()
            .map(TagDef::try_from)
            .collect::<Result<Vec<_>, _>>()?;
        let mut tag_builders = IndexMap::new();
        for tag in tag_defs.clone() {
            let builder = TagBuilder::new(&tag.ty);
            tag_builders.insert(tag, builder);
        }

        Ok(Self {
            header,
            fields,
            tag_defs,
            field_builders,
            tag_builders,
        })
    }

    pub fn header(&self) -> noodles::sam::Header {
        self.header.clone()
    }

    pub fn get_arrow_fields(&self) -> Vec<ArrowField> {
        // fixed fields
        let mut fields: Vec<ArrowField> = self
            .fields
            .iter()
            .map(|field| field.get_arrow_field())
            .collect();

        // tags (optional)
        if !self.tag_defs.is_empty() {
            let nested_fields: Vec<ArrowField> = self
                .tag_defs
                .iter()
                .map(|def| def.get_arrow_field())
                .collect();
            let tag_field = ArrowField::new(
                "tags",
                DataType::Struct(arrow::datatypes::Fields::from(nested_fields)),
                true,
            );
            fields.push(tag_field);
        }

        fields
    }

    pub fn get_arrow_schema(&self) -> Schema {
        Schema::new(self.get_arrow_fields())
    }

    pub fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        // fixed fields
        let mut name_to_array: Vec<(&str, ArrayRef)> = self
            .field_builders
            .iter_mut()
            .map(|(field, builder)| {
                let name = field.name();
                (name, builder.finish())
            })
            .collect();

        // tags (optional)
        if !self.tag_defs.is_empty() {
            let tag_arrays: Vec<(FieldRef, ArrayRef)> = self
                .tag_builders
                .iter_mut()
                .map(|(def, builder)| {
                    let arrow_field = def.get_arrow_field();
                    (Arc::new(arrow_field), builder.finish())
                })
                .collect();
            let tags = StructArray::from(tag_arrays);
            name_to_array.push(("tags", Arc::new(tags)));
        }

        RecordBatch::try_from_iter(name_to_array)
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Append a SAM record to the batch.
impl Push<&noodles::sam::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::sam::Record) -> io::Result<()> {
        // fixed fields
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record, &self.header)?;
        }

        // tags (optional)
        if !self.tag_defs.is_empty() {
            let data = record.data();
            for (def, builder) in self.tag_builders.iter_mut() {
                let tag = Tag::from(def.into_bytes());
                match data.get(&tag) {
                    // tag is present in this record
                    Some(result) => {
                        match result {
                            // tag parsed successfully
                            Ok(value) => {
                                builder.append_value(&value)?;
                            }
                            // tag could not be parsed
                            Err(_) => {
                                eprintln!("Error parsing tag: {:?}", tag);
                                builder.append_null();
                            }
                        }
                    }
                    // tag is not present in this record
                    None => {
                        builder.append_null();
                    }
                };
            }
        }
        Ok(())
    }
}

/// Append a BAM record to the batch.
impl Push<&noodles::bam::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::bam::Record) -> io::Result<()> {
        // fixed fields
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record, &self.header)?;
        }

        // tags (optional)
        if !self.tag_defs.is_empty() {
            let data = record.data();
            for (def, builder) in self.tag_builders.iter_mut() {
                let tag = Tag::from(def.into_bytes());
                match data.get(&tag) {
                    // tag is present in this record
                    Some(result) => {
                        match result {
                            // tag parsed successfully
                            Ok(value) => {
                                builder.append_value(&value)?;
                            }
                            // tag could not be parsed
                            Err(_) => {
                                eprintln!("Error parsing tag: {:?}", tag);
                                builder.append_null();
                            }
                        }
                    }
                    // tag is not present in this record
                    None => {
                        builder.append_null();
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
    fn test_batch_builder_new() {
        let header = noodles::sam::Header::default();
        let field_names = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let tag_defs = Some(vec![("NM".to_string(), "i".to_string())]);
        let capacity = 10;

        let batch_builder =
            BatchBuilder::new(header.clone(), field_names, tag_defs, capacity).unwrap();
        assert_eq!(batch_builder.header(), header);
        assert_eq!(batch_builder.fields.len(), 2);
        assert_eq!(batch_builder.tag_defs.len(), 1);
    }

    #[test]
    fn test_get_arrow_schema() {
        let header = noodles::sam::Header::default();
        let field_names = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let tag_defs = Some(vec![("NM".to_string(), "i".to_string())]);
        let capacity = 10;

        let batch_builder = BatchBuilder::new(header, field_names, tag_defs, capacity).unwrap();
        let schema = batch_builder.get_arrow_schema();
        assert_eq!(schema.fields().len(), 3);
        assert_eq!(schema.field(0).name(), "qname");
        assert_eq!(schema.field(1).name(), "flag");
        assert_eq!(schema.field(2).name(), "tags");
        assert_eq!(
            schema.field(2).data_type(),
            &DataType::Struct(vec![ArrowField::new("NM", DataType::Int64, true)].into())
        );
    }

    #[test]
    fn test_push_sam_record() {
        let header = noodles::sam::Header::default();
        let field_names = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let tag_defs = Some(vec![("NM".to_string(), "i".to_string())]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(header, field_names, tag_defs, capacity).unwrap();

        let record = noodles::sam::Record::default();
        let result = batch_builder.push(&record);
        assert!(result.is_ok());
    }

    #[test]
    fn test_push_bam_record() {
        let header = noodles::sam::Header::default();
        let field_names = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let tag_defs = Some(vec![("NM".to_string(), "i".to_string())]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(header, field_names, tag_defs, capacity).unwrap();

        let record = noodles::bam::Record::default();
        let result = batch_builder.push(&record);
        assert!(result.is_ok());
    }

    #[test]
    fn test_finish() {
        let header = noodles::sam::Header::default();
        let field_names = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let tag_defs = Some(vec![("NM".to_string(), "i".to_string())]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(header, field_names, tag_defs, capacity).unwrap();

        let record = noodles::sam::Record::default();
        batch_builder.push(&record).unwrap();
        let record_batch = batch_builder.finish();
        assert!(record_batch.is_ok());

        let record_batch = record_batch.unwrap();
        assert_eq!(record_batch.num_columns(), 3);
        assert_eq!(record_batch.num_rows(), 1);
    }
}

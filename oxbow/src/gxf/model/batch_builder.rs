use std::io;
use std::sync::Arc;

use arrow::array::{ArrayRef, StructArray};
use arrow::datatypes::FieldRef;
use arrow::datatypes::{DataType, Field as ArrowField, Schema};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use indexmap::IndexMap;

use crate::gxf::model::attribute::{AttributeBuilder, AttributeDef, AttributeValue};
use crate::gxf::model::field::Push as _;
use crate::gxf::model::field::{Field, FieldBuilder, DEFAULT_FIELD_NAMES};

/// A builder for an Arrow record batch of GXF (GTF/GFF) features.
pub struct BatchBuilder {
    fields: Vec<Field>,
    attr_defs: Vec<AttributeDef>,
    field_builders: IndexMap<Field, FieldBuilder>,
    attr_builders: IndexMap<AttributeDef, AttributeBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` for GTF/GFF records.
    pub fn new(
        field_names: Option<Vec<String>>,
        attr_defs: Option<Vec<(String, String)>>,
        capacity: usize,
    ) -> io::Result<Self> {
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
            let builder = FieldBuilder::new(field.clone(), capacity);
            field_builders.insert(field.clone(), builder);
        }

        let attr_defs: Vec<AttributeDef> = attr_defs
            .unwrap_or_default()
            .into_iter()
            .map(AttributeDef::try_from)
            .collect::<Result<Vec<_>, _>>()?;
        let mut attr_builders = IndexMap::new();
        for def in &attr_defs {
            let builder = AttributeBuilder::new(&def.ty);
            attr_builders.insert(def.clone(), builder);
        }

        Ok(Self {
            fields,
            attr_defs,
            field_builders,
            attr_builders,
        })
    }

    pub fn get_arrow_fields(&self) -> Vec<ArrowField> {
        // fixed fields
        let mut fields: Vec<ArrowField> = self
            .fields
            .iter()
            .map(|field| field.get_arrow_field())
            .collect();

        // attributes (optional)
        if !self.attr_defs.is_empty() {
            let nested_fields: Vec<ArrowField> = self
                .attr_defs
                .iter()
                .map(|def| def.get_arrow_field())
                .collect();
            let tag_field = ArrowField::new(
                "attributes",
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

        // attributes (optional)
        if !self.attr_defs.is_empty() {
            let attr_arrays: Vec<(FieldRef, ArrayRef)> = self
                .attr_builders
                .iter_mut()
                .map(|(def, builder)| {
                    let arrow_field = def.get_arrow_field();
                    (Arc::new(arrow_field), builder.finish())
                })
                .collect();
            let attrs = StructArray::from(attr_arrays);
            name_to_array.push(("attributes", Arc::new(attrs)));
        }

        RecordBatch::try_from_iter(name_to_array)
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Append a GFF record to the batch.
impl<'a> Push<&'a noodles::gff::Record<'a>> for BatchBuilder {
    fn push(&mut self, record: &noodles::gff::Record) -> io::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }

        // attributes (optional)
        if !self.attr_defs.is_empty() {
            let attrs = record.attributes();
            for (def, builder) in self.attr_builders.iter_mut() {
                match attrs.get(&def.name) {
                    // attribute is present in this record
                    Some(result) => {
                        match result {
                            // attribute parsed successfully
                            Ok(value) => {
                                let value = AttributeValue::from(&value);
                                builder.append_value(&value)?;
                            }
                            // attribute could not be parsed
                            Err(_) => {
                                eprintln!("Error parsing tag: {:?}", def.name);
                                builder.append_null();
                            }
                        }
                    }
                    // attribute is not present in this record
                    None => {
                        builder.append_null();
                    }
                };
            }
        }
        Ok(())
    }
}

/// Append a GTF record to the batch.
impl Push<&noodles::gtf::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::gtf::Record) -> io::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }

        // attributes (optional)
        if !self.attr_defs.is_empty() {
            let attrs = record.attributes();
            for (def, builder) in self.attr_builders.iter_mut() {
                match attrs.get(&def.name) {
                    Some(value) => {
                        let value = AttributeValue::String(value.to_string());
                        builder.append_value(&value)?;
                    }
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
    use arrow::datatypes::{DataType, Field as ArrowField};

    fn gff_recordbuf_to_line(record_buf: &noodles::gff::RecordBuf) -> noodles::gff::Line {
        let mut writer = noodles::gff::io::Writer::new(Vec::new());
        writer.write_record(record_buf).unwrap();
        let buf = writer.into_inner();
        let mut reader = noodles::gff::Reader::new(std::io::Cursor::new(buf.as_slice()));
        let line = reader.lines().next().unwrap().unwrap();
        line
    }

    #[test]
    fn test_batch_builder_new() {
        let field_names = Some(vec!["seqid".to_string(), "source".to_string()]);
        let attr_defs = Some(vec![("gene_id".to_string(), "String".to_string())]);
        let capacity = 10;

        let batch_builder = BatchBuilder::new(field_names, attr_defs, capacity).unwrap();

        assert_eq!(batch_builder.fields.len(), 2);
        assert_eq!(batch_builder.attr_defs.len(), 1);
    }

    #[test]
    fn test_get_arrow_schema() {
        let field_names = Some(vec!["seqid".to_string(), "source".to_string()]);
        let attr_defs = Some(vec![("gene_id".to_string(), "String".to_string())]);
        let capacity = 10;

        let batch_builder = BatchBuilder::new(field_names, attr_defs, capacity).unwrap();
        let schema = batch_builder.get_arrow_schema();

        assert_eq!(schema.fields().len(), 3); // 2 fields + 1 attributes struct
        assert_eq!(schema.field(0).name(), "seqid");
        assert_eq!(schema.field(1).name(), "source");
        assert_eq!(schema.field(2).name(), "attributes");
        assert_eq!(
            schema.field(2).data_type(),
            &DataType::Struct(vec![ArrowField::new("gene_id", DataType::Utf8, true),].into())
        );
    }

    #[test]
    fn test_push_gff_record() {
        let field_names = Some(vec!["seqid".to_string(), "source".to_string()]);
        let attr_defs = Some(vec![("gene_id".to_string(), "String".to_string())]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(field_names, attr_defs, capacity).unwrap();

        let record_buf = noodles::gff::RecordBuf::default();
        let line = gff_recordbuf_to_line(&record_buf);
        let record = line.as_record().unwrap().unwrap();
        let result = batch_builder.push(&record);
        assert!(result.is_ok());
    }

    #[test]
    fn test_push_gtf_record() {
        let field_names = Some(vec!["seqid".to_string(), "source".to_string()]);
        let attr_defs = Some(vec![("gene_id".to_string(), "String".to_string())]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(field_names, attr_defs, capacity).unwrap();

        let record = noodles::gtf::Record::default();
        let result = batch_builder.push(&record);
        assert!(result.is_ok());
    }

    #[test]
    fn test_finish() {
        let field_names = Some(vec!["seqid".to_string(), "source".to_string()]);
        let attr_defs = Some(vec![("gene_id".to_string(), "String".to_string())]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(field_names, attr_defs, capacity).unwrap();

        let record = noodles::gtf::Record::default();
        batch_builder.push(&record).unwrap();
        let record_batch = batch_builder.finish();
        assert!(record_batch.is_ok());

        let record_batch = record_batch.unwrap();
        assert_eq!(record_batch.num_columns(), 3);
        assert_eq!(record_batch.num_rows(), 1);
    }
}

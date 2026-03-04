use std::io;
use std::sync::Arc;

use arrow::array::{ArrayRef, StructArray};
use arrow::datatypes::FieldRef;
use arrow::datatypes::{DataType, Field as ArrowField, SchemaRef};
use arrow::error::ArrowError;
use arrow::record_batch::{RecordBatch, RecordBatchOptions};
use indexmap::IndexMap;

use crate::gxf::model::attribute::{AttributeBuilder, AttributeDef, AttributeValue};
use crate::gxf::model::field::Push as _;
use crate::gxf::model::field::{Field, FieldBuilder, DEFAULT_FIELD_NAMES};

/// A builder for an Arrow record batch of GXF (GTF/GFF) features.
pub struct BatchBuilder {
    schema: SchemaRef,
    row_count: usize,
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

        // Build schema once
        let mut arrow_fields: Vec<ArrowField> =
            fields.iter().map(|f| f.get_arrow_field()).collect();
        if !attr_defs.is_empty() {
            let nested_fields: Vec<ArrowField> =
                attr_defs.iter().map(|def| def.get_arrow_field()).collect();
            let attr_field = ArrowField::new(
                "attributes",
                DataType::Struct(arrow::datatypes::Fields::from(nested_fields)),
                true,
            );
            arrow_fields.push(attr_field);
        }
        let schema = Arc::new(arrow::datatypes::Schema::new(arrow_fields));

        Ok(Self {
            schema,
            row_count: 0,
            attr_defs,
            field_builders,
            attr_builders,
        })
    }

    pub fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    pub fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        // fixed fields
        let mut columns: Vec<ArrayRef> = self
            .field_builders
            .iter_mut()
            .map(|(_, builder)| builder.finish())
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
            columns.push(Arc::new(attrs));
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
                match attrs.get(def.name.as_bytes()) {
                    Some(result) => match result {
                        Ok(value) => {
                            let value = AttributeValue::from(&value);
                            builder.append_value(&value)?;
                        }
                        Err(_) => {
                            eprintln!("Error parsing attribute: {:?}", def.name);
                            builder.append_null();
                        }
                    },
                    None => {
                        builder.append_null();
                    }
                };
            }
        }
        self.row_count += 1;
        Ok(())
    }
}

/// Append a GTF record to the batch.
impl<'a> Push<&'a noodles::gtf::Record<'a>> for BatchBuilder {
    fn push(&mut self, record: &noodles::gtf::Record) -> io::Result<()> {
        use noodles::gff::feature::record::Attributes as FeatureAttributes;
        use noodles::gtf::record::Attributes as GtfAttributes;

        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }

        // attributes (optional)
        if !self.attr_defs.is_empty() {
            let attrs = record.attributes()?;
            for (def, builder) in self.attr_builders.iter_mut() {
                match <GtfAttributes as FeatureAttributes>::get(&attrs, def.name.as_bytes()) {
                    Some(result) => match result {
                        Ok(value) => {
                            let value = AttributeValue::from(value);
                            builder.append_value(&value)?;
                        }
                        Err(_) => {
                            eprintln!("Error parsing attribute: {:?}", def.name);
                            builder.append_null();
                        }
                    },
                    None => {
                        builder.append_null();
                    }
                };
            }
        }
        self.row_count += 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::datatypes::{DataType, Field as ArrowField};
    use noodles::gff::feature::RecordBuf;

    fn recordbuf_to_gff_line(record_buf: &RecordBuf) -> noodles::gff::Line {
        let mut writer = noodles::gff::io::Writer::new(Vec::new());
        writer.write_record(record_buf).unwrap();
        let buf = writer.into_inner();
        let mut reader = noodles::gff::io::Reader::new(std::io::Cursor::new(buf.as_slice()));
        let line = reader.lines().next().unwrap().unwrap();
        line
    }

    fn recordbuf_to_gtf_line(record_buf: &RecordBuf) -> noodles::gtf::Line {
        let mut writer = noodles::gtf::io::Writer::new(Vec::new());
        writer.write_record(record_buf).unwrap();
        let buf = writer.into_inner();
        let mut reader = noodles::gtf::io::Reader::new(std::io::Cursor::new(buf.as_slice()));
        let line = reader.lines().next().unwrap().unwrap();
        line
    }

    #[test]
    fn test_batch_builder_new() {
        let field_names = Some(vec!["seqid".to_string(), "source".to_string()]);
        let attr_defs = Some(vec![("gene_id".to_string(), "String".to_string())]);
        let capacity = 10;

        let batch_builder = BatchBuilder::new(field_names, attr_defs, capacity).unwrap();

        assert_eq!(batch_builder.attr_defs.len(), 1);
    }

    #[test]
    fn test_schema() {
        let field_names = Some(vec!["seqid".to_string(), "source".to_string()]);
        let attr_defs = Some(vec![("gene_id".to_string(), "String".to_string())]);
        let capacity = 10;

        let batch_builder = BatchBuilder::new(field_names, attr_defs, capacity).unwrap();
        let schema = batch_builder.schema();

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

        let record_buf = RecordBuf::default();
        let line = recordbuf_to_gff_line(&record_buf);
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

        let record_buf = RecordBuf::default();
        let line = recordbuf_to_gtf_line(&record_buf);
        let record = line.as_record().unwrap().unwrap();
        let result = batch_builder.push(&record);
        assert!(result.is_ok());
    }

    #[test]
    fn test_finish() {
        let field_names = Some(vec!["seqid".to_string(), "source".to_string()]);
        let attr_defs = Some(vec![("gene_id".to_string(), "String".to_string())]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(field_names, attr_defs, capacity).unwrap();

        let record_buf = RecordBuf::default();
        let line = recordbuf_to_gtf_line(&record_buf);
        let record = line.as_record().unwrap().unwrap();
        batch_builder.push(&record).unwrap();
        let record_batch = batch_builder.finish();
        assert!(record_batch.is_ok());

        let record_batch = record_batch.unwrap();
        assert_eq!(record_batch.num_columns(), 3);
        assert_eq!(record_batch.num_rows(), 1);
    }
}

use std::sync::Arc;

use arrow::array::{ArrayRef, StructArray};
use arrow::datatypes::FieldRef;
use arrow::datatypes::SchemaRef;
use arrow::error::ArrowError;
use arrow::record_batch::{RecordBatch, RecordBatchOptions};
use indexmap::IndexMap;

use crate::batch::{Push, RecordBatchBuilder};
use crate::gxf::model::attribute::{AttributeBuilder, AttributeDef, AttributeValue};
use crate::gxf::model::field::Push as _;
use crate::gxf::model::field::{Field, FieldBuilder};

use super::Model;

/// A builder for an Arrow record batch of GXF (GTF/GFF) features.
pub struct BatchBuilder {
    schema: SchemaRef,
    row_count: usize,
    has_attributes: bool,
    field_builders: IndexMap<Field, FieldBuilder>,
    attr_builders: IndexMap<AttributeDef, AttributeBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` for GTF/GFF records.
    ///
    /// - `fields`: standard GXF field names. `None` → all 8 standard fields.
    /// - `attr_defs`: `None` → no attributes column. `Some(vec![])` → empty struct.
    pub fn new(
        fields: Option<Vec<String>>,
        attr_defs: Option<Vec<(String, String)>>,
        capacity: usize,
    ) -> crate::Result<Self> {
        let model = Model::new(fields, attr_defs)?;
        Self::from_model(&model, capacity)
    }

    /// Creates a new `BatchBuilder` from a [`Model`].
    pub fn from_model(model: &Model, capacity: usize) -> crate::Result<Self> {
        let mut field_builders = IndexMap::new();
        for field in model.fields() {
            let builder = FieldBuilder::new(field.clone(), capacity);
            field_builders.insert(field.clone(), builder);
        }

        let mut attr_builders = IndexMap::new();
        if let Some(defs) = model.attr_defs() {
            for def in defs {
                let builder = AttributeBuilder::new(&def.ty);
                attr_builders.insert(def.clone(), builder);
            }
        }

        Ok(Self {
            schema: model.schema().clone(),
            row_count: 0,
            has_attributes: model.has_attributes(),
            field_builders,
            attr_builders,
        })
    }
}

impl RecordBatchBuilder for BatchBuilder {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        // fixed fields
        let mut columns: Vec<ArrayRef> = self
            .field_builders
            .iter_mut()
            .map(|(_, builder)| builder.finish())
            .collect();

        // attributes (optional)
        if self.has_attributes {
            if self.attr_builders.is_empty() {
                let attrs = StructArray::new_empty_fields(self.row_count, None);
                columns.push(Arc::new(attrs));
            } else {
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

/// Append a GFF record to the batch.
impl<'a> Push<&'a noodles::gff::Record<'a>> for BatchBuilder {
    fn push(&mut self, record: &noodles::gff::Record) -> crate::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }

        // attributes (optional)
        if self.has_attributes {
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
    fn push(&mut self, record: &noodles::gtf::Record) -> crate::Result<()> {
        use noodles::gff::feature::record::Attributes as FeatureAttributes;
        use noodles::gtf::record::Attributes as GtfAttributes;

        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }

        // attributes (optional)
        if self.has_attributes {
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

        assert!(batch_builder.has_attributes);
        assert_eq!(batch_builder.attr_builders.len(), 1);
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

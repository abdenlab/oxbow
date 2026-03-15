use std::sync::Arc;

use arrow::array::{ArrayRef, StructArray};
use arrow::datatypes::{FieldRef, SchemaRef};
use arrow::error::ArrowError;
use arrow::record_batch::{RecordBatch, RecordBatchOptions};
use indexmap::IndexMap;
use noodles::sam::alignment::record::data::field::Tag;

use crate::batch::{Push, RecordBatchBuilder};

use super::field::Push as _;
use super::field::{Field, FieldBuilder};
use super::tag::{TagBuilder, TagDef};
use super::Model;

/// A builder for an Arrow record batch of alignments.
pub struct BatchBuilder {
    schema: SchemaRef,
    row_count: usize,
    header: noodles::sam::Header,
    has_tags: bool,
    field_builders: IndexMap<Field, FieldBuilder>,
    tag_builders: IndexMap<TagDef, TagBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` for SAM/BAM records.
    ///
    /// - `fields`: standard SAM field names. `None` → all 12 standard fields.
    /// - `tag_defs`: `None` → no tags column. `Some(vec![])` → empty struct.
    pub fn new(
        header: noodles::sam::Header,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        capacity: usize,
    ) -> crate::Result<Self> {
        let model = Model::new(fields, tag_defs)?;
        Self::from_model(&model, header, capacity)
    }

    /// Creates a new `BatchBuilder` from a [`Model`].
    pub fn from_model(
        model: &Model,
        header: noodles::sam::Header,
        capacity: usize,
    ) -> crate::Result<Self> {
        let ref_names: Vec<String> = header
            .reference_sequences()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect();

        let mut field_builders = IndexMap::new();
        for field in model.fields() {
            let builder = match field {
                Field::Rname | Field::Rnext => {
                    FieldBuilder::with_refs(field.clone(), capacity, &ref_names)
                        .map_err(|e| crate::OxbowError::invalid_data(e.to_string()))?
                }
                _ => FieldBuilder::new(field.clone(), capacity),
            };
            field_builders.insert(field.clone(), builder);
        }

        let mut tag_builders = IndexMap::new();
        if let Some(defs) = model.tag_defs() {
            for tag in defs {
                let builder = TagBuilder::new(&tag.ty);
                tag_builders.insert(tag.clone(), builder);
            }
        }

        Ok(Self {
            schema: model.schema().clone(),
            row_count: 0,
            header,
            has_tags: model.has_tags(),
            field_builders,
            tag_builders,
        })
    }

    pub fn header(&self) -> noodles::sam::Header {
        self.header.clone()
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

        // tags (optional)
        if self.has_tags {
            if self.tag_builders.is_empty() {
                // tags column present but no tag defs → empty struct
                let tags = StructArray::new_empty_fields(self.row_count, None);
                columns.push(Arc::new(tags));
            } else {
                let tag_arrays: Vec<(FieldRef, ArrayRef)> = self
                    .tag_builders
                    .iter_mut()
                    .map(|(def, builder)| {
                        let arrow_field = def.get_arrow_field();
                        (Arc::new(arrow_field), builder.finish())
                    })
                    .collect();
                let tags = StructArray::from(tag_arrays);
                columns.push(Arc::new(tags));
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

/// Append a SAM record to the batch.
impl Push<&noodles::sam::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::sam::Record) -> crate::Result<()> {
        // fixed fields
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record, &self.header)?;
        }

        // tags (optional)
        if !self.tag_builders.is_empty() {
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
        self.row_count += 1;
        Ok(())
    }
}

/// Append a BAM record to the batch.
impl Push<&noodles::bam::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::bam::Record) -> crate::Result<()> {
        // fixed fields
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record, &self.header)?;
        }

        // tags (optional)
        if !self.tag_builders.is_empty() {
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
        self.row_count += 1;
        Ok(())
    }
}

/// Append a CRAM record to the batch.
impl Push<&noodles::sam::alignment::RecordBuf> for BatchBuilder {
    fn push(&mut self, record: &noodles::sam::alignment::RecordBuf) -> crate::Result<()> {
        use noodles::sam::alignment::record::data::field::Value;
        // fixed fields
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record, &self.header)?;
        }

        // tags (optional)
        if !self.tag_builders.is_empty() {
            let data = record.data();
            for (def, builder) in self.tag_builders.iter_mut() {
                let tag = Tag::from(def.into_bytes());
                match data.get(&tag) {
                    // tag is present in this record buffer
                    Some(value_buf) => {
                        let value = Value::from(value_buf);
                        builder.append_value(&value)?;
                    }
                    // tag is not present in this record buffer
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
    use arrow::datatypes::{DataType, Field as ArrowField};

    use super::*;

    #[test]
    fn test_batch_builder_new() {
        let header = noodles::sam::Header::default();
        let fields = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let tag_defs = Some(vec![("NM".to_string(), "i".to_string())]);
        let capacity = 10;

        let batch_builder = BatchBuilder::new(header.clone(), fields, tag_defs, capacity).unwrap();
        assert_eq!(batch_builder.header(), header);
        assert_eq!(batch_builder.tag_builders.len(), 1);
        assert!(batch_builder.has_tags);
    }

    #[test]
    fn test_schema() {
        let header = noodles::sam::Header::default();
        let fields = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let tag_defs = Some(vec![("NM".to_string(), "i".to_string())]);
        let capacity = 10;

        let batch_builder = BatchBuilder::new(header, fields, tag_defs, capacity).unwrap();
        let schema = batch_builder.schema();
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
    fn test_no_tags_when_tag_defs_none() {
        let header = noodles::sam::Header::default();
        let fields = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let capacity = 10;

        let batch_builder = BatchBuilder::new(header, fields, None, capacity).unwrap();
        assert!(!batch_builder.has_tags);
        assert!(batch_builder.tag_builders.is_empty());
        assert_eq!(batch_builder.schema().fields().len(), 2);
    }

    #[test]
    fn test_from_model() {
        let model = Model::new(
            Some(vec!["qname".into(), "pos".into()]),
            Some(vec![("NM".into(), "i".into())]),
        )
        .unwrap();
        let header = noodles::sam::Header::default();
        let batch_builder = BatchBuilder::from_model(&model, header, 10).unwrap();
        assert!(batch_builder.has_tags);
        assert_eq!(batch_builder.tag_builders.len(), 1);
        assert_eq!(batch_builder.schema().fields().len(), 3);
    }

    #[test]
    fn test_push_sam_record() {
        let header = noodles::sam::Header::default();
        let fields = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let tag_defs = Some(vec![("NM".to_string(), "i".to_string())]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(header, fields, tag_defs, capacity).unwrap();

        let record = noodles::sam::Record::default();
        let result = batch_builder.push(&record);
        assert!(result.is_ok());
    }

    #[test]
    fn test_push_bam_record() {
        let header = noodles::sam::Header::default();
        let fields = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let tag_defs = Some(vec![("NM".to_string(), "i".to_string())]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(header, fields, tag_defs, capacity).unwrap();

        let record = noodles::bam::Record::default();
        let result = batch_builder.push(&record);
        assert!(result.is_ok());
    }

    #[test]
    fn test_finish() {
        let header = noodles::sam::Header::default();
        let fields = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let tag_defs = Some(vec![("NM".to_string(), "i".to_string())]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(header, fields, tag_defs, capacity).unwrap();

        let record = noodles::sam::Record::default();
        batch_builder.push(&record).unwrap();
        let record_batch = batch_builder.finish();
        assert!(record_batch.is_ok());

        let record_batch = record_batch.unwrap();
        assert_eq!(record_batch.num_columns(), 3);
        assert_eq!(record_batch.num_rows(), 1);
    }

    #[test]
    fn test_finish_empty_tags() {
        let header = noodles::sam::Header::default();
        let fields = Some(vec!["QNAME".to_string()]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(header, fields, Some(vec![]), capacity).unwrap();

        let record = noodles::sam::Record::default();
        batch_builder.push(&record).unwrap();
        let record_batch = batch_builder.finish().unwrap();
        assert_eq!(record_batch.num_columns(), 2);
        assert_eq!(record_batch.num_rows(), 1);
    }

    #[test]
    fn test_finish_no_tags() {
        let header = noodles::sam::Header::default();
        let fields = Some(vec!["QNAME".to_string(), "FLAG".to_string()]);
        let capacity = 10;

        let mut batch_builder = BatchBuilder::new(header, fields, None, capacity).unwrap();

        let record = noodles::sam::Record::default();
        batch_builder.push(&record).unwrap();
        let record_batch = batch_builder.finish().unwrap();
        assert_eq!(record_batch.num_columns(), 2);
        assert_eq!(record_batch.num_rows(), 1);
    }
}

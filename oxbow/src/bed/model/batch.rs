use std::str::FromStr;
use std::sync::Arc;

use arrow::array::ArrayRef;
use arrow::datatypes::SchemaRef;
use arrow::error::ArrowError;
use arrow::record_batch::{RecordBatch, RecordBatchOptions};
use indexmap::IndexMap;

use crate::batch::{Push, RecordBatchBuilder};

use super::field::Push as _;
use super::field::{Field, FieldBuilder};
use super::field_def::{FieldBuilder as GenericFieldBuilder, FieldDef, Push as GenericPush};
use super::schema::BedSchema;

/// A builder for an Arrow record batch of BED features.
pub struct BatchBuilder {
    schema: SchemaRef,
    row_count: usize,
    bed_schema: BedSchema,
    standard_field_builders: IndexMap<Field, FieldBuilder>,
    custom_field_builders: IndexMap<FieldDef, GenericFieldBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` from a [`Model`].
    pub fn from_model(model: &super::Model, capacity: usize) -> crate::Result<Self> {
        Self::new(model.bed_schema(), Some(model.field_names()), capacity)
    }

    /// Creates a new `BatchBuilder` for BED records.
    pub fn new(
        bed_schema: &BedSchema,
        field_names: Option<Vec<String>>,
        capacity: usize,
    ) -> crate::Result<Self> {
        // All the standard and custom field definitions for the given BED schema.
        let standard_field_names = bed_schema.standard_field_names();
        let custom_field_defs = bed_schema.custom_fields();

        // Determine the projected fields and create builders.
        let mut projected_standard_fields = Vec::new();
        let mut projected_custom_defs = Vec::new();
        let mut standard_field_builders = IndexMap::new();
        let mut custom_field_builders = IndexMap::new();
        match &field_names {
            Some(projection) => {
                for name in projection {
                    if let Some(def) = custom_field_defs.iter().find(|d| &d.name == name) {
                        projected_custom_defs.push(def.clone());
                        let builder = GenericFieldBuilder::new(&def.ty, capacity)?;
                        custom_field_builders.insert(def.clone(), builder);
                    } else {
                        let field = Field::from_str(name)?;
                        projected_standard_fields.push(field.clone());
                        let builder = FieldBuilder::new(field.clone(), capacity);
                        standard_field_builders.insert(field.clone(), builder);
                    }
                }
            }
            None => {
                projected_standard_fields.extend(
                    standard_field_names
                        .into_iter()
                        .map(|name| Field::from_str(&name))
                        .collect::<Result<Vec<_>, _>>()?,
                );
                for field in &projected_standard_fields {
                    let builder = FieldBuilder::new(field.clone(), capacity);
                    standard_field_builders.insert(field.clone(), builder);
                }
                for def in custom_field_defs {
                    projected_custom_defs.push(def.clone());
                    let builder = GenericFieldBuilder::new(&def.ty, capacity)?;
                    custom_field_builders.insert(def.clone(), builder);
                }
            }
        };

        // Build schema once
        let mut arrow_fields: Vec<arrow::datatypes::Field> = projected_standard_fields
            .iter()
            .map(|field| field.get_arrow_field())
            .collect();
        for def in &projected_custom_defs {
            arrow_fields.push(def.get_arrow_field());
        }
        let schema = Arc::new(arrow::datatypes::Schema::new(arrow_fields));

        Ok(Self {
            schema,
            row_count: 0,
            bed_schema: bed_schema.clone(),
            standard_field_builders,
            custom_field_builders,
        })
    }

    pub fn bed_schema(&self) -> &BedSchema {
        &self.bed_schema
    }
}

impl RecordBatchBuilder for BatchBuilder {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        // standard fields
        let mut columns: Vec<ArrayRef> = self
            .standard_field_builders
            .iter_mut()
            .map(|(_, builder)| builder.finish())
            .collect();

        // custom fields (optional)
        if !self.custom_field_builders.is_empty() {
            let other: Vec<ArrayRef> = self
                .custom_field_builders
                .iter_mut()
                .map(|(_, builder)| builder.finish())
                .collect();
            columns.extend(other);
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

/// Append a BED record to the batch.
impl Push<&noodles::bed::Record<3>> for BatchBuilder {
    fn push(&mut self, record: &noodles::bed::Record<3>) -> crate::Result<()> {
        // standard fields
        for (_, builder) in self.standard_field_builders.iter_mut() {
            builder.push(record)?;
        }

        // custom fields
        let n = self.bed_schema.standard_field_count();
        if !self.custom_field_builders.is_empty() {
            if self.bed_schema.custom_field_count().is_none() {
                // BEDn+ mode: lump remaining fields into "rest"
                let rest = record
                    .other_fields()
                    .iter()
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>()
                    .join("\t");
                let rest_def =
                    FieldDef::new("rest".to_string(), super::field_def::FieldType::String);
                if let Some(builder) = self.custom_field_builders.get_mut(&rest_def) {
                    builder.push(&rest)?;
                };
            } else {
                // BEDn+m mode: parse each custom field by its FieldDef.
                // Use the full custom field list for positional alignment,
                // only pushing to builders that exist (i.e., projected fields).
                let all_custom_defs = self.bed_schema.custom_fields();
                for (i, value) in record.other_fields().iter().enumerate() {
                    let field_abs_idx = i + 3;
                    let field_rel_idx = match field_abs_idx.checked_sub(n) {
                        Some(i) => i,
                        None => continue,
                    };
                    if let Some(def) = all_custom_defs.get(field_rel_idx) {
                        if let Some(builder) = self.custom_field_builders.get_mut(def) {
                            builder.push(&value.to_string())?;
                        }
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
    use arrow::array::{Int64Array, StringArray};
    use std::io::Cursor;

    fn create_bed_record() -> noodles::bed::Record<3> {
        let mut record = noodles::bed::Record::default();
        let buf = b"chr1\t100\t200\tfoo\tbar\n".to_vec();
        let mut reader = noodles::bed::io::Reader::<3, Cursor<Vec<u8>>>::new(Cursor::new(buf));
        reader.read_record(&mut record).unwrap();
        record
    }

    #[test]
    fn test_batch_builder_new() {
        let bed_schema = BedSchema::new_from_nm(3, Some(2)).unwrap();
        let batch_builder = BatchBuilder::new(&bed_schema, None, 10).unwrap();

        assert_eq!(batch_builder.schema().fields().len(), 5);
    }

    #[test]
    fn test_push_bed_record() {
        let bed_schema = BedSchema::new_from_nm(3, Some(2)).unwrap();
        let mut batch_builder = BatchBuilder::new(&bed_schema, None, 10).unwrap();

        let record = create_bed_record();
        let result = batch_builder.push(&record);
        assert!(result.is_ok());
    }

    #[test]
    fn test_finish_bedn() {
        let record = create_bed_record();

        let bed_schema = BedSchema::new_from_nm(3, Some(0)).unwrap();
        let mut batch_builder = BatchBuilder::new(&bed_schema, None, 10).unwrap();
        batch_builder.push(&record).unwrap();
        let record_batch = batch_builder.finish();
        assert!(record_batch.is_ok());

        let record_batch = record_batch.unwrap();
        assert_eq!(record_batch.num_columns(), 3);
        assert_eq!(record_batch.num_rows(), 1);

        let chrom_array = record_batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(chrom_array.value(0), "chr1");

        let start_array = record_batch
            .column(1)
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        assert_eq!(start_array.value(0), 101);

        let end_array = record_batch
            .column(2)
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        assert_eq!(end_array.value(0), 200);
    }

    #[test]
    fn test_finish_bedn_plus_m() {
        let bed_schema = BedSchema::new_from_nm(3, Some(2)).unwrap();
        let mut batch_builder = BatchBuilder::new(&bed_schema, None, 10).unwrap();

        let record = create_bed_record();
        batch_builder.push(&record).unwrap();
        let record_batch = batch_builder.finish();
        assert!(record_batch.is_ok());

        let record_batch = record_batch.unwrap();
        assert_eq!(record_batch.num_columns(), 5);
        assert_eq!(record_batch.num_rows(), 1);

        let chrom_array = record_batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(chrom_array.value(0), "chr1");

        let start_array = record_batch
            .column(1)
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        assert_eq!(start_array.value(0), 101);

        let end_array = record_batch
            .column(2)
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        assert_eq!(end_array.value(0), 200);

        let custom1_array = record_batch
            .column(3)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(custom1_array.value(0), "foo");
        assert_eq!(
            record_batch.schema().fields().get(3).unwrap().name(),
            "BED3+1"
        );

        let custom2_array = record_batch
            .column(4)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(custom2_array.value(0), "bar");
        assert_eq!(
            record_batch.schema().fields().get(4).unwrap().name(),
            "BED3+2"
        );
    }

    #[test]
    fn test_finish_bedn_plus() {
        let record = create_bed_record();

        let bed_schema = BedSchema::new_from_nm(3, None).unwrap();
        let mut batch_builder = BatchBuilder::new(&bed_schema, None, 10).unwrap();
        batch_builder.push(&record).unwrap();
        let record_batch = batch_builder.finish();
        assert!(record_batch.is_ok());

        let record_batch = record_batch.unwrap();
        assert_eq!(record_batch.num_columns(), 4);
        assert_eq!(record_batch.num_rows(), 1);

        let chrom_array = record_batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(chrom_array.value(0), "chr1");

        let start_array = record_batch
            .column(1)
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        assert_eq!(start_array.value(0), 101);

        let end_array = record_batch
            .column(2)
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap();
        assert_eq!(end_array.value(0), 200);

        let rest_array = record_batch
            .column(3)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(rest_array.value(0), "foo\tbar");
        assert_eq!(
            record_batch.schema().fields().get(3).unwrap().name(),
            "rest"
        );
    }
}

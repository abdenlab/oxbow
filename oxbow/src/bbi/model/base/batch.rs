use std::iter::zip;

use arrow::array::ArrayRef;
use arrow::datatypes::SchemaRef;
use arrow::error::ArrowError;
use arrow::record_batch::{RecordBatch, RecordBatchOptions};
use indexmap::IndexMap;

use crate::batch::{Push, RecordBatchBuilder};

use super::field::Push as _;
pub use super::field::{FieldBuilder, FieldDef, FieldType};
pub use super::{BedSchema, BigBedRecord, BigWigRecord, Model};

/// A builder for an Arrow record batch of BBI records defined by AutoSql.
pub struct BatchBuilder {
    schema: SchemaRef,
    row_count: usize,
    bed_schema: BedSchema,
    bed_schema_field_defs: Vec<FieldDef>,
    builders: IndexMap<FieldDef, FieldBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` for BigWig or BigBed records.
    pub fn new(
        bed_schema: BedSchema,
        fields: Option<Vec<String>>,
        capacity: usize,
    ) -> crate::Result<Self> {
        let model = Model::new(bed_schema, fields)?;
        Self::from_model(&model, capacity)
    }

    /// Creates a new `BatchBuilder` from a [`Model`].
    pub fn from_model(model: &Model, capacity: usize) -> crate::Result<Self> {
        let mut builders = IndexMap::new();
        for def in model.fields() {
            let builder = FieldBuilder::new(&def.ty, capacity)?;
            builders.insert(def.clone(), builder);
        }

        Ok(Self {
            schema: model.schema().clone(),
            row_count: 0,
            bed_schema: model.bed_schema().clone(),
            bed_schema_field_defs: model.bed_schema_field_defs(),
            builders,
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
        let columns: Vec<ArrayRef> = self
            .builders
            .iter_mut()
            .map(|(_, builder)| builder.finish())
            .collect();
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

/// Append a BigBed record to the batch.
impl Push<&BigBedRecord<'_>> for BatchBuilder {
    fn push(&mut self, record: &BigBedRecord) -> crate::Result<()> {
        let mut schema_defs = self.bed_schema_field_defs.iter();

        if let Some(def) = schema_defs.next() {
            if let Some(builder) = self.builders.get_mut(def) {
                match builder {
                    FieldBuilder::String(b) => {
                        b.append_value(record.chrom);
                    }
                    _ => {
                        return Err(crate::OxbowError::invalid_data(
                            "Wrong builder type for chrom",
                        ))
                    }
                }
            }
        }

        if let Some(def) = schema_defs.next() {
            if let Some(builder) = self.builders.get_mut(def) {
                match builder {
                    FieldBuilder::Uint(b) => {
                        b.append_value(record.start);
                    }
                    _ => {
                        return Err(crate::OxbowError::invalid_data(
                            "Wrong builder type for start",
                        ))
                    }
                }
            }
        }

        if let Some(def) = schema_defs.next() {
            if let Some(builder) = self.builders.get_mut(def) {
                match builder {
                    FieldBuilder::Uint(b) => {
                        b.append_value(record.end);
                    }
                    _ => {
                        return Err(crate::OxbowError::invalid_data(
                            "Wrong builder type for end",
                        ))
                    }
                }
            }
        }

        // Parse remaining fields from the tab-separated rest string.
        // schema_defs still has defs for fields 4+ (standard and custom).
        if self.bed_schema.custom_field_count().is_none() {
            // BEDn+ mode: lump all rest into a single "rest" field
            let rest_def = FieldDef::new("rest".to_string(), FieldType::String);
            if let Some(builder) = self.builders.get_mut(&rest_def) {
                match builder {
                    FieldBuilder::String(b) => {
                        b.append_value(record.rest);
                    }
                    _ => {
                        return Err(crate::OxbowError::invalid_data(
                            "Wrong builder type for rest",
                        ))
                    }
                }
            };
        } else if !record.rest.is_empty() {
            // BEDn or BEDn+m: parse each field positionally
            let rest = record.rest.split('\t');
            for (def, value) in zip(schema_defs, rest) {
                if let Some(builder) = self.builders.get_mut(def) {
                    builder.push(value)?;
                }
            }
        }

        self.row_count += 1;
        Ok(())
    }
}

/// Append a BigWig record to the batch.
impl Push<&BigWigRecord<'_>> for BatchBuilder {
    fn push(&mut self, record: &BigWigRecord) -> crate::Result<()> {
        let mut schema_defs = self.bed_schema_field_defs.iter();

        if let Some(def) = schema_defs.next() {
            if let Some(builder) = self.builders.get_mut(def) {
                match builder {
                    FieldBuilder::String(b) => {
                        b.append_value(record.chrom);
                    }
                    _ => {
                        return Err(crate::OxbowError::invalid_data(
                            "Wrong builder type for chrom",
                        ))
                    }
                }
            }
        }

        if let Some(def) = schema_defs.next() {
            if let Some(builder) = self.builders.get_mut(def) {
                match builder {
                    FieldBuilder::Uint(b) => {
                        b.append_value(record.start);
                    }
                    _ => {
                        return Err(crate::OxbowError::invalid_data(
                            "Wrong builder type for start",
                        ))
                    }
                }
            }
        }

        if let Some(def) = schema_defs.next() {
            if let Some(builder) = self.builders.get_mut(def) {
                match builder {
                    FieldBuilder::Uint(b) => {
                        b.append_value(record.end);
                    }
                    _ => {
                        return Err(crate::OxbowError::invalid_data(
                            "Wrong builder type for end",
                        ))
                    }
                }
            }
        }

        if let Some(def) = schema_defs.next() {
            if let Some(builder) = self.builders.get_mut(def) {
                match builder {
                    FieldBuilder::Float(b) => {
                        b.append_value(record.value);
                    }
                    _ => {
                        return Err(crate::OxbowError::invalid_data(
                            "Wrong builder type for value",
                        ))
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
    use arrow::datatypes::DataType;

    fn create_test_bedschema() -> BedSchema {
        let fields = vec![FieldDef::new("value".to_string(), FieldType::Float)];
        BedSchema::new(3, Some(fields)).unwrap()
    }

    #[test]
    fn test_batch_builder_new() {
        let bed_schema = create_test_bedschema();
        let model = Model::new(bed_schema, None).unwrap();
        let builder = BatchBuilder::from_model(&model, 10).unwrap();

        assert_eq!(builder.schema().fields().len(), 4);
        assert_eq!(builder.builders.len(), 4);
    }

    #[test]
    fn test_schema() {
        let bed_schema = create_test_bedschema();
        let model = Model::new(bed_schema, None).unwrap();
        let builder = BatchBuilder::from_model(&model, 10).unwrap();

        let schema = builder.schema();
        assert_eq!(schema.fields().len(), 4);
        assert_eq!(schema.field(0).name(), "chrom");
        assert_eq!(schema.field(1).name(), "start");
        assert_eq!(schema.field(2).name(), "end");
        assert_eq!(schema.field(3).name(), "value");
        assert_eq!(*schema.field(3).data_type(), DataType::Float32);
    }

    #[test]
    fn test_push_bigbed_record() {
        let schema = create_test_bedschema();
        let model = Model::new(schema, None).unwrap();
        let mut builder = BatchBuilder::from_model(&model, 10).unwrap();

        let record = BigBedRecord {
            chrom: "chr1",
            start: 100,
            end: 200,
            rest: &"extra_field".to_string(),
        };
        builder.push(&record).unwrap();
        let batch = builder.finish().unwrap();

        assert_eq!(batch.num_rows(), 1);
        assert_eq!(
            batch
                .column(0)
                .as_any()
                .downcast_ref::<arrow::array::StringArray>()
                .unwrap()
                .value(0),
            "chr1"
        );
        assert_eq!(
            batch
                .column(1)
                .as_any()
                .downcast_ref::<arrow::array::UInt32Array>()
                .unwrap()
                .value(0),
            100
        );
        assert_eq!(
            batch
                .column(2)
                .as_any()
                .downcast_ref::<arrow::array::UInt32Array>()
                .unwrap()
                .value(0),
            200
        );
    }

    #[test]
    fn test_push_bigwig_record() {
        let schema = create_test_bedschema();
        let model = Model::new(schema, None).unwrap();
        let mut builder = BatchBuilder::from_model(&model, 10).unwrap();

        let record = BigWigRecord {
            chrom: "chr1",
            start: 100,
            end: 200,
            value: 1.23,
        };
        builder.push(&record).unwrap();
        let batch = builder.finish().unwrap();

        assert_eq!(batch.num_rows(), 1);
        assert_eq!(
            batch
                .column(0)
                .as_any()
                .downcast_ref::<arrow::array::StringArray>()
                .unwrap()
                .value(0),
            "chr1"
        );
        assert_eq!(
            batch
                .column(1)
                .as_any()
                .downcast_ref::<arrow::array::UInt32Array>()
                .unwrap()
                .value(0),
            100
        );
        assert_eq!(
            batch
                .column(2)
                .as_any()
                .downcast_ref::<arrow::array::UInt32Array>()
                .unwrap()
                .value(0),
            200
        );
        assert_eq!(
            batch
                .column(3)
                .as_any()
                .downcast_ref::<arrow::array::Float32Array>()
                .unwrap()
                .value(0),
            1.23
        );
    }

    #[test]
    fn test_finish_empty_batch() {
        let schema = create_test_bedschema();
        let model = Model::new(schema, None).unwrap();
        let mut builder = BatchBuilder::from_model(&model, 10).unwrap();

        let batch = builder.finish().unwrap();
        assert_eq!(batch.num_rows(), 0);
        assert_eq!(batch.num_columns(), 4);
    }

    #[test]
    fn test_bigbed_bed6_no_custom() {
        // bed6 with no custom fields — standard fields 4-6 are in rest
        let bed_schema: BedSchema = "bed6".parse().unwrap();
        let model = Model::new(bed_schema, None).unwrap();
        let mut builder = BatchBuilder::from_model(&model, 10).unwrap();

        let record = BigBedRecord {
            chrom: "chr1",
            start: 100,
            end: 200,
            rest: &"gene1\t500\t+".to_string(),
        };
        builder.push(&record).unwrap();
        let batch = builder.finish().unwrap();

        assert_eq!(batch.num_rows(), 1);
        assert_eq!(batch.num_columns(), 6);
        // name (field 4) should be parsed from rest
        assert_eq!(
            batch
                .column(3)
                .as_any()
                .downcast_ref::<arrow::array::StringArray>()
                .unwrap()
                .value(0),
            "gene1"
        );
        // score (field 5)
        assert_eq!(
            batch
                .column(4)
                .as_any()
                .downcast_ref::<arrow::array::UInt16Array>()
                .unwrap()
                .value(0),
            500
        );
    }

    #[test]
    fn test_bigbed_bed6_projected() {
        // bed6 projected to chrom + strand — strand is field 6, skipping name/score
        let bed_schema: BedSchema = "bed6".parse().unwrap();
        let model = Model::new(bed_schema, Some(vec!["chrom".into(), "strand".into()])).unwrap();
        let mut builder = BatchBuilder::from_model(&model, 10).unwrap();

        let record = BigBedRecord {
            chrom: "chr1",
            start: 100,
            end: 200,
            rest: &"gene1\t500\t+".to_string(),
        };
        builder.push(&record).unwrap();
        let batch = builder.finish().unwrap();

        assert_eq!(batch.num_rows(), 1);
        assert_eq!(batch.num_columns(), 2);
        assert_eq!(
            batch
                .column(1)
                .as_any()
                .downcast_ref::<arrow::array::StringArray>()
                .unwrap()
                .value(0),
            "+"
        );
    }
}

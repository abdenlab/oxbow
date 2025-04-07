use std::io;
use std::iter::zip;

use arrow::array::ArrayRef;
use arrow::datatypes::{Field as ArrowField, Schema as ArrowSchema};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use indexmap::IndexMap;

use super::field::Push as _;
pub use super::field::{FieldBuilder, FieldDef, FieldType};
pub use super::{BedSchema, BigBedRecord, BigWigRecord};

/// A builder for an Arrow record batch of BBI records defined by AutoSql.
pub struct BatchBuilder {
    schema: BedSchema,
    builders: IndexMap<FieldDef, FieldBuilder>,
}

impl BatchBuilder {
    // Creates a new `BatchBuilder` for BigWig or BigBed records.
    pub fn new(
        schema: BedSchema,
        field_names: Option<Vec<String>>,
        capacity: usize,
    ) -> io::Result<Self> {
        let schema_field_names = schema
            .fields()
            .iter()
            .map(|def| def.name.clone())
            .collect::<Vec<String>>();

        let field_names: Vec<String> = field_names.unwrap_or(schema_field_names);
        let field_defs: Vec<FieldDef> = field_names
            .iter()
            .map(|name| {
                schema
                    .fields()
                    .iter()
                    .find(|field| &field.name == name)
                    .cloned()
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Field '{}' not found in schema", name),
                        )
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut builders = IndexMap::new();
        for def in &field_defs {
            let builder = FieldBuilder::new(&def.ty, capacity)?;
            builders.insert(def.clone(), builder);
        }
        Ok(Self { schema, builders })
    }

    pub fn get_arrow_fields(&self) -> Vec<ArrowField> {
        let arrow_fields: Vec<ArrowField> = self
            .builders
            .iter()
            .map(|(def, _)| def.get_arrow_field())
            .collect();
        arrow_fields
    }

    pub fn get_arrow_schema(&self) -> ArrowSchema {
        ArrowSchema::new(self.get_arrow_fields())
    }

    pub fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        let name_to_array: Vec<(&str, ArrayRef)> = self
            .builders
            .iter_mut()
            .map(|(def, builder)| (def.name.as_str(), builder.finish()))
            .collect();
        RecordBatch::try_from_iter(name_to_array)
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Append a BigBed record to the batch.
impl Push<&BigBedRecord<'_>> for BatchBuilder {
    fn push(&mut self, record: &BigBedRecord) -> io::Result<()> {
        let mut schema_defs = self.schema.fields().iter();

        if let Some(def) = schema_defs.next() {
            if let Some(builder) = self.builders.get_mut(def) {
                match builder {
                    FieldBuilder::String(b) => {
                        b.append_value(record.chrom);
                    }
                    _ => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
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
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
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
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Wrong builder type for end",
                        ))
                    }
                }
            }
        }

        match self.schema.custom_field_count() {
            Some(0) => {}
            Some(_) => {
                let rest = record.rest.split('\t');
                for (def, value) in zip(schema_defs, rest) {
                    if let Some(builder) = self.builders.get_mut(def) {
                        builder.push(value)?;
                    }
                }
            }
            None => {
                let rest_def = FieldDef::new("rest".to_string(), FieldType::String);
                if let Some(builder) = self.builders.get_mut(&rest_def) {
                    match builder {
                        FieldBuilder::String(b) => {
                            b.append_value(record.rest);
                        }
                        _ => {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "Wrong builder type for rest",
                            ))
                        }
                    }
                };
            }
        }

        Ok(())
    }
}

/// Append a BigWig record to the batch.
impl Push<&BigWigRecord<'_>> for BatchBuilder {
    fn push(&mut self, record: &BigWigRecord) -> io::Result<()> {
        let mut schema_defs = self.schema.fields().iter();

        if let Some(def) = schema_defs.next() {
            if let Some(builder) = self.builders.get_mut(def) {
                match builder {
                    FieldBuilder::String(b) => {
                        b.append_value(record.chrom);
                    }
                    _ => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
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
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
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
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
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
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Wrong builder type for value",
                        ))
                    }
                }
            }
        }

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
        let builder = BatchBuilder::new(bed_schema.clone(), None, 10).unwrap();

        assert_eq!(builder.schema, bed_schema);
        assert_eq!(builder.builders.len(), bed_schema.fields().len());
    }

    #[test]
    fn test_get_arrow_schema() {
        let bed_schema = create_test_bedschema();
        let builder = BatchBuilder::new(bed_schema.clone(), None, 10).unwrap();

        let schema = builder.get_arrow_schema();
        assert_eq!(bed_schema.fields().len(), schema.fields().len());
        assert_eq!(schema.field(0).name(), "chrom");
        assert_eq!(schema.field(1).name(), "start");
        assert_eq!(schema.field(2).name(), "end");
        assert_eq!(schema.field(3).name(), "value");
        assert_eq!(*schema.field(3).data_type(), DataType::Float32);
    }

    #[test]
    fn test_push_bigbed_record() {
        let schema = create_test_bedschema();
        let mut builder = BatchBuilder::new(schema, None, 10).unwrap();

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
        let mut builder = BatchBuilder::new(schema, None, 10).unwrap();

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
        let mut builder = BatchBuilder::new(schema, None, 10).unwrap();

        let batch = builder.finish().unwrap();
        assert_eq!(batch.num_rows(), 0);
        assert_eq!(batch.num_columns(), 4);
    }
}

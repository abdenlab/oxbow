use std::io;

use arrow::array::ArrayRef;
use arrow::datatypes::{Field as ArrowField, Schema};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use indexmap::IndexMap;

use super::field::{Field, FieldBuilder, DEFAULT_FIELD_NAMES};
use super::BBIZoomRecord;

/// A builder for an Arrow record batch of BBI zoom level summary statistics.
pub struct BatchBuilder {
    fields: Vec<Field>,
    field_builders: IndexMap<Field, FieldBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` for BBI zoom level summary statistics records.
    pub fn new(
        ref_names: &[String],
        field_names: Option<Vec<String>>,
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
            let builder = FieldBuilder::new(field.clone(), ref_names, capacity)
                .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            field_builders.insert(field.clone(), builder);
        }

        Ok(Self {
            fields,
            field_builders,
        })
    }

    pub fn get_arrow_fields(&self) -> Vec<ArrowField> {
        self.fields
            .iter()
            .map(|field| field.get_arrow_field())
            .collect()
    }

    pub fn get_arrow_schema(&self) -> Schema {
        Schema::new(self.get_arrow_fields())
    }

    pub fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        let name_to_array: Vec<(&str, ArrayRef)> = self
            .fields
            .iter()
            .map(|field| {
                let builder = self.field_builders.get_mut(field).unwrap();
                let name = field.name();
                (name, builder.finish())
            })
            .collect();
        RecordBatch::try_from_iter(name_to_array)
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Append a BBIZoomRecord to the batch.
impl Push<&BBIZoomRecord<'_>> for BatchBuilder {
    fn push(&mut self, record: &BBIZoomRecord) -> io::Result<()> {
        for (_, builder) in &mut self.field_builders {
            match builder {
                FieldBuilder::Chrom(builder) => builder.append_value(record.chrom),
                FieldBuilder::Start(builder) => builder.append_value(record.start),
                FieldBuilder::End(builder) => builder.append_value(record.end),
                FieldBuilder::BasesCovered(builder) => builder.append_value(record.bases_covered),
                FieldBuilder::Min(builder) => builder.append_value(record.min),
                FieldBuilder::Max(builder) => builder.append_value(record.max),
                FieldBuilder::Sum(builder) => builder.append_value(record.sum),
                FieldBuilder::SumSquares(builder) => builder.append_value(record.sum_squares),
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{DictionaryArray, StringArray, UInt32Array};
    use arrow::datatypes::Int32Type;

    #[test]
    fn test_batch_builder_new() {
        let ref_names = vec!["chr1".to_string(), "chr2".to_string()];
        let field_names = Some(vec![
            "chrom".to_string(),
            "start".to_string(),
            "end".to_string(),
        ]);
        let capacity = 10;

        let builder = BatchBuilder::new(&ref_names, field_names, capacity);
        assert!(builder.is_ok());

        let builder = builder.unwrap();
        assert_eq!(builder.fields.len(), 3);
        assert_eq!(builder.field_builders.len(), 3);
    }

    #[test]
    fn test_get_arrow_schema() {
        let ref_names = vec!["chr1".to_string(), "chr2".to_string()];
        let field_names = Some(vec![
            "chrom".to_string(),
            "start".to_string(),
            "end".to_string(),
        ]);
        let capacity = 10;

        let builder = BatchBuilder::new(&ref_names, field_names, capacity).unwrap();
        let schema = builder.get_arrow_schema();

        assert_eq!(schema.fields().len(), 3);
        assert_eq!(schema.field(0).name(), "chrom");
        assert_eq!(schema.field(1).name(), "start");
        assert_eq!(schema.field(2).name(), "end");
    }

    #[test]
    fn test_push_zoom_records() {
        let ref_names = vec!["chr1".to_string(), "chr2".to_string()];
        let field_names = Some(vec![
            "chrom".to_string(),
            "start".to_string(),
            "end".to_string(),
        ]);
        let capacity = 10;

        let mut builder = BatchBuilder::new(&ref_names, field_names, capacity).unwrap();

        let record = BBIZoomRecord {
            chrom: "chr1",
            start: 100,
            end: 200,
            bases_covered: 50,
            min: 1.0,
            max: 5.0,
            sum: 10.0,
            sum_squares: 30.0,
        };
        let result = builder.push(&record);
        assert!(result.is_ok());
    }

    #[test]
    fn test_finish() {
        let ref_names = vec!["chr1".to_string(), "chr2".to_string()];
        let field_names = Some(vec![
            "chrom".to_string(),
            "start".to_string(),
            "end".to_string(),
        ]);
        let capacity = 10;

        let mut builder = BatchBuilder::new(&ref_names, field_names, capacity).unwrap();

        let record1 = BBIZoomRecord {
            chrom: "chr1",
            start: 100,
            end: 200,
            bases_covered: 50,
            min: 1.0,
            max: 5.0,
            sum: 10.0,
            sum_squares: 30.0,
        };

        let record2 = BBIZoomRecord {
            chrom: "chr2",
            start: 300,
            end: 400,
            bases_covered: 60,
            min: 2.0,
            max: 6.0,
            sum: 20.0,
            sum_squares: 40.0,
        };
        builder.push(&record1).unwrap();
        builder.push(&record2).unwrap();

        let batch = builder.finish();
        assert!(batch.is_ok());
        let batch = batch.unwrap();
        assert_eq!(batch.num_columns(), 3);
        assert_eq!(batch.num_rows(), 2);

        let chrom_array = batch
            .column(0)
            .as_any()
            .downcast_ref::<DictionaryArray<Int32Type>>()
            .unwrap()
            .values()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(chrom_array.value(0), "chr1");
        assert_eq!(chrom_array.value(1), "chr2");

        let start_array = batch
            .column(1)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        assert_eq!(start_array.value(0), 100);
        assert_eq!(start_array.value(1), 300);

        let end_array = batch
            .column(2)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        assert_eq!(end_array.value(0), 200);
        assert_eq!(end_array.value(1), 400);
    }
}

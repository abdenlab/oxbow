use arrow::array::ArrayRef;
use arrow::datatypes::SchemaRef;
use arrow::error::ArrowError;
use arrow::record_batch::{RecordBatch, RecordBatchOptions};
use indexmap::IndexMap;

use crate::batch::{Push, RecordBatchBuilder};
use crate::Select;

use super::field::Push as _;
use super::field::{Field, FieldBuilder};
use super::Model;

/// A builder for Arrow record batches of sequence records.
pub struct BatchBuilder {
    schema: SchemaRef,
    row_count: usize,
    field_builders: IndexMap<Field, FieldBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` for FASTQ records.
    pub fn new_fastq(fields: Select<String>, capacity: usize) -> crate::Result<Self> {
        let model = Model::new_fastq(fields)?;
        Self::from_model(&model, capacity)
    }

    /// Creates a new `BatchBuilder` for FASTA records.
    pub fn new_fasta(fields: Select<String>, capacity: usize) -> crate::Result<Self> {
        let model = Model::new_fasta(fields)?;
        Self::from_model(&model, capacity)
    }

    /// Creates a new `BatchBuilder` from a [`Model`].
    pub fn from_model(model: &Model, capacity: usize) -> crate::Result<Self> {
        let mut field_builders = IndexMap::new();
        for field in model.fields() {
            let builder = FieldBuilder::new(field.clone(), capacity);
            field_builders.insert(field.clone(), builder);
        }

        Ok(Self {
            schema: model.schema().clone(),
            row_count: 0,
            field_builders,
        })
    }
}

impl RecordBatchBuilder for BatchBuilder {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        let columns: Vec<ArrayRef> = self
            .field_builders
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

/// Append a FASTA record to the batch.
impl Push<&noodles::fasta::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::fasta::Record) -> crate::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }
        self.row_count += 1;
        Ok(())
    }
}

/// Append a FASTQ record to the batch.
impl Push<&noodles::fastq::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::fastq::Record) -> crate::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }
        self.row_count += 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_fastq_with_default_fields() {
        let batch_builder = BatchBuilder::new_fastq(Select::All, 10).unwrap();
        assert_eq!(batch_builder.schema().fields().len(), 4);
    }

    #[test]
    fn test_new_fasta_with_default_fields() {
        let batch_builder = BatchBuilder::new_fasta(Select::All, 10).unwrap();
        assert_eq!(batch_builder.schema().fields().len(), 3);
    }

    #[test]
    fn test_schema() {
        let batch_builder = BatchBuilder::new_fastq(Select::All, 10).unwrap();
        let schema = batch_builder.schema();
        assert_eq!(schema.fields().len(), 4);
        assert_eq!(schema.field(0).name(), "name");
        assert_eq!(schema.field(3).name(), "quality");
    }

    #[test]
    fn test_push_fasta_record() {
        let capacity = 10;
        let mut batch_builder = BatchBuilder::new_fasta(Select::All, capacity).unwrap();

        let record = noodles::fasta::Record::new(
            noodles::fasta::record::Definition::new(b"s0", Some(b"description".into())),
            noodles::fasta::record::Sequence::from(b"ACGT".to_vec()),
        );
        batch_builder.push(&record).unwrap();

        let record_batch = batch_builder.finish().unwrap();
        assert_eq!(record_batch.num_rows(), 1);
    }

    #[test]
    fn test_push_fastq_record() {
        let capacity = 10;
        let mut batch_builder = BatchBuilder::new_fastq(Select::All, capacity).unwrap();

        let record = noodles::fastq::Record::new(
            noodles::fastq::record::Definition::new(b"s0", b""),
            b"ACGT",
            b"!!!!",
        );
        batch_builder.push(&record).unwrap();

        let record_batch = batch_builder.finish().unwrap();
        assert_eq!(record_batch.num_rows(), 1);
    }

    #[test]
    fn test_finish_empty_batch() {
        let capacity = 10;
        let mut batch_builder = BatchBuilder::new_fastq(Select::All, capacity).unwrap();

        let record_batch = batch_builder.finish().unwrap();
        assert_eq!(record_batch.num_rows(), 0);
        assert_eq!(record_batch.num_columns(), 4);
    }
}

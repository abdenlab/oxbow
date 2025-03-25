use std::io;

use arrow::array::ArrayRef;
use arrow::datatypes::{Field as ArrowField, Schema};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use indexmap::IndexMap;

use super::field::Push as _;
use super::field::{Field, FieldBuilder, FASTA_DEFAULT_FIELD_NAMES, FASTQ_DEFAULT_FIELD_NAMES};

/// A builder for Arrow record batches of sequence records.
pub struct BatchBuilder {
    fields: Vec<Field>,
    field_builders: IndexMap<Field, FieldBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` for FASTQ records.
    ///
    /// # Arguments
    /// * `field_names` - Optional vector of field names to project. If `None`, the default field
    ///     names are used.
    /// * `capacity` - The number of rows to preallocate for a batch.
    ///
    /// # Returns
    /// A `Result` containing the `BatchBuilder` or an `io::Error` if the field names are invalid.
    pub fn new_fastq(field_names: Option<Vec<String>>, capacity: usize) -> io::Result<Self> {
        let default_field_names = FASTQ_DEFAULT_FIELD_NAMES
            .iter()
            .map(|s| s.to_string())
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

        Ok(Self {
            fields,
            field_builders,
        })
    }

    /// Creates a new `BatchBuilder` for FASTA records.
    ///
    /// # Arguments
    /// * `field_names` - Optional vector of field names to project. If `None`, the default field
    ///     names are used.
    /// * `capacity` - The number of rows to preallocate for a batch.
    ///
    /// # Returns
    /// A `Result` containing the `BatchBuilder` or an `io::Error` if the field names are invalid.
    pub fn new_fasta(field_names: Option<Vec<String>>, capacity: usize) -> io::Result<Self> {
        let default_field_names = FASTA_DEFAULT_FIELD_NAMES
            .iter()
            .map(|s| s.to_string())
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
            .field_builders
            .iter_mut()
            .map(|(field, builder)| {
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

/// Append a FASTA record to the batch.
impl Push<&noodles::fasta::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::fasta::Record) -> io::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }
        Ok(())
    }
}

/// Append a FASTQ record to the batch.
impl Push<&noodles::fastq::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::fastq::Record) -> io::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_fastq_with_default_fields() {
        let capacity = 10;
        let batch_builder = BatchBuilder::new_fastq(None, capacity).unwrap();

        assert_eq!(batch_builder.fields.len(), FASTQ_DEFAULT_FIELD_NAMES.len());
        for (field, default_name) in batch_builder.fields.iter().zip(FASTQ_DEFAULT_FIELD_NAMES) {
            assert_eq!(field.name(), default_name);
        }
    }

    #[test]
    fn test_new_fasta_with_default_fields() {
        let capacity = 10;
        let batch_builder = BatchBuilder::new_fasta(None, capacity).unwrap();

        assert_eq!(batch_builder.fields.len(), FASTA_DEFAULT_FIELD_NAMES.len());
        for (field, default_name) in batch_builder.fields.iter().zip(FASTA_DEFAULT_FIELD_NAMES) {
            assert_eq!(field.name(), default_name);
        }
    }

    #[test]
    fn test_get_arrow_schema() {
        let capacity = 10;
        let batch_builder = BatchBuilder::new_fastq(None, capacity).unwrap();

        let schema = batch_builder.get_arrow_schema();
        assert_eq!(schema.fields().len(), FASTQ_DEFAULT_FIELD_NAMES.len());
        for (field, default_name) in schema.fields().iter().zip(FASTQ_DEFAULT_FIELD_NAMES) {
            assert_eq!(field.name(), default_name);
        }
    }

    #[test]
    fn test_push_fasta_record() {
        let capacity = 10;
        let mut batch_builder = BatchBuilder::new_fasta(None, capacity).unwrap();

        let record = noodles::fasta::Record::new(
            noodles::fasta::record::Definition::new(b"s0", Some(b"description".to_vec())),
            noodles::fasta::record::Sequence::from(b"ACGT".to_vec()),
        );
        batch_builder.push(&record).unwrap();

        let record_batch = batch_builder.finish().unwrap();
        assert_eq!(record_batch.num_rows(), 1);
    }

    #[test]
    fn test_push_fastq_record() {
        let capacity = 10;
        let mut batch_builder = BatchBuilder::new_fastq(None, capacity).unwrap();

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
        let mut batch_builder = BatchBuilder::new_fastq(None, capacity).unwrap();

        let record_batch = batch_builder.finish().unwrap();
        assert_eq!(record_batch.num_rows(), 0);
        assert_eq!(record_batch.num_columns(), FASTQ_DEFAULT_FIELD_NAMES.len());
    }
}

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
    /// names are used.
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
            .unwrap_or_else(|| default_field_names)
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
    /// names are used.
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
            .unwrap_or_else(|| default_field_names)
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
        RecordBatch::try_from_iter(name_to_array.into_iter())
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

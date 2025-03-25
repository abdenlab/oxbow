use std::io;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::GenericStringBuilder;
use arrow::datatypes::DataType;

pub const FASTQ_DEFAULT_FIELD_NAMES: [&str; 4] = ["name", "description", "sequence", "quality"];
pub const FASTA_DEFAULT_FIELD_NAMES: [&str; 3] = ["name", "description", "sequence"];

/// A sequence (FASTA/FASTQ) standard field.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum Field {
    Name,
    Descr,
    Seq,
    Qual,
}

impl Field {
    pub fn name(&self) -> &str {
        match self {
            Self::Name => "name",
            Self::Descr => "description",
            Self::Seq => "sequence",
            Self::Qual => "quality",
        }
    }

    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::Name => DataType::Utf8,
            Self::Descr => DataType::Utf8,
            Self::Seq => DataType::Utf8,
            Self::Qual => DataType::Utf8,
        }
    }

    pub fn get_arrow_field(&self) -> arrow::datatypes::Field {
        arrow::datatypes::Field::new(self.name(), self.arrow_type(), true)
    }
}

impl FromStr for Field {
    type Err = std::io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "name" => Ok(Self::Name),
            "description" => Ok(Self::Descr),
            "sequence" => Ok(Self::Seq),
            "quality" => Ok(Self::Qual),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Invalid field name: {}", s),
            )),
        }
    }
}

/// A builder for an Arrow array (column) corresponding to a sequence standard field.
pub enum FieldBuilder {
    Name(GenericStringBuilder<i32>),
    Descr(GenericStringBuilder<i32>),
    Seq(GenericStringBuilder<i32>),
    Qual(GenericStringBuilder<i32>),
}

impl FieldBuilder {
    /// Creates a new `FieldBuilder` for the specified field with the given capacity.
    ///
    /// # Arguments
    /// * `field` - The field for which to create the builder.
    /// * `capacity` - The number of rows to preallocate for a batch.
    ///
    /// # Returns
    /// A `FieldBuilder` for the specified field.
    pub fn new(field: Field, capacity: usize) -> Self {
        // Set a default initial data capacity for the string builders.
        // This is the number of bytes preallocated for all items in the column, not per item.
        let data_capacity = 1024 * 1024;
        match field {
            Field::Name => Self::Name(GenericStringBuilder::<i32>::with_capacity(
                capacity,
                data_capacity,
            )),
            Field::Descr => Self::Descr(GenericStringBuilder::<i32>::with_capacity(
                capacity,
                data_capacity,
            )),
            Field::Seq => Self::Seq(GenericStringBuilder::<i32>::with_capacity(
                capacity,
                data_capacity,
            )),
            Field::Qual => Self::Qual(GenericStringBuilder::<i32>::with_capacity(
                capacity,
                data_capacity,
            )),
        }
    }

    pub fn finish(&mut self) -> arrow::array::ArrayRef {
        match self {
            Self::Name(builder) => Arc::new(builder.finish()),
            Self::Descr(builder) => Arc::new(builder.finish()),
            Self::Seq(builder) => Arc::new(builder.finish()),
            Self::Qual(builder) => Arc::new(builder.finish()),
        }
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Append a field value from a FASTQ record to the column.
impl Push<&noodles::fastq::Record> for FieldBuilder {
    fn push(&mut self, record: &noodles::fastq::Record) -> io::Result<()> {
        match self {
            Self::Name(builder) => builder.append_value(record.name().to_string()),
            Self::Descr(builder) => builder.append_value(record.description().to_string()),
            Self::Seq(builder) => {
                let seq = record.sequence();
                let seq_str = std::str::from_utf8(seq).ok();
                builder.append_option(seq_str);
            }
            Self::Qual(builder) => {
                let qual = record.quality_scores();
                let qual_str = std::str::from_utf8(qual).ok();
                builder.append_option(qual_str);
            }
        }
        Ok(())
    }
}

/// Append a field value from a FASTA record to the column.
impl Push<&noodles::fasta::Record> for FieldBuilder {
    fn push(&mut self, record: &noodles::fasta::Record) -> io::Result<()> {
        match self {
            Self::Name(builder) => {
                let name = record.name();
                let name_str = std::str::from_utf8(name).ok();
                builder.append_option(name_str);
            }
            Self::Descr(builder) => {
                let descr_str = record
                    .description()
                    .and_then(|s| std::str::from_utf8(s).ok());
                builder.append_option(descr_str);
            }
            Self::Seq(builder) => {
                let seq = record.sequence().as_ref();
                let seq_str = std::str::from_utf8(seq).ok();
                builder.append_option(seq_str);
            }
            Self::Qual(_) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "FASTA records do not have quality scores",
                ));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::StringArray;

    #[test]
    fn test_field_arrow_type() {
        for field in [Field::Name, Field::Descr, Field::Seq, Field::Qual] {
            let mut builder = FieldBuilder::new(field.clone(), 10);
            let data_type = builder.finish().data_type().clone();
            assert_eq!(field.arrow_type(), data_type);
        }
    }

    #[test]
    fn test_field_from_str() {
        assert_eq!(Field::from_str("name").unwrap(), Field::Name);
        assert_eq!(Field::from_str("description").unwrap(), Field::Descr);
        assert_eq!(Field::from_str("sequence").unwrap(), Field::Seq);
        assert_eq!(Field::from_str("quality").unwrap(), Field::Qual);
        assert!(Field::from_str("invalid").is_err());
    }

    #[test]
    fn test_field_builder_push_fastq() {
        let record = noodles::fastq::Record::new(
            noodles::fastq::record::Definition::new(b"s0", b""),
            b"ACGT",
            b"!!!!",
        );
        let mut builder = FieldBuilder::new(Field::Name, 10);
        builder.push(&record).unwrap();
        let array = builder.finish();
        let array = array.as_any().downcast_ref::<StringArray>().unwrap();

        assert_eq!(array.value(0), "s0");
    }

    #[test]
    fn test_field_builder_push_fasta() {
        let record = noodles::fasta::Record::new(
            noodles::fasta::record::Definition::new(b"s0", Some(b"description".to_vec())),
            noodles::fasta::record::Sequence::from(b"ACGT".to_vec()),
        );
        let mut builder = FieldBuilder::new(Field::Name, 10);
        builder.push(&record).unwrap();
        let array = builder.finish();
        let array = array.as_any().downcast_ref::<StringArray>().unwrap();

        assert_eq!(array.value(0), "s0");
    }

    #[test]
    fn test_field_builder_push_fasta_no_quality() {
        let record = noodles::fasta::Record::new(
            noodles::fasta::record::Definition::new(b"s0", Some(b"description".to_vec())),
            noodles::fasta::record::Sequence::from(b"ACGT".to_vec()),
        );
        let mut builder = FieldBuilder::new(Field::Qual, 10);

        assert!(builder.push(&record).is_err());
    }
}

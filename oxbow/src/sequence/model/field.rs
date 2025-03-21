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
        match s {
            "name" => Ok(Self::Name),
            "description" => Ok(Self::Descr),
            "sequence" => Ok(Self::Seq),
            "quality" => Ok(Self::Qual),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Invalid field name: {}", s.to_string()),
            )),
        }
    }
}

/// Builds an Arrow array (column) corresponding to a sequence standard field.
pub enum FieldBuilder {
    Name(GenericStringBuilder<i32>),
    Descr(GenericStringBuilder<i32>),
    Seq(GenericStringBuilder<i32>),
    Qual(GenericStringBuilder<i32>),
}

impl FieldBuilder {
    pub fn new(field: Field, capacity: usize, data_capacity: usize) -> Self {
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
                let seq_str = std::str::from_utf8(seq.as_ref()).ok();
                builder.append_option(seq_str);
            }
            Self::Qual(builder) => {
                let qual = record.quality_scores();
                let qual_str = std::str::from_utf8(qual.as_ref()).ok();
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
                let seq = record.sequence();
                let seq_str = std::str::from_utf8(seq.as_ref()).ok();
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

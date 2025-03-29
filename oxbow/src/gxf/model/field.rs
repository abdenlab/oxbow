use std::io;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::{ArrayRef, Float32Builder, GenericStringBuilder, Int32Builder, UInt8Builder};
use arrow::datatypes::DataType;
use arrow::datatypes::Field as ArrowField;

pub const DEFAULT_FIELD_NAMES: [&str; 8] = [
    "seqid", "source", "type", "start", "end", "score", "strand", "frame",
];

/// An GXF feature standard field.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Field {
    SeqId,
    Source,
    Type,
    Start,
    End,
    Score,
    Strand,
    Frame,
}

impl Field {
    pub fn name(&self) -> &str {
        match self {
            Self::SeqId => "seqid",
            Self::Source => "source",
            Self::Type => "type",
            Self::Start => "start",
            Self::End => "end",
            Self::Score => "score",
            Self::Strand => "strand",
            Self::Frame => "frame",
        }
    }

    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::SeqId => DataType::Utf8,
            Self::Source => DataType::Utf8,
            Self::Type => DataType::Utf8,
            Self::Start => DataType::Int32,
            Self::End => DataType::Int32,
            Self::Score => DataType::Float32,
            Self::Strand => DataType::Utf8,
            Self::Frame => DataType::UInt8,
        }
    }

    pub fn get_arrow_field(&self) -> ArrowField {
        ArrowField::new(self.name(), self.arrow_type(), true)
    }
}

impl FromStr for Field {
    type Err = std::io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "seqid" => Ok(Self::SeqId),
            "source" => Ok(Self::Source),
            "type" => Ok(Self::Type),
            "start" => Ok(Self::Start),
            "end" => Ok(Self::End),
            "score" => Ok(Self::Score),
            "strand" => Ok(Self::Strand),
            "frame" => Ok(Self::Frame),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Invalid field name: {}", s),
            )),
        }
    }
}

/// A builder for an Arrow array (column) corresponding to a GXF feature standard field.
pub enum FieldBuilder {
    SeqId(GenericStringBuilder<i32>),
    Source(GenericStringBuilder<i32>),
    Type(GenericStringBuilder<i32>),
    Start(Int32Builder),
    End(Int32Builder),
    Score(Float32Builder),
    Strand(GenericStringBuilder<i32>),
    Frame(UInt8Builder),
}

impl FieldBuilder {
    /// Creates a new `FieldBuilder` for the specified field with the given capacity.
    ///
    /// # Arguments
    /// * `field` - The field to build.
    /// * `capacity` - The number of rows to preallocate for a batch.
    pub fn new(field: Field, capacity: usize) -> Self {
        match field {
            Field::SeqId => Self::SeqId(GenericStringBuilder::<i32>::with_capacity(capacity, 1024)),
            Field::Source => {
                Self::Source(GenericStringBuilder::<i32>::with_capacity(capacity, 1024))
            }
            Field::Type => Self::Type(GenericStringBuilder::<i32>::with_capacity(capacity, 1024)),
            Field::Start => Self::Start(Int32Builder::with_capacity(capacity)),
            Field::End => Self::End(Int32Builder::with_capacity(capacity)),
            Field::Score => Self::Score(Float32Builder::with_capacity(capacity)),
            Field::Strand => {
                Self::Strand(GenericStringBuilder::<i32>::with_capacity(capacity, 1024))
            }
            Field::Frame => Self::Frame(UInt8Builder::with_capacity(capacity)),
        }
    }

    pub fn finish(&mut self) -> ArrayRef {
        match self {
            Self::SeqId(builder) => Arc::new(builder.finish()),
            Self::Source(builder) => Arc::new(builder.finish()),
            Self::Type(builder) => Arc::new(builder.finish()),
            Self::Start(builder) => Arc::new(builder.finish()),
            Self::End(builder) => Arc::new(builder.finish()),
            Self::Score(builder) => Arc::new(builder.finish()),
            Self::Strand(builder) => Arc::new(builder.finish()),
            Self::Frame(builder) => Arc::new(builder.finish()),
        }
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Append a field value from a GFF record to the column.
impl<'a> Push<&'a noodles::gff::Record<'a>> for FieldBuilder {
    fn push(&mut self, record: &noodles::gff::Record) -> io::Result<()> {
        match self {
            Self::SeqId(builder) => {
                let seq_id = record.reference_sequence_name();
                builder.append_value(seq_id);
            }
            Self::Source(builder) => {
                let source = record.source();
                builder.append_value(source);
            }
            Self::Type(builder) => {
                let ty = record.ty();
                builder.append_value(ty);
            }
            Self::Start(builder) => {
                let start = record.start().ok().map(|pos| usize::from(pos) as i32);
                builder.append_option(start);
            }
            Self::End(builder) => {
                let end = record.start().ok().map(|pos| usize::from(pos) as i32);
                builder.append_option(end);
            }
            Self::Score(builder) => {
                let score = record.score().and_then(|score| score.ok());
                builder.append_option(score);
            }
            Self::Strand(builder) => {
                let strand = match record.strand() {
                    Ok(noodles::gff::record::Strand::Forward) => Some("+"),
                    Ok(noodles::gff::record::Strand::Reverse) => Some("-"),
                    _ => None,
                };
                builder.append_option(strand);
            }
            Self::Frame(builder) => {
                let frame = match record.phase() {
                    Some(Ok(noodles::gff::record::Phase::Zero)) => Some(0),
                    Some(Ok(noodles::gff::record::Phase::One)) => Some(1),
                    Some(Ok(noodles::gff::record::Phase::Two)) => Some(2),
                    _ => None,
                };
                builder.append_option(frame);
            }
        }
        Ok(())
    }
}

/// Append a field value from a GTF record to the column.
impl Push<&noodles::gtf::Record> for FieldBuilder {
    fn push(&mut self, record: &noodles::gtf::Record) -> io::Result<()> {
        match self {
            Self::SeqId(builder) => {
                let seq_id = record.reference_sequence_name();
                builder.append_value(seq_id);
            }
            Self::Source(builder) => {
                let source = record.source();
                builder.append_value(source);
            }
            Self::Type(builder) => {
                let ty = record.ty();
                builder.append_value(ty);
            }
            Self::Start(builder) => {
                let start = usize::from(record.start()) as i32;
                builder.append_value(start);
            }
            Self::End(builder) => {
                let end = usize::from(record.end()) as i32;
                builder.append_value(end);
            }
            Self::Score(builder) => {
                builder.append_option(record.score());
            }
            Self::Strand(builder) => {
                let strand = match record.strand() {
                    Some(noodles::gtf::record::Strand::Forward) => Some("+"),
                    Some(noodles::gtf::record::Strand::Reverse) => Some("-"),
                    _ => None,
                };
                builder.append_option(strand);
            }
            Self::Frame(builder) => {
                let frame = record.frame().map(|frame| {
                    let n: u8 = frame.into();
                    n
                });
                builder.append_option(frame);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_field_arrow_type() {
        for field in [
            Field::SeqId,
            Field::Source,
            Field::Type,
            Field::Start,
            Field::End,
            Field::Score,
            Field::Strand,
            Field::Frame,
        ] {
            let mut builder = FieldBuilder::new(field.clone(), 10);
            let data_type = builder.finish().data_type().clone();
            assert_eq!(field.arrow_type(), data_type);
        }
    }

    #[test]
    fn test_field_from_str() {
        assert_eq!(Field::from_str("seqid").unwrap(), Field::SeqId);
        assert_eq!(Field::from_str("type").unwrap(), Field::Type);
        assert_eq!(Field::from_str("source").unwrap(), Field::Source);
        assert_eq!(Field::from_str("start").unwrap(), Field::Start);
        assert_eq!(Field::from_str("end").unwrap(), Field::End);
        assert_eq!(Field::from_str("score").unwrap(), Field::Score);
        assert_eq!(Field::from_str("strand").unwrap(), Field::Strand);
        assert_eq!(Field::from_str("frame").unwrap(), Field::Frame);
        assert!(Field::from_str("invalid").is_err());
    }

    #[test]
    fn test_field_builder_push_gtf() {
        for field in [
            Field::SeqId,
            Field::Source,
            Field::Type,
            Field::Start,
            Field::End,
            Field::Score,
            Field::Strand,
            Field::Frame,
        ] {
            let mut builder = FieldBuilder::new(field, 10);
            let record = noodles::gtf::Record::default();
            assert!(builder.push(&record).is_ok());
        }
    }
}

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
            "seqid" => Ok(Field::SeqId),
            "source" => Ok(Field::Source),
            "type" => Ok(Field::Type),
            "start" => Ok(Field::Start),
            "end" => Ok(Field::End),
            "score" => Ok(Field::Score),
            "strand" => Ok(Field::Strand),
            "frame" => Ok(Field::Frame),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Invalid field name: {}", s.to_string()),
            )),
        }
    }
}

/// Builds an Arrow array (column) corresponding to a GXF feature standard field.
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
    pub fn new(field: Field, capacity: usize) -> Self {
        match field {
            Field::SeqId => Self::SeqId(GenericStringBuilder::<i32>::new()),
            Field::Source => Self::Source(GenericStringBuilder::<i32>::new()),
            Field::Type => Self::Type(GenericStringBuilder::<i32>::new()),
            Field::Start => Self::Start(Int32Builder::with_capacity(capacity)),
            Field::End => Self::End(Int32Builder::with_capacity(capacity)),
            Field::Score => Self::Score(Float32Builder::with_capacity(capacity)),
            Field::Strand => Self::Strand(GenericStringBuilder::<i32>::new()),
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
                let frame = record.frame().and_then(|frame| {
                    let n: u8 = frame.into();
                    Some(n)
                });
                builder.append_option(frame);
            }
        }
        Ok(())
    }
}

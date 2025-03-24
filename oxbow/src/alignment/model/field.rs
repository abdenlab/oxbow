use std::io::{self, Write};
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, GenericStringBuilder, Int32Builder, StringArray, StringDictionaryBuilder,
    UInt16Builder, UInt8Builder,
};
use arrow::datatypes::{DataType, Int32Type};
use arrow::error::ArrowError;

use noodles::sam::alignment::record::{Cigar, QualityScores, Sequence};
use noodles::sam::alignment::Record;

pub const DEFAULT_FIELD_NAMES: [&str; 12] = [
    "qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual",
    "end",
];

/// An alignment (SAM/BAM/CRAM) standard field.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Field {
    Qname,
    Flag,
    Rname,
    Pos,
    Mapq,
    Cigar,
    Rnext,
    Pnext,
    Tlen,
    Seq,
    Qual,
    End,
}

impl Field {
    pub fn name(&self) -> &str {
        match self {
            Self::Qname => "qname",
            Self::Flag => "flag",
            Self::Rname => "rname",
            Self::Pos => "pos",
            Self::Mapq => "mapq",
            Self::Cigar => "cigar",
            Self::Rnext => "rnext",
            Self::Pnext => "pnext",
            Self::Tlen => "tlen",
            Self::Seq => "seq",
            Self::Qual => "qual",
            Self::End => "end",
        }
    }

    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::Qname => DataType::Utf8,
            Self::Flag => DataType::UInt16,
            Self::Rname => {
                DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8))
            }
            Self::Pos => DataType::Int32,
            Self::Mapq => DataType::UInt8,
            Self::Cigar => DataType::Utf8,
            Self::Rnext => {
                DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8))
            }
            Self::Pnext => DataType::Int32,
            Self::Tlen => DataType::Int32,
            Self::Seq => DataType::Utf8,
            Self::Qual => DataType::Utf8,
            Self::End => DataType::Int32,
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
            "qname" => Ok(Self::Qname),
            "flag" => Ok(Self::Flag),
            "rname" => Ok(Self::Rname),
            "pos" => Ok(Self::Pos),
            "mapq" => Ok(Self::Mapq),
            "cigar" => Ok(Self::Cigar),
            "rnext" => Ok(Self::Rnext),
            "pnext" => Ok(Self::Pnext),
            "tlen" => Ok(Self::Tlen),
            "seq" => Ok(Self::Seq),
            "qual" => Ok(Self::Qual),
            "end" => Ok(Self::End),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Invalid field name: {}", s),
            )),
        }
    }
}

/// Builds an Arrow array (column) corresponding to an alignment standard field.
pub enum FieldBuilder {
    Qname(GenericStringBuilder<i32>),
    Flag(UInt16Builder),
    Rname(StringDictionaryBuilder<Int32Type>),
    Pos(Int32Builder),
    Mapq(UInt8Builder),
    Cigar(GenericStringBuilder<i32>),
    Rnext(StringDictionaryBuilder<Int32Type>),
    Pnext(Int32Builder),
    Tlen(Int32Builder),
    Seq(GenericStringBuilder<i32>),
    Qual(GenericStringBuilder<i32>),
    End(Int32Builder),
}

impl FieldBuilder {
    pub fn new(field: Field, capacity: usize, ref_names: &[String]) -> Result<Self, ArrowError> {
        let field = match field {
            Field::Qname => Self::Qname(GenericStringBuilder::<i32>::with_capacity(capacity, 1024)),
            Field::Flag => Self::Flag(UInt16Builder::with_capacity(capacity)),
            Field::Rname => {
                let refs = StringArray::from(ref_names.to_owned());
                Self::Rname(StringDictionaryBuilder::<Int32Type>::new_with_dictionary(
                    capacity, &refs,
                )?)
            }
            Field::Pos => Self::Pos(Int32Builder::with_capacity(capacity)),
            Field::Mapq => Self::Mapq(UInt8Builder::with_capacity(capacity)),
            Field::Cigar => Self::Cigar(GenericStringBuilder::<i32>::with_capacity(capacity, 1024)),
            Field::Rnext => {
                let refs = StringArray::from(ref_names.to_owned());
                Self::Rnext(StringDictionaryBuilder::<Int32Type>::new_with_dictionary(
                    capacity, &refs,
                )?)
            }
            Field::Pnext => Self::Pnext(Int32Builder::with_capacity(capacity)),
            Field::Tlen => Self::Tlen(Int32Builder::with_capacity(capacity)),
            Field::Seq => Self::Seq(GenericStringBuilder::<i32>::with_capacity(capacity, 1024)),
            Field::Qual => Self::Qual(GenericStringBuilder::<i32>::with_capacity(capacity, 1024)),
            Field::End => Self::End(Int32Builder::with_capacity(capacity)),
        };
        Ok(field)
    }

    pub fn finish(&mut self) -> ArrayRef {
        match self {
            Self::Qname(builder) => Arc::new(builder.finish()),
            Self::Flag(builder) => Arc::new(builder.finish()),
            Self::Rname(builder) => Arc::new(builder.finish()),
            Self::Pos(builder) => Arc::new(builder.finish()),
            Self::Mapq(builder) => Arc::new(builder.finish()),
            Self::Cigar(builder) => Arc::new(builder.finish()),
            Self::Rnext(builder) => Arc::new(builder.finish()),
            Self::Pnext(builder) => Arc::new(builder.finish()),
            Self::Tlen(builder) => Arc::new(builder.finish()),
            Self::Seq(builder) => Arc::new(builder.finish()),
            Self::Qual(builder) => Arc::new(builder.finish()),
            Self::End(builder) => Arc::new(builder.finish()),
        }
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T, header: &noodles::sam::Header) -> io::Result<()>;
}

/// Append a field value from a SAM record to the column.
impl Push<&noodles::sam::Record> for FieldBuilder {
    fn push(
        &mut self,
        record: &noodles::sam::Record,
        header: &noodles::sam::Header,
    ) -> io::Result<()> {
        match self {
            Self::Qname(builder) => {
                builder.append_option(record.name().map(|name| name.to_string()));
            }
            Self::Flag(builder) => {
                builder.append_option(record.flags().ok().map(|flags| flags.bits()));
            }
            Self::Rname(builder) => {
                let rname = record
                    .reference_sequence(header)
                    .and_then(|result| result.ok().map(|(name, _)| name.to_string()));
                builder.append_option(rname);
            }
            Self::Pos(builder) => {
                let start = record
                    .alignment_start()
                    .and_then(|result| result.ok().map(|pos| pos.get() as i32));
                builder.append_option(start);
            }
            Self::Mapq(builder) => {
                let mapq = record
                    .mapping_quality()
                    .and_then(|result| result.ok().map(|mq| mq.into()));
                builder.append_option(mapq);
            }
            Self::Cigar(builder) => {
                builder.append_option(get_cigar(record.cigar()));
            }
            Self::Rnext(builder) => {
                let rnext = record
                    .mate_reference_sequence(header)
                    .and_then(|result| result.ok().map(|(name, _)| name.to_string()));
                builder.append_option(rnext);
            }
            Self::Pnext(builder) => {
                let start = record
                    .mate_alignment_start()
                    .and_then(|result| result.ok().map(|pos| pos.get() as i32));
                builder.append_option(start);
            }
            Self::Tlen(builder) => {
                builder.append_option(record.template_length().ok());
            }
            Self::Seq(builder) => {
                let seq = record.sequence().iter().collect();
                builder.append_option(String::from_utf8(seq).ok());
            }
            Self::Qual(builder) => {
                builder.append_option(get_quality_scores(record.quality_scores()));
            }
            Self::End(builder) => {
                let end = record
                    .alignment_end()
                    .and_then(|result| result.ok().map(|pos| pos.get() as i32));
                builder.append_option(end);
            }
        }
        Ok(())
    }
}

/// Append a field value from a BAM record to the column.
impl Push<&noodles::bam::Record> for FieldBuilder {
    fn push(
        &mut self,
        record: &noodles::bam::Record,
        header: &noodles::sam::Header,
    ) -> io::Result<()> {
        match self {
            Self::Qname(builder) => {
                builder.append_option(record.name().map(|name| name.to_string()));
            }
            Self::Flag(builder) => {
                builder.append_value(record.flags().bits());
            }
            Self::Rname(builder) => {
                let rname = record
                    .reference_sequence(header)
                    .and_then(|result| result.ok().map(|(name, _)| name.to_string()));
                builder.append_option(rname);
            }
            Self::Pos(builder) => {
                let start = record
                    .alignment_start()
                    .and_then(|result| result.ok().map(|pos| pos.get() as i32));
                builder.append_option(start);
            }
            Self::Mapq(builder) => {
                builder.append_option(record.mapping_quality().map(|mq| mq.into()));
            }
            Self::Cigar(builder) => {
                builder.append_option(get_cigar(record.cigar()));
            }
            Self::Rnext(builder) => {
                let rnext = record
                    .mate_reference_sequence(header)
                    .and_then(|result| result.ok().map(|(name, _)| name.to_string()));
                builder.append_option(rnext);
            }
            Self::Pnext(builder) => {
                let start = record
                    .mate_alignment_start()
                    .and_then(|result| result.ok().map(|pos| pos.get() as i32));
                builder.append_option(start);
            }
            Self::Tlen(builder) => {
                builder.append_value(record.template_length());
            }
            Self::Seq(builder) => {
                let seq = record.sequence().iter().collect::<Vec<u8>>();
                builder.append_option(String::from_utf8(seq).ok());
            }
            Self::Qual(builder) => {
                builder.append_option(get_quality_scores(record.quality_scores()));
            }
            Self::End(builder) => {
                let end = record
                    .alignment_end()
                    .and_then(|result| result.ok().map(|pos| pos.get() as i32));
                builder.append_option(end);
            }
        }
        Ok(())
    }
}

/// Serialize parsed quality scores to a string.
fn get_quality_scores(quality_scores: impl QualityScores) -> Option<String> {
    const OFFSET: u8 = b'!';
    const MAX_SCORE: u8 = b'~' - OFFSET;

    if quality_scores.is_empty() {
        return None;
    }

    let mut buf = Vec::new();
    for result in quality_scores.iter() {
        let score = match result {
            Ok(score) => score,
            Err(_) => return None,
        };

        if score <= MAX_SCORE {
            // SAFETY: `score` <= 93.
            let m = score + OFFSET;
            buf.write_all(&[m]).ok();
        } else {
            return None;
        }
    }

    String::from_utf8(buf).ok()
}

/// Serialize parsed CIGAR to a string.
pub fn get_cigar(cigar: impl Cigar) -> Option<String> {
    use noodles::sam::alignment::record::cigar::op::Kind;

    if cigar.is_empty() {
        return None;
    }

    let mut buf = Vec::new();
    for result in cigar.iter() {
        let op = match result {
            Ok(op) => op,
            Err(_) => return None,
        };

        buf.write_all(op.len().to_string().as_bytes()).ok();
        let c = match op.kind() {
            Kind::Match => b'M',
            Kind::Insertion => b'I',
            Kind::Deletion => b'D',
            Kind::Skip => b'N',
            Kind::SoftClip => b'S',
            Kind::HardClip => b'H',
            Kind::Pad => b'P',
            Kind::SequenceMatch => b'=',
            Kind::SequenceMismatch => b'X',
        };
        buf.write_all(&[c]).ok();
    }

    String::from_utf8(buf).ok()
}

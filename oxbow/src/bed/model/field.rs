use std::io;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::{FixedSizeListBuilder, ListBuilder};
use arrow::array::{
    Float64Builder, GenericStringBuilder, Int64Builder, StringArray, StringDictionaryBuilder,
    UInt16Builder, UInt8Builder,
};
use arrow::datatypes::{DataType, Field as ArrowField, UInt8Type};

use crate::util::reset_dictarray_builder;

pub const DEFAULT_FIELD_NAMES: [&str; 12] = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts",
];

/// A BED field.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Field {
    // Mandatory fields - coordinates
    Chrom,
    Start,
    End,
    // Optional fields - simple attributes
    Name,
    Score,
    Strand,
    // Optional fields - display attributes
    ThickStart,
    ThickEnd,
    ItemRgb,
    // Optional fields - block attributes
    BlockCount,
    BlockSizes,
    BlockStarts,
    // Non-standard fields
    FloatScore, // Non-standard score field
    Value,      // Fourth column of bedGraph
}

impl Field {
    pub fn name(&self) -> &str {
        match self {
            Self::Chrom => "chrom",
            Self::Start => "start",
            Self::End => "end",
            Self::Name => "name",
            Self::Score => "score",
            Self::Strand => "strand",
            Self::ThickStart => "thickStart",
            Self::ThickEnd => "thickEnd",
            Self::ItemRgb => "itemRgb",
            Self::BlockCount => "blockCount",
            Self::BlockSizes => "blockSizes",
            Self::BlockStarts => "blockStarts",
            Self::FloatScore => "score",
            Self::Value => "value",
        }
    }

    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::Chrom => DataType::Utf8,
            Self::Start => DataType::Int64,
            Self::End => DataType::Int64,
            Self::Name => DataType::Utf8,
            Self::Score => DataType::UInt16,
            Self::Strand => {
                DataType::Dictionary(Box::new(DataType::UInt8), Box::new(DataType::Utf8))
            }
            Self::ThickStart => DataType::Int64,
            Self::ThickEnd => DataType::Int64,
            Self::ItemRgb => {
                let item = ArrowField::new("item", DataType::UInt8, true);
                DataType::FixedSizeList(Arc::new(item), 3)
            }
            Self::BlockCount => DataType::Int64,
            Self::BlockSizes => {
                let item = ArrowField::new("item", DataType::Int64, true);
                DataType::List(Arc::new(item))
            }
            Self::BlockStarts => {
                let item = ArrowField::new("item", DataType::Int64, true);
                DataType::List(Arc::new(item))
            }
            Self::FloatScore => DataType::Float64,
            Self::Value => DataType::Float64,
        }
    }

    pub fn get_arrow_field(&self) -> arrow::datatypes::Field {
        arrow::datatypes::Field::new(self.name(), self.arrow_type(), true)
    }
}

impl FromStr for Field {
    type Err = io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "chrom" => Ok(Self::Chrom),
            "start" => Ok(Self::Start),
            "end" => Ok(Self::End),
            "name" => Ok(Self::Name),
            "score" => Ok(Self::Score),
            "strand" => Ok(Self::Strand),
            "thickStart" => Ok(Self::ThickStart),
            "thickEnd" => Ok(Self::ThickEnd),
            "itemRgb" => Ok(Self::ItemRgb),
            "blockCount" => Ok(Self::BlockCount),
            "blockSizes" => Ok(Self::BlockSizes),
            "blockStarts" => Ok(Self::BlockStarts),
            "score.f64" => Ok(Self::FloatScore),
            "value" => Ok(Self::Value),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid field name: {}", s),
            )),
        }
    }
}

/// A builder for an Arrow array (column) corresponding to a BED field.
pub enum FieldBuilder {
    Chrom(GenericStringBuilder<i32>),
    Start(Int64Builder),
    End(Int64Builder),
    Name(GenericStringBuilder<i32>),
    Score(UInt16Builder),
    Strand(StringDictionaryBuilder<UInt8Type>),
    ThickStart(Int64Builder),
    ThickEnd(Int64Builder),
    ItemRgb(FixedSizeListBuilder<UInt8Builder>),
    BlockCount(Int64Builder),
    BlockSizes(ListBuilder<Int64Builder>),
    BlockStarts(ListBuilder<Int64Builder>),
    FloatScore(Float64Builder),
    Value(Float64Builder),
}

impl FieldBuilder {
    /// Creates a new `FieldBuilder` for the specified field with the given capacity.
    ///
    /// # Arguments
    /// * `field` - The field to build.
    /// * `capacity` - The number of rows to preallocate for a batch.
    pub fn new(field: Field, capacity: usize) -> Self {
        match field {
            Field::Chrom => Self::Chrom(GenericStringBuilder::<i32>::with_capacity(capacity, 1024)),
            Field::Start => Self::Start(Int64Builder::with_capacity(capacity)),
            Field::End => Self::End(Int64Builder::with_capacity(capacity)),
            Field::Name => Self::Name(GenericStringBuilder::<i32>::with_capacity(capacity, 1024)),
            Field::Score => Self::Score(UInt16Builder::with_capacity(capacity)),
            Field::Strand => {
                let strand_values = StringArray::from(vec!["+", "-"]);
                Self::Strand(
                    StringDictionaryBuilder::<UInt8Type>::new_with_dictionary(
                        capacity,
                        &strand_values,
                    )
                    .unwrap(),
                )
            }
            Field::ThickStart => Self::ThickStart(Int64Builder::with_capacity(capacity)),
            Field::ThickEnd => Self::ThickEnd(Int64Builder::with_capacity(capacity)),
            Field::ItemRgb => Self::ItemRgb(FixedSizeListBuilder::new(
                UInt8Builder::with_capacity(capacity),
                3,
            )),
            Field::BlockCount => Self::BlockCount(Int64Builder::with_capacity(capacity)),
            Field::BlockSizes => {
                Self::BlockSizes(ListBuilder::new(Int64Builder::with_capacity(capacity)))
            }
            Field::BlockStarts => {
                Self::BlockStarts(ListBuilder::new(Int64Builder::with_capacity(capacity)))
            }
            Field::FloatScore => Self::FloatScore(Float64Builder::with_capacity(capacity)),
            Field::Value => Self::Value(Float64Builder::with_capacity(capacity)),
        }
    }

    pub fn finish(&mut self) -> arrow::array::ArrayRef {
        match self {
            Self::Chrom(builder) => Arc::new(builder.finish()),
            Self::Start(builder) => Arc::new(builder.finish()),
            Self::End(builder) => Arc::new(builder.finish()),
            Self::Name(builder) => Arc::new(builder.finish()),
            Self::Score(builder) => Arc::new(builder.finish()),
            Self::Strand(builder) => {
                let array = reset_dictarray_builder(builder);
                Arc::new(array)
            }
            Self::ThickStart(builder) => Arc::new(builder.finish()),
            Self::ThickEnd(builder) => Arc::new(builder.finish()),
            Self::ItemRgb(builder) => Arc::new(builder.finish()),
            Self::BlockCount(builder) => Arc::new(builder.finish()),
            Self::BlockSizes(builder) => Arc::new(builder.finish()),
            Self::BlockStarts(builder) => Arc::new(builder.finish()),
            Self::FloatScore(builder) => Arc::new(builder.finish()),
            Self::Value(builder) => Arc::new(builder.finish()),
        }
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Append a field value from a BED record to the column.
impl Push<&noodles::bed::Record<3>> for FieldBuilder {
    fn push(&mut self, record: &noodles::bed::Record<3>) -> io::Result<()> {
        match self {
            Self::Chrom(builder) => {
                builder.append_value(record.reference_sequence_name().to_string());
            }
            Self::Start(builder) => {
                let start = record.feature_start().ok().map(|pos| pos.get() as i64);
                builder.append_option(start);
            }
            Self::End(builder) => {
                let end = record
                    .feature_end()
                    .and_then(|result| result.ok())
                    .map(|pos| pos.get() as i64);
                builder.append_option(end);
            }
            Self::Name(builder) => {
                let other_fields = record.other_fields();
                let name = match other_fields.get(0) {
                    Some(bytes) => match bytes {
                        b"." => None,
                        _ => std::str::from_utf8(bytes).ok().map(|s| s.to_string()),
                    },
                    _ => None,
                };
                builder.append_option(name);
            }
            Self::Value(builder) => {
                let other_fields = record.other_fields();
                let value = match other_fields.get(1) {
                    Some(bytes) => std::str::from_utf8(bytes)
                        .ok()
                        .and_then(|s| s.parse::<f64>().ok()),
                    _ => None,
                };
                builder.append_option(value);
            }
            Self::Score(builder) => {
                let other_fields = record.other_fields();
                let score = match other_fields.get(1) {
                    Some(bytes) => std::str::from_utf8(bytes)
                        .ok()
                        .and_then(|s| s.parse::<u16>().ok()),
                    _ => None,
                };
                builder.append_option(score);
            }
            Self::FloatScore(builder) => {
                let other_fields = record.other_fields();
                let score = match other_fields.get(1) {
                    Some(bytes) => std::str::from_utf8(bytes)
                        .ok()
                        .and_then(|s| s.parse::<f64>().ok()),
                    _ => None,
                };
                builder.append_option(score);
            }
            Self::Strand(builder) => {
                let other_fields = record.other_fields();
                let strand = match other_fields.get(2) {
                    Some(bytes) => std::str::from_utf8(bytes).ok().and_then(|s| match s {
                        "+" => Some("+"),
                        "-" => Some("-"),
                        _ => None,
                    }),
                    _ => None,
                };
                builder.append_option(strand);
            }
            Self::ThickStart(builder) => {
                let other_fields = record.other_fields();
                let thick_start = match other_fields.get(3) {
                    Some(bytes) => std::str::from_utf8(bytes)
                        .ok()
                        .and_then(|s| s.parse::<usize>().ok())
                        .map(|pos| pos as i64),
                    _ => None,
                };
                builder.append_option(thick_start);
            }
            Self::ThickEnd(builder) => {
                let other_fields = record.other_fields();
                let thick_end = match other_fields.get(4) {
                    Some(bytes) => std::str::from_utf8(bytes)
                        .ok()
                        .and_then(|s| s.parse::<usize>().ok())
                        .map(|pos| pos as i64),
                    _ => None,
                };
                builder.append_option(thick_end);
            }
            Self::ItemRgb(builder) => {
                let other_fields = record.other_fields();
                let rgb = match other_fields.get(5) {
                    Some(bytes) => match std::str::from_utf8(bytes).ok() {
                        Some("0") => Some([0; 3]),
                        Some(s) => {
                            let mut rgb = [0; 3];
                            for (i, value) in s.split(',').enumerate() {
                                rgb[i] = value.parse::<u8>().map_err(|_| {
                                    io::Error::new(io::ErrorKind::InvalidData, "Invalid RGB value")
                                })?;
                            }
                            Some(rgb)
                        }
                        _ => None,
                    },
                    _ => None,
                };
                match rgb {
                    Some(rgb) => {
                        for value in rgb.iter() {
                            builder.values().append_value(*value);
                        }
                        builder.append(true);
                    }
                    _ => builder.append(false),
                }
            }
            Self::BlockCount(builder) => {
                let other_fields = record.other_fields();
                let block_count = match other_fields.get(6) {
                    Some(bytes) => std::str::from_utf8(bytes)
                        .ok()
                        .and_then(|s| s.parse::<usize>().ok())
                        .map(|pos| pos as i64),
                    _ => None,
                };
                builder.append_option(block_count);
            }
            Self::BlockSizes(builder) => {
                let other_fields = record.other_fields();
                let block_sizes = match other_fields.get(7) {
                    Some(bytes) => std::str::from_utf8(bytes)
                        .ok()
                        .and_then(|s| s.strip_suffix(","))
                        .map(|s| {
                            s.split(",")
                                .map(|s| s.parse::<usize>().map(|v| v as i64))
                                .collect::<Result<Vec<i64>, _>>()
                                .map_err(|_| {
                                    io::Error::new(
                                        io::ErrorKind::InvalidData,
                                        "Invalid block sizes",
                                    )
                                })
                        })
                        .transpose()?,
                    _ => None,
                };
                match block_sizes {
                    Some(block_sizes) => {
                        for value in block_sizes.iter() {
                            builder.values().append_value(*value);
                        }
                        builder.append(true);
                    }
                    _ => builder.append(false),
                };
            }
            Self::BlockStarts(builder) => {
                let other_fields = record.other_fields();
                let block_starts = match other_fields.get(8) {
                    Some(bytes) => std::str::from_utf8(bytes)
                        .ok()
                        .and_then(|s| s.strip_suffix(","))
                        .map(|s| {
                            s.split(",")
                                .map(|s| s.parse::<usize>().map(|v| v as i64))
                                .collect::<Result<Vec<i64>, _>>()
                                .map_err(|_| {
                                    io::Error::new(
                                        io::ErrorKind::InvalidData,
                                        "Invalid block starts",
                                    )
                                })
                        })
                        .transpose()?,
                    _ => None,
                };
                match block_starts {
                    Some(block_starts) => {
                        for value in block_starts.iter() {
                            builder.values().append_value(*value);
                        }
                        builder.append(true);
                    }
                    _ => builder.append(false),
                };
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
            Field::Chrom,
            Field::Start,
            Field::End,
            Field::Name,
            Field::Score,
            Field::Strand,
            Field::ThickStart,
            Field::ThickEnd,
            Field::ItemRgb,
            Field::BlockCount,
            Field::BlockSizes,
            Field::BlockStarts,
            Field::FloatScore,
            Field::Value,
        ] {
            let mut builder = FieldBuilder::new(field.clone(), 10);
            let data_type = builder.finish().data_type().clone();
            assert_eq!(field.arrow_type(), data_type);
        }
    }

    #[test]
    fn test_field_from_str() {
        assert_eq!(Field::from_str("chrom").unwrap(), Field::Chrom);
        assert_eq!(Field::from_str("start").unwrap(), Field::Start);
        assert_eq!(Field::from_str("end").unwrap(), Field::End);
        assert_eq!(Field::from_str("name").unwrap(), Field::Name);
        assert_eq!(Field::from_str("score").unwrap(), Field::Score);
        assert_eq!(Field::from_str("strand").unwrap(), Field::Strand);
        assert_eq!(Field::from_str("thickStart").unwrap(), Field::ThickStart);
        assert_eq!(Field::from_str("thickEnd").unwrap(), Field::ThickEnd);
        assert_eq!(Field::from_str("itemRgb").unwrap(), Field::ItemRgb);
        assert_eq!(Field::from_str("blockCount").unwrap(), Field::BlockCount);
        assert_eq!(Field::from_str("blockSizes").unwrap(), Field::BlockSizes);
        assert_eq!(Field::from_str("blockStarts").unwrap(), Field::BlockStarts);
        assert_eq!(Field::from_str("value").unwrap(), Field::Value);
        assert!(Field::from_str("unknown").is_err());
    }

    #[test]
    fn test_field_builder_push() {
        for field in [
            Field::Chrom,
            Field::Start,
            Field::End,
            Field::Name,
            Field::Score,
            Field::Strand,
            Field::ThickStart,
            Field::ThickEnd,
            Field::ItemRgb,
            Field::BlockCount,
            Field::BlockSizes,
            Field::BlockStarts,
            Field::FloatScore,
            Field::Value,
        ] {
            let mut builder = FieldBuilder::new(field, 10);
            let record = noodles::bed::Record::default();
            assert!(builder.push(&record).is_ok());
        }
    }
}

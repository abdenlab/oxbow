use std::io;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::{
    Float64Builder, StringArray, StringDictionaryBuilder, UInt32Builder, UInt64Builder,
};
use arrow::datatypes::{DataType, Int32Type};
use arrow::error::ArrowError;

use crate::util::reset_dictarray_builder;

pub const DEFAULT_FIELD_NAMES: [&str; 8] = [
    "chrom",
    "start",
    "end",
    "bases_covered",
    "min",
    "max",
    "sum",
    "sum_squares",
];

/// A BBI zoom summary field.
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub enum Field {
    Chrom,
    Start,
    End,
    BasesCovered, // Count of bases with actual data
    Min,          // Minimum value of items
    Max,          // Maximum value of items
    Sum,          // Sum of values for each base
    SumSquares,   // Sum of squares for each base
}

impl Field {
    pub fn name(&self) -> &str {
        match self {
            Self::Chrom => "chrom",
            Self::Start => "start",
            Self::End => "end",
            Self::BasesCovered => "bases_covered",
            Self::Min => "min",
            Self::Max => "max",
            Self::Sum => "sum",
            Self::SumSquares => "sum_squares",
        }
    }

    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::Chrom => {
                DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8))
            }
            Self::Start => DataType::UInt32,
            Self::End => DataType::UInt32,
            Self::BasesCovered => DataType::UInt64,
            Self::Min => DataType::Float64,
            Self::Max => DataType::Float64,
            Self::Sum => DataType::Float64,
            Self::SumSquares => DataType::Float64,
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
            "bases_covered" => Ok(Self::BasesCovered),
            "min" => Ok(Self::Min),
            "max" => Ok(Self::Max),
            "sum" => Ok(Self::Sum),
            "sum_squares" => Ok(Self::SumSquares),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid field name: {}", s),
            )),
        }
    }
}

/// A builder for an Arrow array (column) corresponding to BBI zoom field.
pub enum FieldBuilder {
    Chrom(StringDictionaryBuilder<Int32Type>),
    Start(UInt32Builder),
    End(UInt32Builder),
    BasesCovered(UInt64Builder),
    Min(Float64Builder),
    Max(Float64Builder),
    Sum(Float64Builder),
    SumSquares(Float64Builder),
}

impl FieldBuilder {
    /// Creates a new `FieldBuilder` for the specified field with the given capacity.
    ///
    /// # Arguments
    /// * `field` - The field to build.
    /// * `capacity` - The number of rows to preallocate for a batch.
    pub fn new(field: Field, ref_names: &[String], capacity: usize) -> Result<Self, ArrowError> {
        let builder = match field {
            Field::Chrom => {
                let refs = StringArray::from(ref_names.to_owned());
                Self::Chrom(StringDictionaryBuilder::<Int32Type>::new_with_dictionary(
                    capacity, &refs,
                )?)
            }
            Field::Start => Self::Start(UInt32Builder::with_capacity(capacity)),
            Field::End => Self::End(UInt32Builder::with_capacity(capacity)),
            Field::BasesCovered => Self::BasesCovered(UInt64Builder::with_capacity(capacity)),
            Field::Min => Self::Min(Float64Builder::with_capacity(capacity)),
            Field::Max => Self::Max(Float64Builder::with_capacity(capacity)),
            Field::Sum => Self::Sum(Float64Builder::with_capacity(capacity)),
            Field::SumSquares => Self::SumSquares(Float64Builder::with_capacity(capacity)),
        };
        Ok(builder)
    }

    pub fn finish(&mut self) -> arrow::array::ArrayRef {
        match self {
            Self::Chrom(builder) => {
                let array = reset_dictarray_builder(builder);
                Arc::new(array)
            }
            Self::Start(builder) => Arc::new(builder.finish()),
            Self::End(builder) => Arc::new(builder.finish()),
            Self::BasesCovered(builder) => Arc::new(builder.finish()),
            Self::Min(builder) => Arc::new(builder.finish()),
            Self::Max(builder) => Arc::new(builder.finish()),
            Self::Sum(builder) => Arc::new(builder.finish()),
            Self::SumSquares(builder) => Arc::new(builder.finish()),
        }
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
            Field::BasesCovered,
            Field::Min,
            Field::Max,
            Field::Sum,
            Field::SumSquares,
        ] {
            let ref_names = vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()];
            let mut builder = FieldBuilder::new(field.clone(), &ref_names, 10).unwrap();
            let data_type = builder.finish().data_type().clone();
            assert_eq!(field.arrow_type(), data_type);
        }
    }

    #[test]
    fn test_field_from_str() {
        assert_eq!(Field::from_str("chrom").unwrap(), Field::Chrom);
        assert_eq!(Field::from_str("start").unwrap(), Field::Start);
        assert_eq!(Field::from_str("end").unwrap(), Field::End);
        assert_eq!(
            Field::from_str("bases_covered").unwrap(),
            Field::BasesCovered
        );
        assert_eq!(Field::from_str("min").unwrap(), Field::Min);
        assert_eq!(Field::from_str("max").unwrap(), Field::Max);
        assert_eq!(Field::from_str("sum").unwrap(), Field::Sum);
        assert_eq!(Field::from_str("sum_squares").unwrap(), Field::SumSquares);
        assert!(Field::from_str("invalid_field").is_err());
    }
}

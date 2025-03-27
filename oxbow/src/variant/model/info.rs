use std::io;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, BooleanBuilder, FixedSizeListBuilder, Float32Builder, GenericStringBuilder,
    Int32Builder, ListBuilder,
};
use arrow::datatypes::{DataType, Field as ArrowField};

use noodles::vcf::header::record::value::map::info::Number;
use noodles::vcf::header::record::value::map::info::Type;
use noodles::vcf::variant::record::info::field::value::array::Array as Values;
use noodles::vcf::variant::record::info::field::Value;

/// A variant (VCF/BCF) INFO field definition.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct InfoDef {
    pub name: String,
    pub ty: InfoType,
}

impl InfoDef {
    pub fn new(name: String, number: &Number, ty: &Type) -> Self {
        let ty = InfoType::from_noodles(ty, number);
        Self { name, ty }
    }

    pub fn get_arrow_field(&self) -> ArrowField {
        ArrowField::new(&self.name, self.ty.arrow_type(), true)
    }
}

impl TryFrom<(String, String, String)> for InfoDef {
    type Error = io::Error;

    fn try_from(def: (String, String, String)) -> Result<Self, Self::Error> {
        let (name, number, ty) = def;
        let number = match number.as_str() {
            "A" => Number::AlternateBases,
            "R" => Number::ReferenceAlternateBases,
            "G" => Number::Samples,
            "." => Number::Unknown,
            _ => {
                if let Ok(n) = number.parse::<usize>() {
                    Number::Count(n)
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Invalid number parameter for INFO field '{}': {}",
                            name, number
                        ),
                    ));
                }
            }
        };
        let ty: Type = ty.parse().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid type for INFO field '{}': {}", name, ty),
            )
        })?;
        Ok(Self::new(name, &number, &ty))
    }
}

/// A mapping of native INFO types to Arrow data types.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum InfoType {
    Character,
    CharacterList,
    CharacterFixedSizeList(usize),
    String,
    StringList,
    StringFixedSizeList(usize),
    Integer,
    IntegerList,
    IntegerFixedSizeList(usize),
    Float,
    FloatList,
    FloatFixedSizeList(usize),
    Flag,
}

impl InfoType {
    pub fn from_noodles(ty: &Type, number: &Number) -> Self {
        match ty {
            Type::Character => {
                match number {
                    Number::Count(1) => InfoType::Character,
                    Number::Count(n) => InfoType::CharacterFixedSizeList(*n),
                    _ => InfoType::CharacterList, // "A", "R", "G", "."
                }
            }
            Type::String => {
                match number {
                    Number::Count(1) => InfoType::String,
                    Number::Count(n) => InfoType::StringFixedSizeList(*n),
                    _ => InfoType::StringList, // "A", "R", "G", "."
                }
            }
            Type::Integer => {
                match number {
                    Number::Count(1) => InfoType::Integer,
                    Number::Count(n) => InfoType::IntegerFixedSizeList(*n),
                    _ => InfoType::IntegerList, // "A", "R", "G", "."
                }
            }
            Type::Float => {
                match number {
                    Number::Count(1) => InfoType::Float,
                    Number::Count(n) => InfoType::FloatFixedSizeList(*n),
                    _ => InfoType::FloatList, // "A", "R", "G", "."
                }
            }
            Type::Flag => InfoType::Flag,
        }
    }

    pub fn number(&self) -> Option<usize> {
        match self {
            Self::Flag => Some(0),
            Self::Character => Some(1),
            Self::String => Some(1),
            Self::Integer => Some(1),
            Self::Float => Some(1),
            Self::CharacterFixedSizeList(n) => Some(*n),
            Self::StringFixedSizeList(n) => Some(*n),
            Self::IntegerFixedSizeList(n) => Some(*n),
            Self::FloatFixedSizeList(n) => Some(*n),
            Self::CharacterList => None,
            Self::StringList => None,
            Self::IntegerList => None,
            Self::FloatList => None,
        }
    }

    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::Character => DataType::Utf8,
            Self::CharacterList => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                DataType::List(Arc::new(item))
            }
            Self::CharacterFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }

            Self::String => DataType::Utf8,
            Self::StringList => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                DataType::List(Arc::new(item))
            }
            Self::StringFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }

            Self::Integer => DataType::Int32,
            Self::IntegerList => {
                let item = ArrowField::new("item", DataType::Int32, true);
                DataType::List(Arc::new(item))
            }
            Self::IntegerFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::Int32, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }

            Self::Float => DataType::Float32,
            Self::FloatList => {
                let item = ArrowField::new("item", DataType::Float32, true);
                DataType::List(Arc::new(item))
            }
            Self::FloatFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::Float32, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }

            Self::Flag => DataType::Boolean,
        }
    }
}

/// A builder for Arrow arrays (columns) based on variant INFO fields.
pub enum InfoBuilder {
    Character(GenericStringBuilder<i32>),
    CharacterList(ListBuilder<GenericStringBuilder<i32>>),
    CharacterFixedSizeList(FixedSizeListBuilder<GenericStringBuilder<i32>>),
    String(GenericStringBuilder<i32>),
    StringList(ListBuilder<GenericStringBuilder<i32>>),
    StringFixedSizeList(FixedSizeListBuilder<GenericStringBuilder<i32>>),
    Integer(Int32Builder),
    IntegerList(ListBuilder<Int32Builder>),
    IntegerFixedSizeList(FixedSizeListBuilder<Int32Builder>),
    Float(Float32Builder),
    FloatList(ListBuilder<Float32Builder>),
    FloatFixedSizeList(FixedSizeListBuilder<Float32Builder>),
    Flag(BooleanBuilder),
}

impl InfoBuilder {
    pub fn new(ty: &InfoType) -> Self {
        match ty {
            InfoType::Character => InfoBuilder::Character(GenericStringBuilder::<i32>::new()),
            InfoType::CharacterList => {
                InfoBuilder::CharacterList(ListBuilder::new(GenericStringBuilder::<i32>::new()))
            }
            InfoType::CharacterFixedSizeList(n) => InfoBuilder::CharacterFixedSizeList(
                FixedSizeListBuilder::new(GenericStringBuilder::<i32>::new(), *n as i32),
            ),

            InfoType::String => InfoBuilder::String(GenericStringBuilder::<i32>::new()),
            InfoType::StringList => {
                InfoBuilder::StringList(ListBuilder::new(GenericStringBuilder::<i32>::new()))
            }
            InfoType::StringFixedSizeList(n) => InfoBuilder::StringFixedSizeList(
                FixedSizeListBuilder::new(GenericStringBuilder::<i32>::new(), *n as i32),
            ),

            InfoType::Integer => InfoBuilder::Integer(Int32Builder::new()),
            InfoType::IntegerList => {
                InfoBuilder::IntegerList(ListBuilder::new(Int32Builder::new()))
            }
            InfoType::IntegerFixedSizeList(n) => InfoBuilder::IntegerFixedSizeList(
                FixedSizeListBuilder::new(Int32Builder::new(), *n as i32),
            ),

            InfoType::Float => InfoBuilder::Float(Float32Builder::new()),
            InfoType::FloatList => InfoBuilder::FloatList(ListBuilder::new(Float32Builder::new())),
            InfoType::FloatFixedSizeList(n) => InfoBuilder::FloatFixedSizeList(
                FixedSizeListBuilder::new(Float32Builder::new(), *n as i32),
            ),

            InfoType::Flag => InfoBuilder::Flag(BooleanBuilder::new()),
        }
    }

    pub fn append_null(&mut self) {
        match self {
            Self::Character(builder) => builder.append_null(),
            Self::CharacterList(builder) => builder.append(false),
            Self::CharacterFixedSizeList(builder) => {
                for _ in 0..builder.value_length() {
                    builder.values().append_null();
                }
                builder.append(false);
            }
            Self::String(builder) => builder.append_null(),
            Self::StringList(builder) => builder.append(false),
            Self::StringFixedSizeList(builder) => {
                for _ in 0..builder.value_length() {
                    builder.values().append_null();
                }
                builder.append(false);
            }
            Self::Integer(builder) => builder.append_null(),
            Self::IntegerList(builder) => builder.append(false),
            Self::IntegerFixedSizeList(builder) => {
                for _ in 0..builder.value_length() {
                    builder.values().append_null();
                }
                builder.append(false);
            }
            Self::Float(builder) => builder.append_null(),
            Self::FloatList(builder) => builder.append(false),
            Self::FloatFixedSizeList(builder) => {
                for _ in 0..builder.value_length() {
                    builder.values().append_null();
                }
                builder.append(false);
            }
            Self::Flag(builder) => builder.append_null(),
        }
    }

    pub fn append_value(&mut self, value: &Value) -> io::Result<()> {
        match value {
            Value::Character(c) => match self {
                Self::Character(builder) => {
                    builder.append_value(c.to_string());
                }
                Self::CharacterList(builder) => {
                    builder.values().append_value(c.to_string());
                    builder.append(true);
                }
                Self::CharacterFixedSizeList(builder) => {
                    builder.values().append_value(c.to_string());
                    builder.append(true);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing INFO field: type mismatch",
                    ));
                }
            },
            Value::String(s) => match self {
                Self::String(builder) => {
                    builder.append_value(s);
                }
                Self::StringList(builder) => {
                    builder.values().append_value(s);
                    builder.append(true);
                }
                Self::StringFixedSizeList(builder) => {
                    builder.values().append_value(s);
                    builder.append(true);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing INFO field: type mismatch",
                    ));
                }
            },
            Value::Integer(n) => match self {
                Self::Integer(builder) => {
                    builder.append_value(*n);
                }
                Self::IntegerList(builder) => {
                    builder.values().append_value(*n);
                    builder.append(true);
                }
                Self::IntegerFixedSizeList(builder) => {
                    builder.values().append_value(*n);
                    builder.append(true);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing INFO field: type mismatch",
                    ));
                }
            },
            Value::Float(f) => match self {
                Self::Float(builder) => {
                    builder.append_value(*f);
                }
                Self::FloatList(builder) => {
                    builder.values().append_value(*f);
                    builder.append(true);
                }
                Self::FloatFixedSizeList(builder) => {
                    builder.values().append_value(*f);
                    builder.append(true);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing INFO field: type mismatch",
                    ));
                }
            },
            Value::Flag => match self {
                Self::Flag(builder) => {
                    builder.append_value(true);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing INFO field: type mismatch",
                    ));
                }
            },
            Value::Array(array) => {
                self.append_values(array)?;
            }
        };
        Ok(())
    }

    fn append_values(&mut self, array: &Values) -> io::Result<()> {
        match array {
            Values::Character(values) => match self {
                Self::CharacterList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(c) => {
                                non_null = true;
                                builder.values().append_value(c.to_string())
                            }
                            None => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                Self::CharacterFixedSizeList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(c) => {
                                non_null = true;
                                builder.values().append_value(c.to_string())
                            }
                            None => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing INFO field: type mismatch",
                    ));
                }
            },
            Values::String(values) => match self {
                Self::StringList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(s) => {
                                non_null = true;
                                builder.values().append_value(s)
                            }
                            None => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                Self::StringFixedSizeList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(s) => {
                                non_null = true;
                                builder.values().append_value(s)
                            }
                            None => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing INFO field: type mismatch",
                    ));
                }
            },
            Values::Integer(values) => match self {
                Self::IntegerList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(n) => {
                                non_null = true;
                                builder.values().append_value(n);
                            }
                            _ => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                Self::IntegerFixedSizeList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(n) => {
                                non_null = true;
                                builder.values().append_value(n);
                            }
                            _ => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing INFO field: type mismatch",
                    ));
                }
            },
            Values::Float(values) => match self {
                Self::FloatList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(f) => {
                                non_null = true;
                                builder.values().append_value(f);
                            }
                            _ => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                Self::FloatFixedSizeList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(f) => {
                                non_null = true;
                                builder.values().append_value(f);
                            }
                            _ => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing INFO field: type mismatch",
                    ));
                }
            },
        };
        Ok(())
    }

    pub fn finish(&mut self) -> ArrayRef {
        match self {
            Self::Character(builder) => Arc::new(builder.finish()),
            Self::CharacterList(builder) => Arc::new(builder.finish()),
            Self::CharacterFixedSizeList(builder) => Arc::new(builder.finish()),
            Self::String(builder) => Arc::new(builder.finish()),
            Self::StringList(builder) => Arc::new(builder.finish()),
            Self::StringFixedSizeList(builder) => Arc::new(builder.finish()),
            Self::Integer(builder) => Arc::new(builder.finish()),
            Self::IntegerList(builder) => Arc::new(builder.finish()),
            Self::IntegerFixedSizeList(builder) => Arc::new(builder.finish()),
            Self::Float(builder) => Arc::new(builder.finish()),
            Self::FloatList(builder) => Arc::new(builder.finish()),
            Self::FloatFixedSizeList(builder) => Arc::new(builder.finish()),
            Self::Flag(builder) => Arc::new(builder.finish()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::vcf::header::record::value::map::info::{Number, Type};

    #[test]
    fn test_infodef_new() {
        let name = "DP".to_string();
        let number = Number::Count(1);
        let ty = Type::Integer;

        let infodef = InfoDef::new(name.clone(), &number, &ty);

        assert_eq!(infodef.name, name);
        assert_eq!(infodef.ty, InfoType::Integer);
    }

    #[test]
    fn test_infodef_try_from_valid() {
        let def = ("DP".to_string(), "1".to_string(), "Integer".to_string());
        let infodef = InfoDef::try_from(def).unwrap();

        assert_eq!(infodef.name, "DP");
        assert_eq!(infodef.ty, InfoType::Integer);
    }

    #[test]
    fn test_infodef_try_from_invalid_number() {
        let def = (
            "DP".to_string(),
            "invalid".to_string(),
            "Integer".to_string(),
        );
        let result = InfoDef::try_from(def);
        assert!(result.is_err());

        let def = ("DP".to_string(), "1".to_string(), "InvalidType".to_string());
        let result = InfoDef::try_from(def);
        assert!(result.is_err());
    }

    #[test]
    fn test_infobuilder_arrow_type() {
        for ty in [
            InfoType::Character,
            InfoType::CharacterList,
            InfoType::CharacterFixedSizeList(10),
            InfoType::String,
            InfoType::StringList,
            InfoType::StringFixedSizeList(10),
            InfoType::Integer,
            InfoType::IntegerList,
            InfoType::IntegerFixedSizeList(10),
            InfoType::Float,
            InfoType::FloatList,
            InfoType::FloatFixedSizeList(10),
            InfoType::Flag,
        ] {
            let mut builder = InfoBuilder::new(&ty);
            let data_type = builder.finish().data_type().clone();
            assert_eq!(ty.arrow_type(), data_type);
        }
    }

    #[test]
    fn test_infobuilder_append_null() {
        let mut builder = InfoBuilder::new(&InfoType::Integer);
        builder.append_null();
        let array = builder.finish();
        assert!(array.is_nullable());
        assert!(array.is_null(0));
    }

    #[test]
    fn test_infobuilder_append_value() {
        let mut builder = InfoBuilder::new(&InfoType::Integer);
        let value = Value::Integer(42);
        builder.append_value(&value).unwrap();
        let array = builder.finish();
        let int_array = array
            .as_any()
            .downcast_ref::<arrow::array::Int32Array>()
            .unwrap();
        assert_eq!(int_array.len(), 1);
        assert_eq!(int_array.value(0), 42);
    }
}

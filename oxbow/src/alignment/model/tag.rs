use std::collections::BTreeMap;
use std::io;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, Float32Builder, GenericStringBuilder, Int16Builder, Int32Builder, Int8Builder,
    ListBuilder, UInt16Builder, UInt32Builder, UInt8Builder,
};
use arrow::datatypes::{DataType, Field as ArrowField};

use bstr::ByteSlice;
use noodles::sam::alignment::record::data::field::value::array::Subtype;
use noodles::sam::alignment::record::data::field::value::Array;
use noodles::sam::alignment::record::data::field::{Tag, Type, Value};
use noodles::sam::alignment::record::Data;

/// An alignment (SAM/BAM/CRAM) tag definition.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct TagDef {
    pub name: String,
    pub ty: TagType,
}

impl TagDef {
    pub fn into_bytes(&self) -> [u8; 2] {
        self.name.as_bytes().try_into().unwrap()
    }

    pub fn get_arrow_field(&self) -> ArrowField {
        ArrowField::new(&self.name, self.ty.arrow_type(), true)
    }
}

impl TryFrom<(String, String)> for TagDef {
    type Error = io::Error;

    fn try_from(def: (String, String)) -> Result<Self, Self::Error> {
        let (name, code) = def;
        if name.len() != 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Tag name must be 2 characters",
            ));
        }
        let ty: TagType = code.parse().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid tag type: {}", code),
            )
        })?;
        Ok(Self { name, ty })
    }
}

/// A mapping of native alignment tag types to Arrow data types.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum TagType {
    Character,
    String,
    Hex,
    Int8,
    UInt8,
    Int16,
    UInt16,
    Int32,
    UInt32,
    Float,
    ArrayInt8,
    ArrayUInt8,
    ArrayInt16,
    ArrayUInt16,
    ArrayInt32,
    ArrayUInt32,
    ArrayFloat,
}

impl TagType {
    pub fn code(&self) -> &str {
        match self {
            Self::Character => "A",
            Self::String => "Z",
            Self::Hex => "H",
            Self::Int8 => "c",
            Self::UInt8 => "C",
            Self::Int16 => "s",
            Self::UInt16 => "S",
            Self::Int32 => "i",
            Self::UInt32 => "I",
            Self::Float => "f",
            Self::ArrayInt8 => "Bc",
            Self::ArrayUInt8 => "BC",
            Self::ArrayInt16 => "Bs",
            Self::ArrayUInt16 => "BS",
            Self::ArrayInt32 => "Bi",
            Self::ArrayUInt32 => "BI",
            Self::ArrayFloat => "Bf",
        }
    }

    pub fn noodles_type(&self) -> Type {
        match self {
            Self::Character => Type::Character,
            Self::String => Type::String,
            Self::Hex => Type::Hex,
            Self::Int8 => Type::Int8,
            Self::UInt8 => Type::UInt8,
            Self::Int16 => Type::Int16,
            Self::UInt16 => Type::UInt16,
            Self::Int32 => Type::Int32,
            Self::UInt32 => Type::UInt32,
            Self::Float => Type::Float,
            Self::ArrayInt8 => Type::Array,
            Self::ArrayUInt8 => Type::Array,
            Self::ArrayInt16 => Type::Array,
            Self::ArrayUInt16 => Type::Array,
            Self::ArrayInt32 => Type::Array,
            Self::ArrayUInt32 => Type::Array,
            Self::ArrayFloat => Type::Array,
        }
    }

    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::Character => DataType::Utf8,
            Self::String => DataType::Utf8,
            Self::Hex => DataType::Utf8,
            Self::Int8 => DataType::Int8,
            Self::UInt8 => DataType::UInt8,
            Self::Int16 => DataType::Int16,
            Self::UInt16 => DataType::UInt16,
            Self::Int32 => DataType::Int32,
            Self::UInt32 => DataType::UInt32,
            Self::Float => DataType::Float32,
            Self::ArrayInt8 => {
                let item = ArrowField::new("item", DataType::Int8, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayUInt8 => {
                let item = ArrowField::new("item", DataType::UInt8, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayInt16 => {
                let item = ArrowField::new("item", DataType::Int16, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayUInt16 => {
                let item = ArrowField::new("item", DataType::UInt16, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayInt32 => {
                let item = ArrowField::new("item", DataType::Int32, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayUInt32 => {
                let item = ArrowField::new("item", DataType::UInt32, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayFloat => {
                let item = ArrowField::new("item", DataType::Float32, true);
                DataType::List(Arc::new(item))
            }
        }
    }
}

impl From<&Value<'_>> for TagType {
    fn from(value: &Value) -> Self {
        match value {
            Value::Character(_) => TagType::Character,
            Value::String(_) => TagType::String,
            Value::Hex(_) => TagType::Hex,
            Value::Int8(_) => TagType::Int8,
            Value::UInt8(_) => TagType::UInt8,
            Value::Int16(_) => TagType::Int16,
            Value::UInt16(_) => TagType::UInt16,
            Value::Int32(_) => TagType::Int32,
            Value::UInt32(_) => TagType::UInt32,
            Value::Float(_) => TagType::Float,
            Value::Array(array) => match array.subtype() {
                Subtype::Int8 => TagType::ArrayInt8,
                Subtype::UInt8 => TagType::ArrayUInt8,
                Subtype::Int16 => TagType::ArrayInt16,
                Subtype::UInt16 => TagType::ArrayUInt16,
                Subtype::Int32 => TagType::ArrayInt32,
                Subtype::UInt32 => TagType::ArrayUInt32,
                Subtype::Float => TagType::ArrayFloat,
            },
        }
    }
}

impl FromStr for TagType {
    type Err = std::io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" => Ok(TagType::Character),
            "Z" => Ok(TagType::String),
            "H" => Ok(TagType::Hex),
            "c" => Ok(TagType::Int8),
            "C" => Ok(TagType::UInt8),
            "s" => Ok(TagType::Int16),
            "S" => Ok(TagType::UInt16),
            "i" => Ok(TagType::Int32),
            "I" => Ok(TagType::UInt32),
            "f" => Ok(TagType::Float),
            "Bc" => Ok(TagType::ArrayInt8),
            "BC" => Ok(TagType::ArrayUInt8),
            "Bs" => Ok(TagType::ArrayInt16),
            "BS" => Ok(TagType::ArrayUInt16),
            "Bi" => Ok(TagType::ArrayInt32),
            "BI" => Ok(TagType::ArrayUInt32),
            "Bf" => Ok(TagType::ArrayFloat),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Invalid tag type code: {}", s),
            )),
        }
    }
}

/// A builder for Arrow arrays (columns) based on alignment tags.
#[derive(Debug)]
pub enum TagBuilder {
    Character(GenericStringBuilder<i32>),
    String(GenericStringBuilder<i32>),
    Hex(GenericStringBuilder<i32>),
    Int8(Int8Builder),
    UInt8(UInt8Builder),
    Int16(Int16Builder),
    UInt16(UInt16Builder),
    Int32(Int32Builder),
    UInt32(UInt32Builder),
    Float(Float32Builder),
    ArrayInt8(ListBuilder<Int8Builder>),
    ArrayUInt8(ListBuilder<UInt8Builder>),
    ArrayInt16(ListBuilder<Int16Builder>),
    ArrayUInt16(ListBuilder<UInt16Builder>),
    ArrayInt32(ListBuilder<Int32Builder>),
    ArrayUInt32(ListBuilder<UInt32Builder>),
    ArrayFloat(ListBuilder<Float32Builder>),
}

impl TagBuilder {
    pub fn new(tag_type: &TagType) -> Self {
        match tag_type {
            TagType::Character => Self::Character(GenericStringBuilder::<i32>::new()),
            TagType::String => Self::String(GenericStringBuilder::<i32>::new()),
            TagType::Hex => Self::Hex(GenericStringBuilder::<i32>::new()),
            TagType::Int8 => Self::Int8(Int8Builder::new()),
            TagType::UInt8 => Self::UInt8(UInt8Builder::new()),
            TagType::Int16 => Self::Int16(Int16Builder::new()),
            TagType::UInt16 => Self::UInt16(UInt16Builder::new()),
            TagType::Int32 => Self::Int32(Int32Builder::new()),
            TagType::UInt32 => Self::UInt32(UInt32Builder::new()),
            TagType::Float => Self::Float(Float32Builder::new()),
            TagType::ArrayInt8 => Self::ArrayInt8(ListBuilder::new(Int8Builder::new())),
            TagType::ArrayUInt8 => Self::ArrayUInt8(ListBuilder::new(UInt8Builder::new())),
            TagType::ArrayInt16 => Self::ArrayInt16(ListBuilder::new(Int16Builder::new())),
            TagType::ArrayUInt16 => Self::ArrayUInt16(ListBuilder::new(UInt16Builder::new())),
            TagType::ArrayInt32 => Self::ArrayInt32(ListBuilder::new(Int32Builder::new())),
            TagType::ArrayUInt32 => Self::ArrayUInt32(ListBuilder::new(UInt32Builder::new())),
            TagType::ArrayFloat => Self::ArrayFloat(ListBuilder::new(Float32Builder::new())),
        }
    }

    pub fn append_null(&mut self) {
        match self {
            Self::Character(builder) => builder.append_null(),
            Self::String(builder) => builder.append_null(),
            Self::Hex(builder) => builder.append_null(),
            Self::Int8(builder) => builder.append_null(),
            Self::UInt8(builder) => builder.append_null(),
            Self::Int16(builder) => builder.append_null(),
            Self::UInt16(builder) => builder.append_null(),
            Self::Int32(builder) => builder.append_null(),
            Self::UInt32(builder) => builder.append_null(),
            Self::Float(builder) => builder.append_null(),
            Self::ArrayInt8(builder) => builder.append(false),
            Self::ArrayUInt8(builder) => builder.append(false),
            Self::ArrayInt16(builder) => builder.append(false),
            Self::ArrayUInt16(builder) => builder.append(false),
            Self::ArrayInt32(builder) => builder.append(false),
            Self::ArrayUInt32(builder) => builder.append(false),
            Self::ArrayFloat(builder) => builder.append(false),
        }
    }

    pub fn append_value(&mut self, value: &Value) -> io::Result<()> {
        match value {
            Value::Character(v) => match self {
                Self::Character(builder) => {
                    builder.append_value(v.to_string());
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        value, self
                    ),
                )),
            },
            Value::String(v) => match self {
                Self::String(builder) => {
                    match v.to_str() {
                        Ok(s) => builder.append_value(s),
                        Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, e)),
                    }
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        value, self
                    ),
                )),
            },
            Value::Hex(v) => match self {
                Self::Hex(builder) => {
                    match v.to_str() {
                        Ok(s) => builder.append_value(s),
                        Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, e)),
                    }
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        value, self
                    ),
                )),
            },
            Value::Int8(v) => match self {
                Self::Int8(builder) => {
                    builder.append_value(*v);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        value, self
                    ),
                )),
            },
            Value::UInt8(v) => match self {
                Self::UInt8(builder) => {
                    builder.append_value(*v);
                    Ok(())
                }
                Self::Int8(builder) => {
                    if *v > i8::MAX as u8 {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("Value {} exceeds i8 range", v),
                        ));
                    }
                    builder.append_value(*v as i8);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        value, self
                    ),
                )),
            },
            Value::Int16(v) => match self {
                Self::Int16(builder) => {
                    builder.append_value(*v);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        value, self
                    ),
                )),
            },
            Value::UInt16(v) => match self {
                Self::UInt16(builder) => {
                    builder.append_value(*v);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        value, self
                    ),
                )),
            },
            Value::Int32(v) => match self {
                Self::Int32(builder) => {
                    builder.append_value(*v);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        value, self
                    ),
                )),
            },
            Value::UInt32(v) => match self {
                Self::UInt32(builder) => {
                    builder.append_value(*v);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        value, self
                    ),
                )),
            },
            Value::Float(v) => match self {
                Self::Float(builder) => {
                    builder.append_value(*v);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        value, self
                    ),
                )),
            },
            Value::Array(array) => self.append_values(array),
        }
    }

    fn append_values(&mut self, array: &Array) -> io::Result<()> {
        match array {
            Array::Int8(values) => match self {
                Self::ArrayInt8(builder) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(v),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        array, self
                    ),
                )),
            },
            Array::UInt8(values) => match self {
                Self::ArrayUInt8(builder) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(v),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        array, self
                    ),
                )),
            },
            Array::Int16(values) => match self {
                Self::ArrayInt16(builder) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(v),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        array, self
                    ),
                )),
            },
            Array::UInt16(values) => match self {
                Self::ArrayUInt16(builder) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(v),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        array, self
                    ),
                )),
            },
            Array::Int32(values) => match self {
                Self::ArrayInt32(builder) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(v),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        array, self
                    ),
                )),
            },
            Array::UInt32(values) => match self {
                Self::ArrayUInt32(builder) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(v),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        array, self
                    ),
                )),
            },
            Array::Float(values) => match self {
                Self::ArrayFloat(builder) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(v),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected builder for {:?}, got {:?}",
                        array, self
                    ),
                )),
            },
        }
    }

    pub fn finish(&mut self) -> ArrayRef {
        match self {
            Self::Character(builder) => Arc::new(builder.finish()),
            Self::String(builder) => Arc::new(builder.finish()),
            Self::Hex(builder) => Arc::new(builder.finish()),
            Self::Int8(builder) => Arc::new(builder.finish()),
            Self::UInt8(builder) => Arc::new(builder.finish()),
            Self::Int16(builder) => Arc::new(builder.finish()),
            Self::UInt16(builder) => Arc::new(builder.finish()),
            Self::Int32(builder) => Arc::new(builder.finish()),
            Self::UInt32(builder) => Arc::new(builder.finish()),
            Self::Float(builder) => Arc::new(builder.finish()),
            Self::ArrayInt8(builder) => Arc::new(builder.finish()),
            Self::ArrayUInt8(builder) => Arc::new(builder.finish()),
            Self::ArrayInt16(builder) => Arc::new(builder.finish()),
            Self::ArrayUInt16(builder) => Arc::new(builder.finish()),
            Self::ArrayInt32(builder) => Arc::new(builder.finish()),
            Self::ArrayUInt32(builder) => Arc::new(builder.finish()),
            Self::ArrayFloat(builder) => Arc::new(builder.finish()),
        }
    }
}

/// A scanner to collect unique tag definitions from a stream of alignment records.
pub struct TagScanner {
    tags: BTreeMap<Tag, String>,
}

impl Default for TagScanner {
    fn default() -> Self {
        Self::new()
    }
}

impl TagScanner {
    pub fn new() -> Self {
        Self {
            tags: BTreeMap::new(),
        }
    }

    pub fn push(&mut self, record: &impl noodles::sam::alignment::Record) {
        let data = record.data();
        data.iter().for_each(|result| match result {
            Ok((tag, value)) => {
                let ty = TagType::from(&value);
                let type_code = ty.code();
                self.tags
                    .entry(tag)
                    .or_insert_with(|| type_code.to_string());
            }
            Err(e) => {
                eprintln!("Error processing tag: {}", e);
            }
        });
    }

    pub fn collect(&self) -> Vec<(String, String)> {
        let mut tag_defs: Vec<(String, String)> = self
            .tags
            .iter()
            .map(|(tag, type_code)| {
                let tag_name = format!("{}", &tag.as_ref().as_bstr());
                (tag_name, type_code.clone())
            })
            .collect();
        tag_defs.sort_by(|a, b| a.0.cmp(&b.0));
        tag_defs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tag_def_try_from_valid() {
        let tag_def = TagDef::try_from(("NM".to_string(), "i".to_string())).unwrap();
        assert_eq!(tag_def.name, "NM");
        assert_eq!(tag_def.ty, TagType::Int32);
    }

    #[test]
    fn test_tag_def_try_from_invalid() {
        let result = TagDef::try_from(("N".to_string(), "i".to_string()));
        assert!(result.is_err());
        let result = TagDef::try_from(("NM".to_string(), "x".to_string()));
        assert!(result.is_err());
    }

    #[test]
    fn test_tag_def_into_bytes() {
        let tag_def = TagDef {
            name: "NM".to_string(),
            ty: TagType::Int32,
        };
        assert_eq!(tag_def.into_bytes(), [78, 77]); // ASCII for 'N' and 'M'
    }

    #[test]
    fn test_tag_arrow_type() {
        for ty in [
            TagType::Character,
            TagType::String,
            TagType::Hex,
            TagType::Int8,
            TagType::UInt8,
            TagType::Int16,
            TagType::UInt16,
            TagType::Int32,
            TagType::UInt32,
            TagType::Float,
        ] {
            let mut builder = TagBuilder::new(&ty);
            let data_type = builder.finish().data_type().clone();
            assert_eq!(ty.arrow_type(), data_type);
        }
    }

    #[test]
    fn test_tag_builder_append_null() {
        let mut builder = TagBuilder::new(&TagType::Int32);
        builder.append_null();
        let array = builder.finish();
        assert!(array.is_nullable());
        assert!(array.is_null(0));
    }

    #[test]
    fn test_tag_builder_append_value() {
        let mut builder = TagBuilder::new(&TagType::Int32);
        builder.append_value(&Value::Int32(42)).unwrap();
        let array = builder.finish();
        let int_array = array.as_any().downcast_ref::<arrow::array::Int32Array>().unwrap();
        assert_eq!(int_array.value(0), 42);
    }
}

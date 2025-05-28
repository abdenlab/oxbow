use std::collections::BTreeMap;
use std::io;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::{ArrayRef, Float32Builder, GenericStringBuilder, Int64Builder, ListBuilder};
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
            Self::Int8 => DataType::Int64,
            Self::UInt8 => DataType::Int64,
            Self::Int16 => DataType::Int64,
            Self::UInt16 => DataType::Int64,
            Self::Int32 => DataType::Int64,
            Self::UInt32 => DataType::Int64,
            Self::Float => DataType::Float32,
            Self::ArrayInt8 => {
                let item = ArrowField::new("item", DataType::Int64, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayUInt8 => {
                let item = ArrowField::new("item", DataType::Int64, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayInt16 => {
                let item = ArrowField::new("item", DataType::Int64, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayUInt16 => {
                let item = ArrowField::new("item", DataType::Int64, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayInt32 => {
                let item = ArrowField::new("item", DataType::Int64, true);
                DataType::List(Arc::new(item))
            }
            Self::ArrayUInt32 => {
                let item = ArrowField::new("item", DataType::Int64, true);
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
            Value::Character(_) => Self::Character,
            Value::String(_) => Self::String,
            Value::Hex(_) => Self::Hex,
            Value::Int8(_) => Self::Int8,
            Value::UInt8(_) => Self::UInt8,
            Value::Int16(_) => Self::Int16,
            Value::UInt16(_) => Self::UInt16,
            Value::Int32(_) => Self::Int32,
            Value::UInt32(_) => Self::UInt32,
            Value::Float(_) => Self::Float,
            Value::Array(array) => match array.subtype() {
                Subtype::Int8 => Self::ArrayInt8,
                Subtype::UInt8 => Self::ArrayUInt8,
                Subtype::Int16 => Self::ArrayInt16,
                Subtype::UInt16 => Self::ArrayUInt16,
                Subtype::Int32 => Self::ArrayInt32,
                Subtype::UInt32 => Self::ArrayUInt32,
                Subtype::Float => Self::ArrayFloat,
            },
        }
    }
}

impl FromStr for TagType {
    type Err = std::io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" => Ok(Self::Character),
            "Z" => Ok(Self::String),
            "H" => Ok(Self::Hex),
            "c" => Ok(Self::Int8),
            "C" => Ok(Self::UInt8),
            "s" => Ok(Self::Int16),
            "S" => Ok(Self::UInt16),
            "i" => Ok(Self::Int32),
            "I" => Ok(Self::UInt32),
            "f" => Ok(Self::Float),
            "Bc" => Ok(Self::ArrayInt8),
            "BC" => Ok(Self::ArrayUInt8),
            "Bs" => Ok(Self::ArrayInt16),
            "BS" => Ok(Self::ArrayUInt16),
            "Bi" => Ok(Self::ArrayInt32),
            "BI" => Ok(Self::ArrayUInt32),
            "Bf" => Ok(Self::ArrayFloat),
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("Invalid tag type code: {}", s),
            )),
        }
    }
}

/// A builder for Arrow arrays (columns) based on alignment tags.
///
/// The `String` builder also serves as a sink for any tag value in cases where there are
/// unresolvable mismatches between tag types for the same tag name across different records.
#[derive(Debug)]
pub enum TagBuilder {
    Character(GenericStringBuilder<i32>),
    String(GenericStringBuilder<i32>),
    Hex(GenericStringBuilder<i32>),
    Integer(Int64Builder),
    Float(Float32Builder),
    IntegerArray(ListBuilder<Int64Builder>),
    FloatArray(ListBuilder<Float32Builder>),
}

impl TagBuilder {
    pub fn new(tag_type: &TagType) -> Self {
        match tag_type {
            TagType::Character => Self::Character(GenericStringBuilder::<i32>::new()),
            TagType::String => Self::String(GenericStringBuilder::<i32>::new()),
            TagType::Hex => Self::Hex(GenericStringBuilder::<i32>::new()),
            TagType::Int8 => Self::Integer(Int64Builder::new()),
            TagType::UInt8 => Self::Integer(Int64Builder::new()),
            TagType::Int16 => Self::Integer(Int64Builder::new()),
            TagType::UInt16 => Self::Integer(Int64Builder::new()),
            TagType::Int32 => Self::Integer(Int64Builder::new()),
            TagType::UInt32 => Self::Integer(Int64Builder::new()),
            TagType::Float => Self::Float(Float32Builder::new()),
            TagType::ArrayInt8 => Self::IntegerArray(ListBuilder::new(Int64Builder::new())),
            TagType::ArrayUInt8 => Self::IntegerArray(ListBuilder::new(Int64Builder::new())),
            TagType::ArrayInt16 => Self::IntegerArray(ListBuilder::new(Int64Builder::new())),
            TagType::ArrayUInt16 => Self::IntegerArray(ListBuilder::new(Int64Builder::new())),
            TagType::ArrayInt32 => Self::IntegerArray(ListBuilder::new(Int64Builder::new())),
            TagType::ArrayUInt32 => Self::IntegerArray(ListBuilder::new(Int64Builder::new())),
            TagType::ArrayFloat => Self::FloatArray(ListBuilder::new(Float32Builder::new())),
        }
    }

    pub fn append_null(&mut self) {
        match self {
            Self::Character(builder) => builder.append_null(),
            Self::String(builder) => builder.append_null(),
            Self::Hex(builder) => builder.append_null(),
            Self::Integer(builder) => builder.append_null(),
            Self::Float(builder) => builder.append_null(),
            Self::IntegerArray(builder) => builder.append(false),
            Self::FloatArray(builder) => builder.append(false),
        }
    }

    pub fn append_value(&mut self, value: &Value) -> io::Result<()> {
        match self {
            Self::Character(builder) => {
                if let Value::Character(v) = value {
                    builder.append_value(v.to_string());
                    Ok(())
                } else {
                    Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("Type mismatch: expected Character, got {:?}", value),
                    ))
                }
            }
            Self::Hex(builder) => {
                if let Value::Hex(v) = value {
                    match v.to_str() {
                        Ok(s) => {
                            builder.append_value(s);
                            Ok(())
                        }
                        Err(e) => Err(io::Error::new(io::ErrorKind::InvalidData, e)),
                    }
                } else {
                    Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("Type mismatch: expected Hex, got {:?}", value),
                    ))
                }
            }
            Self::Integer(builder) => match value.as_int() {
                Some(v) => {
                    builder.append_value(v);
                    Ok(())
                }
                None => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Type mismatch: expected Integer, got {:?}", value),
                )),
            },
            Self::Float(builder) => match value {
                Value::Float(v) => {
                    builder.append_value(*v);
                    Ok(())
                }
                Value::Int8(v) => {
                    builder.append_value(f32::from(*v));
                    Ok(())
                }
                Value::UInt8(v) => {
                    builder.append_value(f32::from(*v));
                    Ok(())
                }
                Value::Int16(v) => {
                    builder.append_value(f32::from(*v));
                    Ok(())
                }
                Value::UInt16(v) => {
                    builder.append_value(f32::from(*v));
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected Float or compatible integer, got {:?}",
                        value
                    ),
                )),
            },
            Self::String(builder) => match value {
                Value::String(v) => match v.to_str() {
                    Ok(s) => {
                        builder.append_value(s);
                        Ok(())
                    }
                    Err(e) => Err(io::Error::new(io::ErrorKind::InvalidData, e)),
                },
                Value::Character(v) => {
                    builder.append_value(v.to_string());
                    Ok(())
                }
                Value::Hex(v) => {
                    builder.append_value(v.to_string());
                    Ok(())
                }
                Value::Int8(v) => {
                    builder.append_value(v.to_string());
                    Ok(())
                }
                Value::UInt8(v) => {
                    builder.append_value(v.to_string());
                    Ok(())
                }
                Value::Int16(v) => {
                    builder.append_value(v.to_string());
                    Ok(())
                }
                Value::UInt16(v) => {
                    builder.append_value(v.to_string());
                    Ok(())
                }
                Value::Int32(v) => {
                    builder.append_value(v.to_string());
                    Ok(())
                }
                Value::UInt32(v) => {
                    builder.append_value(v.to_string());
                    Ok(())
                }
                Value::Float(v) => {
                    builder.append_value(v.to_string());
                    Ok(())
                }
                Value::Array(array) => {
                    let s = array_to_string(array)?;
                    builder.append_value(s);
                    Ok(())
                }
            },
            _ => {
                if let Value::Array(array) = value {
                    self.append_values(array)?;
                    Ok(())
                } else {
                    Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Type mismatch: expected an Array, got non-array {:?}",
                            value
                        ),
                    ))
                }
            }
        }
    }

    fn append_values(&mut self, array: &Array) -> io::Result<()> {
        match self {
            Self::IntegerArray(builder) => match array {
                Array::Int8(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(i64::from(v)),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                Array::UInt8(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(i64::from(v)),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                Array::Int16(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(i64::from(v)),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                Array::UInt16(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(i64::from(v)),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                Array::Int32(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(i64::from(v)),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                Array::UInt32(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(i64::from(v)),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected Array<Int8/16/32/UInt8/16/32>, got {:?}",
                        array
                    ),
                )),
            },
            Self::FloatArray(builder) => match array {
                Array::Float(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(v),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                Array::Int8(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(f32::from(v)),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                Array::UInt8(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(f32::from(v)),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                Array::Int16(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(f32::from(v)),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                Array::UInt16(values) => {
                    for value in values.iter() {
                        match value {
                            Ok(v) => builder.values().append_value(f32::from(v)),
                            Err(e) => return Err(e),
                        }
                    }
                    builder.append(true);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected Array<Float> or compatible integers, got {:?}",
                        array
                    ),
                )),
            },
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Type mismatch: expected an Array, got non-array {:?}",
                    array,
                ),
            )),
        }
    }

    pub fn finish(&mut self) -> ArrayRef {
        match self {
            Self::Character(builder) => Arc::new(builder.finish()),
            Self::String(builder) => Arc::new(builder.finish()),
            Self::Hex(builder) => Arc::new(builder.finish()),
            Self::Integer(builder) => Arc::new(builder.finish()),
            Self::Float(builder) => Arc::new(builder.finish()),
            Self::IntegerArray(builder) => Arc::new(builder.finish()),
            Self::FloatArray(builder) => Arc::new(builder.finish()),
        }
    }
}

fn array_to_string(array: &Array) -> io::Result<String> {
    match array {
        Array::Int8(values) => {
            let values: Result<Vec<String>, _> =
                values.iter().map(|v| v.map(|x| x.to_string())).collect();
            let values = values.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(format!("[{}]", values.join(", ")))
        }
        Array::UInt8(values) => {
            let values: Result<Vec<String>, _> =
                values.iter().map(|v| v.map(|x| x.to_string())).collect();
            let values = values.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(format!("[{}]", values.join(", ")))
        }
        Array::Int16(values) => {
            let values: Result<Vec<String>, _> =
                values.iter().map(|v| v.map(|x| x.to_string())).collect();
            let values = values.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(format!("[{}]", values.join(", ")))
        }
        Array::UInt16(values) => {
            let values: Result<Vec<String>, _> =
                values.iter().map(|v| v.map(|x| x.to_string())).collect();
            let values = values.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(format!("[{}]", values.join(", ")))
        }
        Array::Int32(values) => {
            let values: Result<Vec<String>, _> =
                values.iter().map(|v| v.map(|x| x.to_string())).collect();
            let values = values.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(format!("[{}]", values.join(", ")))
        }
        Array::UInt32(values) => {
            let values: Result<Vec<String>, _> =
                values.iter().map(|v| v.map(|x| x.to_string())).collect();
            let values = values.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(format!("[{}]", values.join(", ")))
        }
        Array::Float(values) => {
            let values: Result<Vec<String>, _> =
                values.iter().map(|v| v.map(|x| x.to_string())).collect();
            let values = values.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            Ok(format!("[{}]", values.join(", ")))
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
    fn test_tagdef_try_from_valid() {
        let tag_def = TagDef::try_from(("NM".to_string(), "i".to_string())).unwrap();
        assert_eq!(tag_def.name, "NM");
        assert_eq!(tag_def.ty, TagType::Int32);
    }

    #[test]
    fn test_tagdef_try_from_invalid() {
        let result = TagDef::try_from(("N".to_string(), "i".to_string()));
        assert!(result.is_err());
        let result = TagDef::try_from(("NM".to_string(), "x".to_string()));
        assert!(result.is_err());
    }

    #[test]
    fn test_tagdef_into_bytes() {
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
    fn test_tagbuilder_append_null() {
        let mut builder = TagBuilder::new(&TagType::Int32);
        builder.append_null();
        let array = builder.finish();
        assert!(array.is_nullable());
        assert!(array.is_null(0));
    }

    #[test]
    fn test_tagbuilder_append_value() {
        let mut builder = TagBuilder::new(&TagType::Int32);
        builder.append_value(&Value::Int32(42)).unwrap();
        let array = builder.finish();
        let int_array = array
            .as_any()
            .downcast_ref::<arrow::array::Int64Array>()
            .unwrap();
        assert_eq!(int_array.len(), 1);
        assert_eq!(int_array.value(0), 42);
    }
}

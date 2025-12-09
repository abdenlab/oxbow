use std::collections::BTreeMap;
use std::io;
use std::sync::Arc;

use arrow::array::{ArrayRef, GenericStringBuilder, ListBuilder};
use arrow::datatypes::{DataType, Field as ArrowField};
use noodles::gff::feature::record::attributes::field::Value as FeatureAttributeValue;
use noodles::gff::feature::record::Attributes as FeatureAttributes;
use noodles::gff::record::attributes::field::Value as GffAttributeValue;
use noodles::gtf::record::Attributes as GtfAttributes;

/// An GXF attribute definition.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct AttributeDef {
    pub name: String,
    pub ty: AttributeType,
}

impl AttributeDef {
    pub fn get_arrow_field(&self) -> ArrowField {
        ArrowField::new(&self.name, self.ty.arrow_type(), true)
    }
}

impl TryFrom<(String, String)> for AttributeDef {
    type Error = io::Error;

    fn try_from(def: (String, String)) -> Result<Self, Self::Error> {
        let (name, ty) = def;
        let ty = match ty.to_lowercase().as_str() {
            "string" => Ok(AttributeType::String),
            "array" => Ok(AttributeType::Array),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Invalid attribute type: '{}'. Must be 'String' or 'Array'.",
                    ty
                ),
            )),
        }?;
        Ok(Self { name, ty })
    }
}

/// A mapping of native attribute field types to Arrow data types.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum AttributeType {
    String,
    Array,
}

impl AttributeType {
    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::String => DataType::Utf8,
            Self::Array => DataType::List(Arc::new(ArrowField::new("item", DataType::Utf8, true))),
        }
    }
}

impl std::fmt::Display for AttributeType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::String => "String",
                Self::Array => "Array",
            }
        )
    }
}

/// Harmonizes GTF and GFF attribute values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AttributeValue {
    String(String),
    Array(Vec<String>),
}

impl From<FeatureAttributeValue<'_>> for AttributeValue {
    fn from(value: FeatureAttributeValue) -> Self {
        match value {
            FeatureAttributeValue::String(s) => Self::String(s.to_string()),
            FeatureAttributeValue::Array(a) => {
                Self::Array(a.iter().map(|s| s.unwrap().to_string()).collect())
            }
        }
    }
}

impl<'a> From<&'a GffAttributeValue<'a>> for AttributeValue {
    fn from(value: &'a GffAttributeValue<'a>) -> Self {
        match value {
            GffAttributeValue::String(s) => Self::String(s.to_string()),
            GffAttributeValue::Array(a) => Self::Array(a.iter().map(|s| s.to_string()).collect()),
        }
    }
}

/// A builder for Arrow arrays (columns) based on GXF attributes.
#[derive(Debug)]
pub enum AttributeBuilder {
    String(GenericStringBuilder<i32>),
    Array(ListBuilder<GenericStringBuilder<i32>>),
}

impl AttributeBuilder {
    pub fn new(attr_type: &AttributeType) -> Self {
        match attr_type {
            AttributeType::String => Self::String(GenericStringBuilder::<i32>::new()),
            AttributeType::Array => Self::Array(ListBuilder::<GenericStringBuilder<i32>>::new(
                GenericStringBuilder::<i32>::new(),
            )),
        }
    }

    pub fn append_null(&mut self) {
        match self {
            Self::String(builder) => builder.append_null(),
            Self::Array(builder) => builder.append(false),
        }
    }

    pub fn append_value(&mut self, value: &AttributeValue) -> io::Result<()> {
        match self {
            Self::String(builder) => match value {
                AttributeValue::String(v) => {
                    builder.append_value(v);
                    Ok(())
                }
                _ => Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "Type mismatch: expected an attribute of type String, got {:?}",
                        value
                    ),
                )),
            },
            Self::Array(builder) => match value {
                AttributeValue::Array(array) => {
                    for value in array.iter() {
                        builder.values().append_value(value);
                    }
                    builder.append(true);
                    Ok(())
                }
                AttributeValue::String(v) => {
                    // interpret a scalar as a length-1 array
                    builder.values().append_value(v);
                    builder.append(true);
                    Ok(())
                }
            },
        }
    }

    pub fn finish(&mut self) -> ArrayRef {
        match self {
            Self::String(builder) => Arc::new(builder.finish()),
            Self::Array(builder) => Arc::new(builder.finish()),
        }
    }
}

/// A scanner to collect unique attribute definitions from a stream of GXF records.
pub struct AttributeScanner {
    attrs: BTreeMap<String, AttributeType>,
}

impl Default for AttributeScanner {
    fn default() -> Self {
        Self::new()
    }
}

impl AttributeScanner {
    pub fn new() -> Self {
        Self {
            attrs: BTreeMap::new(),
        }
    }

    pub fn collect(&self) -> Vec<(String, String)> {
        let defs: Vec<(String, String)> = self
            .attrs
            .iter()
            .map(|(key, attr_type)| (key.clone(), attr_type.to_string()))
            .collect();
        defs
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T);
}

impl Push<noodles::gtf::Record<'_>> for AttributeScanner {
    fn push(&mut self, record: noodles::gtf::Record) {
        let attrs = match record.attributes() {
            Ok(attrs) => attrs,
            Err(_) => return,
        };
        <GtfAttributes as FeatureAttributes>::iter(&attrs).for_each(|result| {
            if let Ok((key, value)) = result {
                let attr_type = match value {
                    FeatureAttributeValue::String(_) => AttributeType::String,
                    FeatureAttributeValue::Array(_) => AttributeType::Array,
                };
                self.attrs
                    .entry(key.to_string())
                    .and_modify(|e| {
                        if attr_type == AttributeType::Array && *e == AttributeType::String {
                            *e = attr_type.clone();
                        }
                    })
                    .or_insert(attr_type);
            };
        });
    }
}

impl Push<noodles::gff::Record<'_>> for AttributeScanner {
    fn push(&mut self, record: noodles::gff::Record) {
        let attrs = record.attributes();
        attrs.iter().for_each(|result| {
            if let Ok((key, value)) = result {
                let attr_type = match value {
                    GffAttributeValue::String(_) => AttributeType::String,
                    GffAttributeValue::Array(_) => AttributeType::Array,
                };
                self.attrs
                    .entry(key.to_string())
                    .and_modify(|e| {
                        if attr_type == AttributeType::Array && *e == AttributeType::String {
                            *e = attr_type.clone();
                        }
                    })
                    .or_insert(attr_type);
            };
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{Array, StringArray};

    #[test]
    fn test_attribute_def_try_from_valid() {
        let def = AttributeDef::try_from(("name".to_string(), "String".to_string())).unwrap();
        assert_eq!(def.name, "name");
        assert_eq!(def.ty, AttributeType::String);

        let def = AttributeDef::try_from(("name".to_string(), "Array".to_string())).unwrap();
        assert_eq!(def.name, "name");
        assert_eq!(def.ty, AttributeType::Array);
    }

    #[test]
    fn test_attribute_def_try_from_invalid() {
        let result = AttributeDef::try_from(("name".to_string(), "InvalidType".to_string()));
        assert!(result.is_err());
    }

    #[test]
    fn test_attribute_arrow_type() {
        for ty in [AttributeType::String, AttributeType::Array] {
            let mut builder = AttributeBuilder::new(&ty);
            let data_type = builder.finish().data_type().clone();
            assert_eq!(ty.arrow_type(), data_type);
        }
    }

    #[test]
    fn test_attribute_builder_append_null() {
        let mut builder = AttributeBuilder::new(&AttributeType::String);
        builder.append_null();
        let array = builder.finish();
        assert!(array.is_nullable());
        assert!(array.is_null(0));
    }

    #[test]
    fn test_attribute_builder_append_value_string() {
        let mut builder = AttributeBuilder::new(&AttributeType::String);
        builder
            .append_value(&AttributeValue::String("value".to_string()))
            .unwrap();
        let array = builder.finish();
        let string_array = array.as_any().downcast_ref::<StringArray>().unwrap();
        assert_eq!(string_array.len(), 1);
        assert_eq!(string_array.value(0), "value");
    }

    #[test]
    fn test_attribute_builder_append_value_array() {
        let mut builder = AttributeBuilder::new(&AttributeType::Array);
        builder
            .append_value(&AttributeValue::Array(vec![
                "value1".to_string(),
                "value2".to_string(),
            ]))
            .unwrap();
        let array = builder.finish();
        let list_array = array
            .as_any()
            .downcast_ref::<arrow::array::ListArray>()
            .unwrap();
        let values = list_array.value(0);
        let item = values.as_any().downcast_ref::<StringArray>().unwrap();
        assert_eq!(item.value(0), "value1");
        assert_eq!(item.value(1), "value2");
    }
}

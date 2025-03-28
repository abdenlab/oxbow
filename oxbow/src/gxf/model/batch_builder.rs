use std::io;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::{ArrayRef, StructArray};
use arrow::datatypes::FieldRef;
use arrow::datatypes::{DataType, Field as ArrowField, Schema};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use indexmap::IndexMap;

use crate::gxf::model::attribute::{AttributeBuilder, AttributeDef, AttributeValue};
use crate::gxf::model::field::Push as _;
use crate::gxf::model::field::{Field, FieldBuilder, DEFAULT_FIELD_NAMES};

/// A builder for an Arrow record batch of GXF (GTF/GFF) features.
pub struct BatchBuilder {
    fields: Vec<Field>,
    attr_defs: Vec<AttributeDef>,
    field_builders: IndexMap<Field, FieldBuilder>,
    attr_builders: IndexMap<AttributeDef, AttributeBuilder>,
}

impl BatchBuilder {
    /// Creates a new `BatchBuilder` for GTF/GFF records.
    pub fn new(
        field_names: Option<Vec<String>>,
        attr_defs: Option<Vec<(String, String)>>,
        capacity: usize,
    ) -> io::Result<Self> {
        let default_field_names: Vec<String> = DEFAULT_FIELD_NAMES
            .into_iter()
            .map(|name| name.to_string())
            .collect();
        let fields: Vec<Field> = field_names
            .unwrap_or_else(|| default_field_names)
            .into_iter()
            .map(|name| Field::from_str(&name))
            .collect::<Result<Vec<_>, _>>()?;
        let mut field_builders = IndexMap::new();
        for field in &fields {
            let builder = FieldBuilder::new(field.clone(), capacity);
            field_builders.insert(field.clone(), builder);
        }

        let attr_defs: Vec<AttributeDef> = attr_defs
            .unwrap_or_else(|| vec![])
            .into_iter()
            .map(AttributeDef::try_from)
            .collect::<Result<Vec<_>, _>>()?;
        let mut attr_builders = IndexMap::new();
        for def in &attr_defs {
            let builder = AttributeBuilder::new(&def.ty);
            attr_builders.insert(def.clone(), builder);
        }

        Ok(Self {
            fields,
            attr_defs,
            field_builders,
            attr_builders,
        })
    }

    pub fn get_arrow_fields(&self) -> Vec<ArrowField> {
        // fixed fields
        let mut fields: Vec<ArrowField> = self
            .fields
            .iter()
            .map(|field| field.get_arrow_field())
            .collect();

        // attributes (optional)
        if !self.attr_defs.is_empty() {
            let nested_fields: Vec<ArrowField> = self
                .attr_defs
                .iter()
                .map(|def| def.get_arrow_field())
                .collect();
            let tag_field = ArrowField::new(
                "attributes",
                DataType::Struct(arrow::datatypes::Fields::from(nested_fields)),
                true,
            );
            fields.push(tag_field);
        }

        fields
    }

    pub fn get_arrow_schema(&self) -> Schema {
        Schema::new(self.get_arrow_fields())
    }

    pub fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        // fixed fields
        let mut name_to_array: Vec<(&str, ArrayRef)> = self
            .field_builders
            .iter_mut()
            .map(|(field, builder)| {
                let name = field.name();
                (name, builder.finish())
            })
            .collect();

        // attributes (optional)
        if !self.attr_defs.is_empty() {
            let attr_arrays: Vec<(FieldRef, ArrayRef)> = self
                .attr_builders
                .iter_mut()
                .map(|(def, builder)| {
                    let arrow_field = def.get_arrow_field();
                    (Arc::new(arrow_field), builder.finish())
                })
                .collect();
            let attrs = StructArray::from(attr_arrays);
            name_to_array.push(("attributes", Arc::new(attrs)));
        }

        RecordBatch::try_from_iter(name_to_array.into_iter())
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Append a GFF record to the batch.
impl<'a> Push<&'a noodles::gff::Record<'a>> for BatchBuilder {
    fn push(&mut self, record: &noodles::gff::Record) -> io::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }

        // attributes (optional)
        if !self.attr_defs.is_empty() {
            let attrs = record.attributes();
            for (def, builder) in self.attr_builders.iter_mut() {
                match attrs.get(&def.name) {
                    // attribute is present in this record
                    Some(result) => {
                        match result {
                            // attribute parsed successfully
                            Ok(value) => {
                                let value = AttributeValue::from(&value);
                                builder.append_value(&value)?;
                            }
                            // attribute could not be parsed
                            Err(_) => {
                                eprintln!("Error parsing tag: {:?}", def.name);
                                builder.append_null();
                            }
                        }
                    }
                    // attribute is not present in this record
                    None => {
                        builder.append_null();
                    }
                };
            }
        }
        Ok(())
    }
}

/// Append a GTF record to the batch.
impl Push<&noodles::gtf::Record> for BatchBuilder {
    fn push(&mut self, record: &noodles::gtf::Record) -> io::Result<()> {
        for (_, builder) in self.field_builders.iter_mut() {
            builder.push(record)?;
        }

        // attributes (optional)
        if !self.attr_defs.is_empty() {
            let attrs = record.attributes();
            for (def, builder) in self.attr_builders.iter_mut() {
                match attrs.get(&def.name) {
                    Some(value) => {
                        let value = AttributeValue::String(value.to_string());
                        builder.append_value(&value)?;
                    }
                    None => {
                        builder.append_null();
                    }
                };
            }
        }
        Ok(())
    }
}

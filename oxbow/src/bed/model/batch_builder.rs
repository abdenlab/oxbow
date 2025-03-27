use std::io;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::ArrayRef;
use arrow::array::GenericStringBuilder;
use arrow::datatypes::{DataType, Field as ArrowField, Schema};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use indexmap::IndexMap;

use super::field::Push as _;
use super::field::{Field, FieldBuilder, DEFAULT_FIELD_NAMES};
use super::schema::BedSchema;

/// A builder for an Arrow record batch of BED features.
pub struct BatchBuilder {
    standard_fields: Vec<Field>,
    custom_field_names: Vec<String>,
    standard_field_builders: IndexMap<Field, FieldBuilder>,
    custom_field_builders: IndexMap<String, GenericStringBuilder<i32>>,
}

impl BatchBuilder {
    pub fn new(field_names: Option<Vec<String>>, schema: &BedSchema, capacity: usize) -> Self {
        let n = schema.standard_field_count();
        let default_field_names: Vec<String> = DEFAULT_FIELD_NAMES
            .into_iter()
            .take(n)
            .map(|name| name.to_string())
            .collect();
        let standard_fields: Vec<Field> = field_names
            .unwrap_or_else(|| default_field_names)
            .into_iter()
            .map(|name| Field::from_str(&name))
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let mut standard_field_builders = IndexMap::new();
        for field in &standard_fields {
            let builder = FieldBuilder::new(field.clone(), capacity);
            standard_field_builders.insert(field.clone(), builder);
        }

        let mut custom_field_builders = IndexMap::new();
        let mut custom_field_names = Vec::new();
        match schema.custom_field_count() {
            Some(m) => {
                for i in 1..=m {
                    let name = format!("BED{}+{}", n, i);
                    let builder = GenericStringBuilder::<i32>::new();
                    custom_field_builders.insert(name.clone(), builder);
                    custom_field_names.push(name.clone());
                }
            }
            None => {
                panic!("Extended BED schemas are only supported with a known number of additional fields");
                // let builder = GenericStringBuilder::<i32>::new();
                // let name = "rest".to_string();
                // custom_field_builders.insert(name.clone(), builder);
                // custom_field_names.push(name.clone());
            }
        }

        Self {
            standard_fields,
            custom_field_names,
            standard_field_builders,
            custom_field_builders,
        }
    }

    pub fn get_arrow_fields(&self) -> Vec<ArrowField> {
        // standard fields
        let mut fields: Vec<ArrowField> = self
            .standard_fields
            .iter()
            .map(|field| field.get_arrow_field())
            .collect();

        // custom fields (optional)
        if !self.custom_field_builders.is_empty() {
            let other_fields: Vec<ArrowField> = self
                .custom_field_names
                .iter()
                .map(|name| ArrowField::new(name, DataType::Utf8, true))
                .collect();
            fields.extend(other_fields);
        }

        fields
    }

    pub fn get_arrow_schema(&self) -> Schema {
        Schema::new(self.get_arrow_fields())
    }

    pub fn finish(&mut self) -> Result<RecordBatch, ArrowError> {
        // standard fields
        let mut name_to_array: Vec<(&str, ArrayRef)> = self
            .standard_field_builders
            .iter_mut()
            .map(|(field, builder)| {
                let name = field.name();
                (name, builder.finish())
            })
            .collect();

        // custom fields (optional)
        if !self.custom_field_builders.is_empty() {
            let other: Vec<(&str, ArrayRef)> = self
                .custom_field_builders
                .iter_mut()
                .map(|(name, builder)| (name.as_str(), Arc::new(builder.finish()) as ArrayRef))
                .collect();
            name_to_array.extend(other);
        }

        RecordBatch::try_from_iter(name_to_array.into_iter())
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Appends a BED record to the batch.
impl Push<&noodles::bed::Record<3>> for BatchBuilder {
    fn push(&mut self, record: &noodles::bed::Record<3>) -> io::Result<()> {
        // standard fields
        for (_, builder) in self.standard_field_builders.iter_mut() {
            builder.push(record)?;
        }

        // custom fields
        let n = self.standard_fields.len();
        if !self.custom_field_builders.is_empty() {
            for (i, value) in record.other_fields().iter().enumerate() {
                let field_abs_idx = i + 3;
                let field_rel_idx = match field_abs_idx.checked_sub(n) {
                    Some(i) => i,
                    None => continue,
                };
                let name = format!("BED{}+{}", n, field_rel_idx + 1);
                let builder = match self.custom_field_builders.get_mut(&name) {
                    Some(builder) => builder,
                    None => continue, // skip if we are not collecting this field
                };
                builder.append_value(value.to_string());
            }
        }

        Ok(())
    }
}

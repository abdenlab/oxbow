pub mod batch;
pub mod field;
pub mod field_def;
pub mod schema;

pub use batch::BatchBuilder;
pub use field_def::{
    bed_standard_fields, FieldBuilder as GenericFieldBuilder, FieldDef, FieldType,
};
pub use schema::BedSchema;

use std::sync::Arc;

use arrow::datatypes::{Field as ArrowField, Schema, SchemaRef};

use crate::OxbowError;

/// A data model for BED/BBI records.
///
/// Wraps a [`BedSchema`] (which defines the parsing interpretation) with an
/// optional field projection to select which columns to include in output.
///
/// - `bed_schema` defines how to parse each record (standard + custom fields).
/// - `fields` projects which of those fields become Arrow columns.
///   `None` → all fields from the bed_schema.
///
/// # Examples
///
/// ```
/// use oxbow::bed::model::{Model, BedSchema};
///
/// // BED6 with all fields.
/// let bed_schema: BedSchema = "bed6".parse().unwrap();
/// let model = Model::new(bed_schema, None).unwrap();
/// assert_eq!(model.field_names().len(), 6);
///
/// // BED6 projected to 3 fields.
/// let bed_schema: BedSchema = "bed6".parse().unwrap();
/// let model = Model::new(bed_schema, Some(vec!["chrom".into(), "start".into(), "end".into()])).unwrap();
/// assert_eq!(model.field_names().len(), 3);
/// ```
#[derive(Clone, Debug)]
pub struct Model {
    bed_schema: BedSchema,
    fields: Vec<FieldDef>,
    schema: SchemaRef,
}

impl Model {
    /// Create a new BED model.
    ///
    /// - `bed_schema`: the parsing interpretation.
    /// - `fields`: column names to project. `None` → all fields from the schema.
    pub fn new(bed_schema: BedSchema, fields: Option<Vec<String>>) -> crate::Result<Self> {
        let projected = match fields {
            None => bed_schema.fields().clone(),
            Some(names) => {
                let available = bed_schema.fields();
                let mut projected = Vec::new();
                for name in &names {
                    let def = available
                        .iter()
                        .find(|d| d.name.eq_ignore_ascii_case(name))
                        .ok_or_else(|| {
                            OxbowError::invalid_input(format!(
                                "Field '{}' not in BED schema. Available: {:?}",
                                name,
                                bed_schema.field_names()
                            ))
                        })?;
                    projected.push(def.clone());
                }
                projected
            }
        };

        let arrow_fields: Vec<ArrowField> = projected.iter().map(|d| d.get_arrow_field()).collect();
        let schema = Arc::new(Schema::new(arrow_fields));

        Ok(Self {
            bed_schema,
            fields: projected,
            schema,
        })
    }

    /// The BED schema (parsing interpretation).
    pub fn bed_schema(&self) -> &BedSchema {
        &self.bed_schema
    }

    /// The projected field definitions.
    pub fn fields(&self) -> &[FieldDef] {
        &self.fields
    }

    /// The projected field names.
    pub fn field_names(&self) -> Vec<String> {
        self.fields.iter().map(|d| d.name.clone()).collect()
    }

    /// The Arrow schema for the projected fields.
    pub fn schema(&self) -> &SchemaRef {
        &self.schema
    }

    /// All column names in the Arrow schema.
    pub fn column_names(&self) -> Vec<String> {
        self.schema
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect()
    }

    /// Create a projected model containing only the specified columns.
    ///
    /// Returns an error if any column name is not in this model's schema.
    pub fn project(&self, columns: &[String]) -> crate::Result<Self> {
        let available = self.column_names();
        let unknown: Vec<&str> = columns
            .iter()
            .filter(|c| !available.iter().any(|a| a.eq_ignore_ascii_case(c)))
            .map(|c| c.as_str())
            .collect();
        if !unknown.is_empty() {
            return Err(OxbowError::invalid_input(format!(
                "Unknown columns: {:?}. Available: {:?}",
                unknown, available
            )));
        }

        let projected: Vec<String> = self
            .fields
            .iter()
            .filter(|d| columns.iter().any(|c| c.eq_ignore_ascii_case(&d.name)))
            .map(|d| d.name.clone())
            .collect();

        Self::new(self.bed_schema.clone(), Some(projected))
    }
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.bed_schema == other.bed_schema && self.fields == other.fields
    }
}

impl Eq for Model {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed6_all_fields() {
        let bed_schema: BedSchema = "bed6".parse().unwrap();
        let model = Model::new(bed_schema, None).unwrap();
        assert_eq!(model.field_names().len(), 6);
        assert_eq!(model.schema().fields().len(), 6);
    }

    #[test]
    fn test_bed6_projected() {
        let bed_schema: BedSchema = "bed6".parse().unwrap();
        let model = Model::new(
            bed_schema,
            Some(vec!["chrom".into(), "start".into(), "end".into()]),
        )
        .unwrap();
        assert_eq!(model.field_names(), vec!["chrom", "start", "end"]);
        assert_eq!(model.schema().fields().len(), 3);
    }

    #[test]
    fn test_bed3_plus() {
        let bed_schema: BedSchema = "bed3+".parse().unwrap();
        let model = Model::new(bed_schema, None).unwrap();
        assert_eq!(model.field_names(), vec!["chrom", "start", "end", "rest"]);
    }

    #[test]
    fn test_bedgraph() {
        let bed_schema = BedSchema::new_bedgraph().unwrap();
        let model = Model::new(bed_schema, None).unwrap();
        assert_eq!(model.field_names(), vec!["chrom", "start", "end", "value"]);
    }

    #[test]
    fn test_custom_schema() {
        let defs = vec![
            FieldDef::new("chrom".into(), FieldType::String),
            FieldDef::new("start".into(), FieldType::Uint),
            FieldDef::new("end".into(), FieldType::Uint),
            FieldDef::new("signalValue".into(), FieldType::Float),
        ];
        let bed_schema = BedSchema::from_defs(defs).unwrap();
        let model = Model::new(bed_schema, None).unwrap();
        assert_eq!(
            model.field_names(),
            vec!["chrom", "start", "end", "signalValue"]
        );
    }

    #[test]
    fn test_project() {
        let bed_schema: BedSchema = "bed6+3".parse().unwrap();
        let model = Model::new(bed_schema, None).unwrap();
        let projected = model
            .project(&["chrom".into(), "end".into(), "BED6+1".into()])
            .unwrap();
        assert_eq!(projected.field_names(), vec!["chrom", "end", "BED6+1"]);
    }

    #[test]
    fn test_project_unknown() {
        let bed_schema: BedSchema = "bed3".parse().unwrap();
        let model = Model::new(bed_schema, None).unwrap();
        let result = model.project(&["nonexistent".into()]);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_field_name() {
        let bed_schema: BedSchema = "bed3".parse().unwrap();
        let result = Model::new(bed_schema, Some(vec!["nonexistent".into()]));
        assert!(result.is_err());
    }
}

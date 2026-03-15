pub mod batch;
pub mod field;
pub mod field_def;
pub mod schema;

pub use batch::BatchBuilder;
pub use field_def::{
    bed_standard_fields, FieldBuilder as GenericFieldBuilder, FieldDef, FieldType,
};
pub use schema::BedSchema;

use std::str::FromStr;
use std::sync::Arc;

use arrow::datatypes::{Field as ArrowField, Schema, SchemaRef};

use crate::OxbowError;
use field::Field;

/// A data model for BED records.
///
/// Wraps a [`BedSchema`] (which defines the parsing interpretation) with an
/// optional field projection to select which columns to include in output.
///
/// Uses BED-specific Arrow types for standard fields (e.g., Int64 for
/// positions) and FieldDef types for custom fields.
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
    field_names: Vec<String>,
    schema: SchemaRef,
}

impl Model {
    /// Create a new BED model.
    ///
    /// - `bed_schema`: the parsing interpretation.
    /// - `fields`: column names to project. `None` → all fields from the schema.
    pub fn new(bed_schema: BedSchema, fields: Option<Vec<String>>) -> crate::Result<Self> {
        let available_names = bed_schema.field_names();
        let projected_names = match fields {
            None => available_names.clone(),
            Some(names) => {
                for name in &names {
                    if !available_names.iter().any(|a| a.eq_ignore_ascii_case(name)) {
                        return Err(OxbowError::invalid_input(format!(
                            "Field '{}' not in BED schema. Available: {:?}",
                            name, available_names
                        )));
                    }
                }
                names
            }
        };

        let standard_names = bed_schema.standard_field_names();
        let custom_fields = bed_schema.custom_fields();
        let arrow_fields: Vec<ArrowField> = projected_names
            .iter()
            .map(|name| {
                // Standard fields: use the specialized BED Field types
                if standard_names.iter().any(|s| s.eq_ignore_ascii_case(name)) {
                    if let Ok(field) = Field::from_str(name) {
                        return field.get_arrow_field();
                    }
                }
                // Custom fields: use FieldDef types
                if let Some(def) = custom_fields
                    .iter()
                    .find(|d| d.name.eq_ignore_ascii_case(name))
                {
                    return def.get_arrow_field();
                }
                // Fallback (shouldn't happen after validation)
                ArrowField::new(name, arrow::datatypes::DataType::Utf8, true)
            })
            .collect();
        let schema = Arc::new(Schema::new(arrow_fields));

        Ok(Self {
            bed_schema,
            field_names: projected_names,
            schema,
        })
    }

    /// The BED schema (parsing interpretation).
    pub fn bed_schema(&self) -> &BedSchema {
        &self.bed_schema
    }

    /// The projected field names.
    pub fn field_names(&self) -> Vec<String> {
        self.field_names.clone()
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
            .field_names
            .iter()
            .filter(|n| columns.iter().any(|c| c.eq_ignore_ascii_case(n)))
            .cloned()
            .collect();

        Self::new(self.bed_schema.clone(), Some(projected))
    }
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.bed_schema == other.bed_schema && self.field_names == other.field_names
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
            FieldDef::new("signalValue".into(), FieldType::Float),
            FieldDef::new("pValue".into(), FieldType::Float),
        ];
        let bed_schema = BedSchema::new(3, Some(defs)).unwrap();
        let model = Model::new(bed_schema, None).unwrap();
        assert_eq!(
            model.field_names(),
            vec!["chrom", "start", "end", "signalValue", "pValue"]
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

    #[test]
    fn test_bed3_projected_subset() {
        let bed_schema: BedSchema = "bed3".parse().unwrap();
        let model = Model::new(bed_schema, Some(vec!["chrom".into(), "end".into()])).unwrap();
        assert_eq!(model.field_names(), vec!["chrom", "end"]);
    }

    #[test]
    fn test_bed9_noncontiguous_projection() {
        let bed_schema: BedSchema = "bed9".parse().unwrap();
        let model = Model::new(
            bed_schema,
            Some(vec!["chrom".into(), "strand".into(), "itemRgb".into()]),
        )
        .unwrap();
        assert_eq!(model.field_names(), vec!["chrom", "strand", "itemRgb"]);
        assert_eq!(model.schema().fields().len(), 3);
    }

    #[test]
    fn test_bed12_with_custom_mixed_projection() {
        let defs = vec![
            FieldDef::new("extra1".into(), FieldType::Float),
            FieldDef::new("extra2".into(), FieldType::String),
        ];
        let bed_schema = BedSchema::new(12, Some(defs)).unwrap();
        let model = Model::new(bed_schema, None).unwrap();
        assert_eq!(model.field_names().len(), 14);

        let projected = model
            .project(&["chrom".into(), "blockSizes".into(), "extra1".into()])
            .unwrap();
        assert_eq!(
            projected.field_names(),
            vec!["chrom", "blockSizes", "extra1"]
        );
    }

    #[test]
    fn test_bedgraph_arrow_types() {
        use arrow::datatypes::DataType;

        let bed_schema = BedSchema::new_bedgraph().unwrap();
        let model = Model::new(bed_schema, None).unwrap();
        // Standard fields use BED types (Int64 for positions)
        assert_eq!(model.schema().field(1).data_type(), &DataType::Int64);
        // Custom "value" field uses FieldDef type (Float32)
        assert_eq!(model.schema().field(3).data_type(), &DataType::Float32);
    }
}

pub mod batch;
pub mod field;

pub use batch::BatchBuilder;

use std::sync::Arc;

use arrow::datatypes::{Field as ArrowField, Schema, SchemaRef};

use crate::{OxbowError, Select};
use field::{Field, FASTA_DEFAULT_FIELD_NAMES, FASTQ_DEFAULT_FIELD_NAMES};

/// A data model for sequence records (FASTA/FASTQ).
///
/// Encapsulates the schema-defining parameters for a sequence projection:
/// which fields to include.
///
/// - `fields` selects which fields become Arrow columns.
///   `All` → format-specific defaults (3 for FASTA, 4 for FASTQ).
///   `Omit` → no fields. `Some(vec)` → specific fields.
///
/// # Examples
///
/// ```
/// use oxbow::sequence::model::Model;
/// use oxbow::Select;
///
/// // FASTA defaults: name, description, sequence.
/// let model = Model::new_fasta(Select::All).unwrap();
/// assert_eq!(model.field_names().len(), 3);
///
/// // FASTQ defaults: name, description, sequence, quality.
/// let model = Model::new_fastq(Select::All).unwrap();
/// assert_eq!(model.field_names().len(), 4);
///
/// // Custom field selection.
/// let model = Model::new_fastq(Select::Some(vec!["name".into(), "sequence".into()])).unwrap();
/// assert_eq!(model.field_names(), vec!["name", "sequence"]);
/// ```
#[derive(Clone, Debug)]
pub struct Model {
    fields: Vec<Field>,
    schema: SchemaRef,
}

impl Model {
    /// Create a new FASTA model.
    ///
    /// `fields`: `All` → `["name", "description", "sequence"]`. `Omit` → no
    /// fields. `Some(vec)` → specific fields.
    pub fn new_fasta(fields: Select<String>) -> crate::Result<Self> {
        let defaults = || {
            FASTA_DEFAULT_FIELD_NAMES
                .iter()
                .map(|&s| s.to_string())
                .collect()
        };
        let field_names = match fields {
            Select::All => defaults(),
            Select::Some(names) => names,
            Select::Omit => Vec::new(),
        };
        Self::new(field_names)
    }

    /// Create a new FASTQ model.
    ///
    /// `fields`: `All` → `["name", "description", "sequence", "quality"]`.
    /// `Omit` → no fields. `Some(vec)` → specific fields.
    pub fn new_fastq(fields: Select<String>) -> crate::Result<Self> {
        let defaults = || {
            FASTQ_DEFAULT_FIELD_NAMES
                .iter()
                .map(|&s| s.to_string())
                .collect()
        };
        let field_names = match fields {
            Select::All => defaults(),
            Select::Some(names) => names,
            Select::Omit => Vec::new(),
        };
        Self::new(field_names)
    }

    fn new(field_names: Vec<String>) -> crate::Result<Self> {
        let mut parsed_fields = Vec::new();
        for name in &field_names {
            let field: Field = name
                .parse()
                .map_err(|e: std::io::Error| OxbowError::invalid_input(e.to_string()))?;
            parsed_fields.push(field);
        }

        let arrow_fields: Vec<ArrowField> =
            parsed_fields.iter().map(|f| f.get_arrow_field()).collect();
        let schema = Arc::new(Schema::new(arrow_fields));

        Ok(Self {
            fields: parsed_fields,
            schema,
        })
    }

    /// The validated fields.
    pub fn fields(&self) -> &[Field] {
        &self.fields
    }

    /// The field names.
    pub fn field_names(&self) -> Vec<String> {
        self.fields.iter().map(|f| f.name().to_string()).collect()
    }

    /// The Arrow schema for this model.
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
            .filter(|f| columns.iter().any(|c| c.eq_ignore_ascii_case(f.name())))
            .map(|f| f.name().to_string())
            .collect();

        Self::new(projected)
    }
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.fields == other.fields
    }
}

impl Eq for Model {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_defaults() {
        let model = Model::new_fasta(Select::All).unwrap();
        assert_eq!(model.field_names(), vec!["name", "description", "sequence"]);
        assert_eq!(model.schema().fields().len(), 3);
    }

    #[test]
    fn test_fastq_defaults() {
        let model = Model::new_fastq(Select::All).unwrap();
        assert_eq!(
            model.field_names(),
            vec!["name", "description", "sequence", "quality"]
        );
        assert_eq!(model.schema().fields().len(), 4);
    }

    #[test]
    fn test_custom_fields() {
        let model = Model::new_fastq(Select::Some(vec!["name".into(), "sequence".into()])).unwrap();
        assert_eq!(model.field_names(), vec!["name", "sequence"]);
        assert_eq!(model.schema().fields().len(), 2);
    }

    #[test]
    fn test_invalid_field() {
        let result = Model::new_fasta(Select::Some(vec!["invalid".into()]));
        assert!(result.is_err());
    }

    #[test]
    fn test_project() {
        let model = Model::new_fastq(Select::All).unwrap();
        let projected = model.project(&["name".into(), "quality".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["name", "quality"]);
    }

    #[test]
    fn test_project_unknown() {
        let model = Model::new_fasta(Select::All).unwrap();
        let result = model.project(&["nonexistent".into()]);
        assert!(result.is_err());
    }
}

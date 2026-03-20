pub mod batch;
pub mod field;

pub use batch::BatchBuilder;

use std::sync::Arc;

use arrow::datatypes::{Field as ArrowField, Schema, SchemaRef};

use crate::{CoordSystem, OxbowError, Select};
use field::{Field, DEFAULT_FIELD_NAMES};

pub struct BBIZoomRecord<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
    pub bases_covered: u64,
    pub min: f64,
    pub max: f64,
    pub sum: f64,
    pub sum_squares: f64,
}

impl<'a> BBIZoomRecord<'a> {
    pub fn new(chrom: &'a str, entry: &'a bigtools::ZoomRecord) -> Self {
        Self {
            chrom,
            start: entry.start,
            end: entry.end,
            bases_covered: entry.summary.bases_covered,
            min: entry.summary.min_val,
            max: entry.summary.max_val,
            sum: entry.summary.sum,
            sum_squares: entry.summary.sum_squares,
        }
    }
}

/// A data model for BBI zoom level summary records.
///
/// Fixed schema: 8 fields (chrom, start, end, bases_covered, min, max,
/// sum, sum_squares). The `fields` parameter projects which to include.
///
/// # Examples
///
/// ```
/// use oxbow::bbi::model::zoom::Model;
/// use oxbow::{CoordSystem, Select};
///
/// let model = Model::new(Select::All, CoordSystem::ZeroHalfOpen).unwrap();
/// assert_eq!(model.field_names().len(), 8);
///
/// let model = Model::new(Select::Some(vec!["chrom".into(), "start".into(), "end".into(), "sum".into()]), CoordSystem::ZeroHalfOpen).unwrap();
/// assert_eq!(model.field_names().len(), 4);
/// ```
#[derive(Clone, Debug)]
pub struct Model {
    fields: Vec<Field>,
    coord_system: CoordSystem,
    schema: SchemaRef,
}

impl Model {
    /// Create a new BBI zoom model.
    ///
    /// `fields`: field names. `None` → all 8 default fields.
    pub fn new(fields: Select<String>, coord_system: CoordSystem) -> crate::Result<Self> {
        let field_names = match fields {
            Select::All => DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect(),
            Select::Some(names) => names,
            Select::Omit => Vec::new(),
        };

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
            coord_system,
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

    /// The output coordinate system.
    pub fn coord_system(&self) -> CoordSystem {
        self.coord_system
    }

    /// The Arrow schema.
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
            .fields
            .iter()
            .filter(|f| columns.iter().any(|c| c.eq_ignore_ascii_case(f.name())))
            .map(|f| f.name().to_string())
            .collect();

        Self::new(Select::Some(projected), self.coord_system)
    }
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.fields == other.fields && self.coord_system == other.coord_system
    }
}

impl Eq for Model {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_defaults() {
        let model = Model::new(Select::All, CoordSystem::ZeroHalfOpen).unwrap();
        assert_eq!(model.field_names().len(), 8);
        assert_eq!(model.schema().fields().len(), 8);
    }

    #[test]
    fn test_custom() {
        let model = Model::new(
            Select::Some(vec!["chrom".into(), "start".into(), "sum".into()]),
            CoordSystem::ZeroHalfOpen,
        )
        .unwrap();
        assert_eq!(model.field_names(), vec!["chrom", "start", "sum"]);
    }

    #[test]
    fn test_project() {
        let model = Model::new(Select::All, CoordSystem::ZeroHalfOpen).unwrap();
        let projected = model
            .project(&["chrom".into(), "min".into(), "max".into()])
            .unwrap();
        assert_eq!(projected.field_names(), vec!["chrom", "min", "max"]);
    }

    #[test]
    fn test_invalid_field() {
        let result = Model::new(
            Select::Some(vec!["invalid".into()]),
            CoordSystem::ZeroHalfOpen,
        );
        assert!(result.is_err());
    }
}

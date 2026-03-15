pub mod attribute;
pub mod batch;
pub mod field;

pub use batch::BatchBuilder;

use std::sync::Arc;

use arrow::datatypes::{DataType, Field as ArrowField, Schema, SchemaRef};

use crate::OxbowError;
use attribute::AttributeDef;
use field::{Field, DEFAULT_FIELD_NAMES};

/// A data model for GXF (GFF/GTF) feature records.
///
/// Encapsulates the schema-defining parameters for a feature projection:
/// which standard fields to include and which attributes (with their types)
/// to materialize.
///
/// - `fields` selects which standard GXF fields become Arrow columns.
///   `None` → all 8 standard fields.
/// - `attr_defs` controls the attributes struct column independently.
///   `None` → no attributes column. `Some(vec![])` → empty struct column.
///   `Some(vec![...])` → struct column with the specified sub-fields.
///
/// The model can produce an Arrow schema independently of any file content.
///
/// # Examples
///
/// ```
/// use oxbow::gxf::model::Model;
///
/// // Default: all 8 standard fields, no attributes column.
/// let model = Model::new(None, None).unwrap();
/// assert_eq!(model.field_names().len(), 8);
/// assert!(!model.has_attributes());
///
/// // Custom: selected fields with attributes.
/// let model = Model::new(
///     Some(vec!["seqid".into(), "start".into(), "end".into()]),
///     Some(vec![("gene_id".into(), "String".into())]),
/// ).unwrap();
/// assert_eq!(model.field_names(), vec!["seqid", "start", "end"]);
/// assert!(model.has_attributes());
/// assert_eq!(model.attr_defs().unwrap().len(), 1);
/// // Schema: 4 columns (seqid, start, end, attributes{gene_id})
/// assert_eq!(model.schema().fields().len(), 4);
/// ```
#[derive(Clone, Debug)]
pub struct Model {
    fields: Vec<Field>,
    attr_defs: Option<Vec<AttributeDef>>,
    schema: SchemaRef,
}

impl Model {
    /// Create a new GXF model.
    ///
    /// - `fields`: standard GXF field names. `None` → all 8 standard fields.
    /// - `attr_defs`: attribute definitions as `(name, type)` pairs. `None` →
    ///   no attributes column. `Some(vec![])` → attributes column with empty
    ///   struct.
    pub fn new(
        fields: Option<Vec<String>>,
        attr_defs: Option<Vec<(String, String)>>,
    ) -> crate::Result<Self> {
        let field_names =
            fields.unwrap_or_else(|| DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect());

        let mut parsed_fields = Vec::new();
        for name in &field_names {
            let field: Field = name
                .parse()
                .map_err(|e: std::io::Error| OxbowError::invalid_input(e.to_string()))?;
            parsed_fields.push(field);
        }

        let attr_defs: Option<Vec<AttributeDef>> = attr_defs
            .map(|defs| {
                defs.into_iter()
                    .map(AttributeDef::try_from)
                    .collect::<crate::Result<Vec<_>>>()
            })
            .transpose()?;

        let schema = Self::build_schema(&parsed_fields, attr_defs.as_deref());
        Ok(Self {
            fields: parsed_fields,
            attr_defs,
            schema,
        })
    }

    /// Create a model with all 8 default standard fields and no attributes.
    pub fn default_fields() -> Self {
        Self::new(None, None).expect("default fields are always valid")
    }

    fn build_schema(fields: &[Field], attr_defs: Option<&[AttributeDef]>) -> SchemaRef {
        let mut arrow_fields: Vec<ArrowField> =
            fields.iter().map(|f| f.get_arrow_field()).collect();

        if let Some(defs) = attr_defs {
            let nested: Vec<ArrowField> = defs.iter().map(|d| d.get_arrow_field()).collect();
            arrow_fields.push(ArrowField::new(
                "attributes",
                DataType::Struct(arrow::datatypes::Fields::from(nested)),
                true,
            ));
        }

        Arc::new(Schema::new(arrow_fields))
    }

    /// The validated standard fields.
    pub fn fields(&self) -> &[Field] {
        &self.fields
    }

    /// The standard field names.
    pub fn field_names(&self) -> Vec<String> {
        self.fields.iter().map(|f| f.name().to_string()).collect()
    }

    /// The validated attribute definitions, if attributes are included.
    pub fn attr_defs(&self) -> Option<&[AttributeDef]> {
        self.attr_defs.as_deref()
    }

    /// Whether the attributes struct column is included.
    pub fn has_attributes(&self) -> bool {
        self.attr_defs.is_some()
    }

    /// The Arrow schema for this model.
    pub fn schema(&self) -> &SchemaRef {
        &self.schema
    }

    /// All top-level column names in the Arrow schema.
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

        let projected_fields: Vec<String> = self
            .fields
            .iter()
            .filter(|f| columns.iter().any(|c| c.eq_ignore_ascii_case(f.name())))
            .map(|f| f.name().to_string())
            .collect();

        let include_attrs =
            self.has_attributes() && columns.iter().any(|c| c.eq_ignore_ascii_case("attributes"));

        let attr_defs = if include_attrs {
            self.attr_defs
                .as_ref()
                .map(|defs| defs.iter().map(|d| d.to_tuple()).collect())
        } else {
            None
        };

        Self::new(Some(projected_fields), attr_defs)
    }
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.fields == other.fields && self.attr_defs == other.attr_defs
    }
}

impl Eq for Model {}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::datatypes::DataType;

    #[test]
    fn test_default_model() {
        let model = Model::new(None, None).unwrap();
        assert_eq!(model.field_names().len(), 8);
        assert!(!model.has_attributes());
        assert!(model.attr_defs().is_none());
        assert_eq!(model.schema().fields().len(), 8);
    }

    #[test]
    fn test_default_fields_constructor() {
        let model = Model::default_fields();
        assert_eq!(model, Model::new(None, None).unwrap());
    }

    #[test]
    fn test_custom_fields_no_attrs() {
        let model = Model::new(
            Some(vec!["seqid".into(), "start".into(), "end".into()]),
            None,
        )
        .unwrap();
        assert_eq!(model.field_names(), vec!["seqid", "start", "end"]);
        assert!(!model.has_attributes());
        assert_eq!(model.schema().fields().len(), 3);
    }

    #[test]
    fn test_fields_with_attrs() {
        let model = Model::new(
            Some(vec!["seqid".into(), "start".into()]),
            Some(vec![("gene_id".into(), "String".into())]),
        )
        .unwrap();
        assert_eq!(model.field_names(), vec!["seqid", "start"]);
        assert!(model.has_attributes());
        assert_eq!(model.attr_defs().unwrap().len(), 1);
        assert_eq!(model.schema().fields().len(), 3);
        assert_eq!(model.column_names(), vec!["seqid", "start", "attributes"]);
    }

    #[test]
    fn test_attrs_empty_defs_is_empty_struct() {
        let model = Model::new(Some(vec!["seqid".into()]), Some(vec![])).unwrap();
        assert!(model.has_attributes());
        assert!(model.attr_defs().unwrap().is_empty());
        assert_eq!(model.schema().fields().len(), 2);
        let attr_field = model.schema().field_with_name("attributes").unwrap();
        match attr_field.data_type() {
            DataType::Struct(fields) => assert!(fields.is_empty()),
            other => panic!("Expected Struct, got {:?}", other),
        }
    }

    #[test]
    fn test_no_attrs_when_attr_defs_none() {
        let model = Model::new(Some(vec!["seqid".into(), "start".into()]), None).unwrap();
        assert!(!model.has_attributes());
        assert!(model.attr_defs().is_none());
        assert_eq!(model.schema().fields().len(), 2);
    }

    #[test]
    fn test_invalid_field() {
        let result = Model::new(Some(vec!["invalid".into()]), None);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_attr_type() {
        let result = Model::new(None, Some(vec![("gene_id".into(), "InvalidType".into())]));
        assert!(result.is_err());
    }

    #[test]
    fn test_attr_defs_to_tuple() {
        let model = Model::new(
            None,
            Some(vec![
                ("gene_id".into(), "String".into()),
                ("tag".into(), "Array".into()),
            ]),
        )
        .unwrap();
        let tuples: Vec<_> = model
            .attr_defs()
            .unwrap()
            .iter()
            .map(|d| d.to_tuple())
            .collect();
        assert_eq!(
            tuples,
            vec![
                ("gene_id".to_string(), "String".to_string()),
                ("tag".to_string(), "Array".to_string()),
            ]
        );
    }

    #[test]
    fn test_project() {
        let model = Model::new(
            Some(vec!["seqid".into(), "start".into(), "end".into()]),
            Some(vec![("gene_id".into(), "String".into())]),
        )
        .unwrap();

        let projected = model.project(&["seqid".into(), "end".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["seqid", "end"]);
        assert!(!projected.has_attributes());
    }

    #[test]
    fn test_project_with_attrs() {
        let model = Model::new(
            Some(vec!["seqid".into(), "start".into()]),
            Some(vec![("gene_id".into(), "String".into())]),
        )
        .unwrap();

        let projected = model
            .project(&["seqid".into(), "attributes".into()])
            .unwrap();
        assert_eq!(projected.field_names(), vec!["seqid"]);
        assert!(projected.has_attributes());
        assert_eq!(projected.attr_defs().unwrap().len(), 1);
    }

    #[test]
    fn test_project_unknown_column() {
        let model = Model::default_fields();
        let result = model.project(&["nonexistent".into()]);
        assert!(result.is_err());
    }
}

pub mod batch;
pub mod field;
pub mod tag;

pub use batch::BatchBuilder;

use std::fmt;
use std::str::FromStr;
use std::sync::Arc;

use arrow::datatypes::{DataType, Field as ArrowField, Schema, SchemaRef};

use crate::{CoordSystem, OxbowError, Select};
use field::{Field, DEFAULT_FIELD_NAMES};
use tag::TagDef;

/// A data model for alignment records (SAM/BAM/CRAM).
///
/// Encapsulates the schema-defining parameters for an alignment projection:
/// which standard fields to include, which auxiliary tags (with their
/// types) to materialize, and the output coordinate system.
///
/// - `fields` selects which standard SAM fields become Arrow columns.
///   `None` → all 12 standard fields.
/// - `tag_defs` controls the tags struct column independently.
///   `None` → no tags column. `Some(vec![])` → empty struct column.
///   `Some(vec![...])` → struct column with the specified sub-fields.
/// - `coord_system` controls the coordinate system of position columns
///   (`pos`, `pnext`). `None` → 1-based closed (`"11"`), the SAM convention.
///   End coordinates are not affected.
///
/// The model can produce an Arrow schema independently of any file header.
///
/// # Examples
///
/// ```
/// use oxbow::alignment::model::Model;
/// use oxbow::{CoordSystem, Select};
///
/// // Default: all 12 standard fields, no tags column, 1-based coordinates.
/// let model = Model::new(Select::All, None, CoordSystem::OneClosed).unwrap();
/// assert_eq!(model.field_names().len(), 12);
/// assert!(!model.has_tags());
/// assert_eq!(model.coord_system(), CoordSystem::OneClosed);
///
/// // Custom: selected fields with tags, 0-based coordinates.
/// let model = Model::new(
///     Select::Some(vec!["qname".into(), "pos".into()]),
///     Some(vec![("NM".into(), "i".into()), ("MD".into(), "Z".into())]),
///     CoordSystem::ZeroHalfOpen,
/// ).unwrap();
/// assert_eq!(model.field_names(), vec!["qname", "pos"]);
/// assert!(model.has_tags());
/// assert_eq!(model.tag_defs().unwrap().len(), 2);
/// // Schema: 3 columns (qname, pos, tags{NM, MD})
/// assert_eq!(model.schema().fields().len(), 3);
/// ```
#[derive(Clone, Debug)]
pub struct Model {
    fields: Vec<Field>,
    tag_defs: Option<Vec<TagDef>>,
    coord_system: CoordSystem,
    schema: SchemaRef,
}

impl Model {
    /// Create a new alignment model.
    ///
    /// - `fields`: standard SAM field selection. `All` → all 12 standard
    ///   fields. `Select(vec)` → specific fields. `Omit` → no fields.
    /// - `tag_defs`: tag definitions as `(name, type_code)` pairs. `None` →
    ///   no tags column. `Some(vec![])` → tags column with empty struct.
    /// - `coord_system`: output coordinate system for position columns
    ///   (`pos`, `pnext`). `None` defaults to [`CoordSystem::OneClosed`]
    ///   (1-based, matching the SAM convention).
    pub fn new(
        fields: Select<String>,
        tag_defs: Option<Vec<(String, String)>>,
        coord_system: CoordSystem,
    ) -> crate::Result<Self> {
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

        let tag_defs: Option<Vec<TagDef>> = tag_defs
            .map(|defs| {
                defs.into_iter()
                    .map(TagDef::try_from)
                    .collect::<crate::Result<Vec<_>>>()
            })
            .transpose()?;
        let schema = Self::build_schema(&parsed_fields, tag_defs.as_deref());
        Ok(Self {
            fields: parsed_fields,
            tag_defs,
            coord_system,
            schema,
        })
    }

    fn build_schema(fields: &[Field], tag_defs: Option<&[TagDef]>) -> SchemaRef {
        let mut arrow_fields: Vec<ArrowField> =
            fields.iter().map(|f| f.get_arrow_field()).collect();

        if let Some(defs) = tag_defs {
            let nested: Vec<ArrowField> = defs.iter().map(|d| d.get_arrow_field()).collect();
            arrow_fields.push(ArrowField::new(
                "tags",
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

    /// The standard field names (lowercase).
    pub fn field_names(&self) -> Vec<String> {
        self.fields.iter().map(|f| f.to_string()).collect()
    }

    /// The validated tag definitions, if tags are included.
    pub fn tag_defs(&self) -> Option<&[TagDef]> {
        self.tag_defs.as_deref()
    }

    /// Whether the tags struct column is included.
    pub fn has_tags(&self) -> bool {
        self.tag_defs.is_some()
    }

    /// The output coordinate system for position columns.
    pub fn coord_system(&self) -> CoordSystem {
        self.coord_system
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
            .filter(|f| {
                columns
                    .iter()
                    .any(|c| c.eq_ignore_ascii_case(&f.to_string()))
            })
            .map(|f| f.to_string())
            .collect();

        let include_tags =
            self.has_tags() && columns.iter().any(|c| c.eq_ignore_ascii_case("tags"));

        let tag_defs = if include_tags {
            self.tag_defs
                .as_ref()
                .map(|defs| defs.iter().map(|d| d.to_tuple()).collect())
        } else {
            None
        };

        Self::new(Select::Some(projected_fields), tag_defs, self.coord_system)
    }
}

impl Default for Model {
    fn default() -> Self {
        Self::new(Select::All, None, CoordSystem::OneClosed)
            .expect("default fields are always valid")
    }
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.fields == other.fields
            && self.tag_defs == other.tag_defs
            && self.coord_system == other.coord_system
    }
}

impl Eq for Model {}

impl fmt::Display for Model {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let default_names: Vec<String> =
            DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()).collect();
        let current_names = self.field_names();

        if current_names == default_names {
            write!(f, "fields=*")?;
        } else {
            write!(f, "fields={}", current_names.join(","))?;
        }

        if let Some(defs) = &self.tag_defs {
            if defs.is_empty() {
                write!(f, ";tags")?;
            } else {
                let tags: Vec<String> = defs.iter().map(|d| d.to_string()).collect();
                write!(f, ";tags={}", tags.join(","))?;
            }
        }

        if self.coord_system != CoordSystem::OneClosed {
            write!(f, ";coords={}", self.coord_system)?;
        }

        Ok(())
    }
}

impl FromStr for Model {
    type Err = OxbowError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields: Option<Vec<String>> = None;
        let mut tag_defs: Option<Vec<(String, String)>> = None;
        let mut coord_system: Option<CoordSystem> = None;

        for part in s.split(';') {
            let part = part.trim();
            if let Some(value) = part.strip_prefix("fields=") {
                let mut names: Vec<String> = Vec::new();
                for name in value.split(',') {
                    let name = name.trim();
                    if name == "*" {
                        names.extend(DEFAULT_FIELD_NAMES.iter().map(|&s| s.to_string()));
                    } else if !name.is_empty() {
                        names.push(name.to_string());
                    }
                }
                fields = Some(names);
            } else if part == "tags" {
                // tags with no definitions
                tag_defs = Some(Vec::new());
            } else if let Some(value) = part.strip_prefix("tags=") {
                let defs: crate::Result<Vec<(String, String)>> = value
                    .split(',')
                    .map(|s| {
                        let s = s.trim();
                        let colon = s.find(':').ok_or_else(|| {
                            OxbowError::invalid_input(format!(
                                "Invalid tag def '{}': expected NAME:TYPE",
                                s
                            ))
                        })?;
                        Ok((s[..colon].to_string(), s[colon + 1..].to_string()))
                    })
                    .collect();
                tag_defs = Some(defs?);
            } else if let Some(value) = part.strip_prefix("coords=") {
                coord_system = Some(value.parse()?);
            } else {
                return Err(OxbowError::invalid_input(format!(
                    "Invalid Model segment: '{}'",
                    part
                )));
            }
        }

        let fields = match fields {
            Some(names) => Select::Some(names),
            None => Select::All,
        };
        Self::new(
            fields,
            tag_defs,
            coord_system.unwrap_or(CoordSystem::OneClosed),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const CS: CoordSystem = CoordSystem::OneClosed;

    #[test]
    fn test_default_model() {
        let model = Model::new(Select::All, None, CS).unwrap();
        assert_eq!(model.field_names().len(), 12);
        assert!(!model.has_tags());
        assert!(model.tag_defs().is_none());
        assert_eq!(model.schema().fields().len(), 12);
        assert_eq!(model.coord_system(), CoordSystem::OneClosed);
    }

    #[test]
    fn test_default_fields_constructor() {
        let model = Model::default();
        assert_eq!(model, Model::new(Select::All, None, CS).unwrap());
    }

    #[test]
    fn test_custom_fields_no_tags() {
        let model = Model::new(
            Select::Some(vec!["qname".into(), "flag".into(), "pos".into()]),
            None,
            CS,
        )
        .unwrap();
        assert_eq!(model.field_names(), vec!["qname", "flag", "pos"]);
        assert!(!model.has_tags());
        assert_eq!(model.schema().fields().len(), 3);
    }

    #[test]
    fn test_fields_with_tags() {
        let model = Model::new(
            Select::Some(vec!["qname".into(), "pos".into()]),
            Some(vec![("NM".into(), "i".into()), ("MD".into(), "Z".into())]),
            CS,
        )
        .unwrap();
        assert_eq!(model.field_names(), vec!["qname", "pos"]);
        assert!(model.has_tags());
        assert_eq!(model.tag_defs().unwrap().len(), 2);
        // 2 standard fields + 1 tags struct
        assert_eq!(model.schema().fields().len(), 3);
        assert_eq!(model.column_names(), vec!["qname", "pos", "tags"]);
    }

    #[test]
    fn test_tags_empty_defs_is_empty_struct() {
        let model = Model::new(Select::Some(vec!["qname".into()]), Some(vec![]), CS).unwrap();
        assert!(model.has_tags());
        assert!(model.tag_defs().unwrap().is_empty());
        assert_eq!(model.schema().fields().len(), 2);
        let tags_field = model.schema().field_with_name("tags").unwrap();
        match tags_field.data_type() {
            DataType::Struct(fields) => assert!(fields.is_empty()),
            other => panic!("Expected Struct, got {:?}", other),
        }
    }

    #[test]
    fn test_no_tags_when_tag_defs_none() {
        let model = Model::new(Select::Some(vec!["qname".into(), "pos".into()]), None, CS).unwrap();
        assert!(!model.has_tags());
        assert!(model.tag_defs().is_none());
        assert_eq!(model.schema().fields().len(), 2);
    }

    #[test]
    fn test_invalid_field() {
        let result = Model::new(Select::Some(vec!["invalid".into()]), None, CS);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_tag_name() {
        let result = Model::new(Select::All, Some(vec![("X".into(), "i".into())]), CS);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_tag_type() {
        let result = Model::new(Select::All, Some(vec![("NM".into(), "Q".into())]), CS);
        assert!(result.is_err());
    }

    #[test]
    fn test_project() {
        let model = Model::new(
            Select::Some(vec!["qname".into(), "flag".into(), "pos".into()]),
            Some(vec![("NM".into(), "i".into())]),
            CS,
        )
        .unwrap();

        let projected = model.project(&["qname".into(), "pos".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["qname", "pos"]);
        assert!(!projected.has_tags());
    }

    #[test]
    fn test_project_with_tags() {
        let model = Model::new(
            Select::Some(vec!["qname".into(), "pos".into()]),
            Some(vec![("NM".into(), "i".into())]),
            CS,
        )
        .unwrap();

        let projected = model.project(&["qname".into(), "tags".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["qname"]);
        assert!(projected.has_tags());
        assert_eq!(projected.tag_defs().unwrap().len(), 1);
    }

    #[test]
    fn test_project_unknown_column() {
        let model = Model::default();
        let result = model.project(&["nonexistent".into()]);
        assert!(result.is_err());
    }

    #[test]
    fn test_project_propagates_coord_system() {
        let model = Model::new(
            Select::Some(vec!["qname".into(), "pos".into()]),
            None,
            CoordSystem::ZeroHalfOpen,
        )
        .unwrap();
        let projected = model.project(&["pos".into()]).unwrap();
        assert_eq!(projected.coord_system(), CoordSystem::ZeroHalfOpen);
    }

    #[test]
    fn test_display_defaults() {
        let model = Model::default();
        assert_eq!(model.to_string(), "fields=*");
    }

    #[test]
    fn test_display_custom_with_tags() {
        let model = Model::new(
            Select::Some(vec!["qname".into(), "pos".into()]),
            Some(vec![("NM".into(), "i".into()), ("MD".into(), "Z".into())]),
            CS,
        )
        .unwrap();
        assert_eq!(model.to_string(), "fields=qname,pos;tags=NM:i,MD:Z");
    }

    #[test]
    fn test_display_zero_half_open() {
        let model = Model::new(
            Select::Some(vec!["qname".into(), "pos".into()]),
            None,
            CoordSystem::ZeroHalfOpen,
        )
        .unwrap();
        assert_eq!(model.to_string(), "fields=qname,pos;coords=01");
    }

    #[test]
    fn test_display_one_closed_omitted() {
        // OneClosed is the default; should not be emitted.
        let model = Model::new(
            Select::Some(vec!["pos".into()]),
            None,
            CoordSystem::OneClosed,
        )
        .unwrap();
        assert_eq!(model.to_string(), "fields=pos");
    }

    #[test]
    fn test_display_tags_no_defs() {
        let model = Model::new(Select::Some(vec!["qname".into()]), Some(vec![]), CS).unwrap();
        assert_eq!(model.to_string(), "fields=qname;tags");
    }

    #[test]
    fn test_from_str_defaults() {
        let model: Model = "fields=*".parse().unwrap();
        assert_eq!(model, Model::default());
    }

    #[test]
    fn test_from_str_roundtrip() {
        let model = Model::new(
            Select::Some(vec!["qname".into(), "pos".into()]),
            Some(vec![("NM".into(), "i".into()), ("MD".into(), "Z".into())]),
            CS,
        )
        .unwrap();
        let s = model.to_string();
        let parsed: Model = s.parse().unwrap();
        assert_eq!(model, parsed);
    }

    #[test]
    fn test_from_str_roundtrip_defaults() {
        let model = Model::default();
        let s = model.to_string();
        let parsed: Model = s.parse().unwrap();
        assert_eq!(model, parsed);
    }

    #[test]
    fn test_from_str_roundtrip_empty_tags() {
        let model = Model::new(Select::Some(vec!["qname".into()]), Some(vec![]), CS).unwrap();
        let s = model.to_string();
        let parsed: Model = s.parse().unwrap();
        assert_eq!(model, parsed);
    }

    #[test]
    fn test_from_str_roundtrip_coord_system() {
        let model = Model::new(
            Select::Some(vec!["qname".into(), "pos".into()]),
            None,
            CoordSystem::ZeroHalfOpen,
        )
        .unwrap();
        let s = model.to_string();
        let parsed: Model = s.parse().unwrap();
        assert_eq!(model, parsed);
        assert_eq!(parsed.coord_system(), CoordSystem::ZeroHalfOpen);
    }

    #[test]
    fn test_clone_eq() {
        let model = Model::new(
            Select::Some(vec!["qname".into()]),
            Some(vec![("NM".into(), "i".into())]),
            CS,
        )
        .unwrap();
        let cloned = model.clone();
        assert_eq!(model, cloned);
    }

    #[test]
    fn test_schema_independence() {
        // Schema should not depend on any file header content.
        let m1 = Model::new(
            Select::Some(vec!["qname".into(), "rname".into(), "pos".into()]),
            None,
            CS,
        )
        .unwrap();
        let m2 = Model::new(
            Select::Some(vec!["qname".into(), "rname".into(), "pos".into()]),
            None,
            CS,
        )
        .unwrap();
        assert_eq!(m1.schema().as_ref(), m2.schema().as_ref());
    }
}

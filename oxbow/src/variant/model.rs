pub mod batch;
pub mod field;
pub mod genotype;
pub mod info;

pub use batch::{BatchBuilder, GenotypeBy};

use std::sync::Arc;

use arrow::datatypes::{DataType, Field as ArrowField, Schema, SchemaRef};

use crate::{OxbowError, Select};
use field::{Field, DEFAULT_FIELD_NAMES};
use genotype::GenotypeDef;
use info::InfoDef;

/// A data model for variant records (VCF/BCF).
///
/// Encapsulates the schema-defining parameters for a variant projection:
/// which standard fields, INFO sub-fields, FORMAT/genotype fields, samples,
/// and genotype layout mode to use.
///
/// - `fields` selects which standard VCF fields become Arrow columns.
///   `None` → all 7 standard fields.
/// - `info_defs` controls the INFO struct column independently.
///   `None` → no info column. `Some(vec![])` → empty struct.
/// - `genotype_defs` + `samples` control per-sample/per-field genotype columns.
///   Both must be `Some` (and non-empty) to produce genotype columns.
/// - `genotype_by` controls the layout: `Sample` (default) or `Field`.
/// - `samples_nested` controls whether genotype columns are wrapped in a
///   single `"samples"` struct column (`true`) or are top-level (`false`,
///   default).
///
/// The model can produce an Arrow schema independently of any file header.
/// Use `from_header()` to derive definitions from a VCF header.
#[derive(Clone, Debug)]
pub struct Model {
    fields: Vec<Field>,
    info_defs: Option<Vec<InfoDef>>,
    genotype_defs: Option<Vec<GenotypeDef>>,
    genotype_by: GenotypeBy,
    samples: Option<Vec<String>>,
    samples_nested: bool,
    schema: SchemaRef,
}

impl Model {
    /// Create a new variant model from validated definitions.
    ///
    /// - `fields`: standard VCF field selection. `All` → all 7 standard fields.
    /// - `info_defs`: validated INFO definitions. `None` → no info column.
    /// - `genotype_defs`: validated FORMAT definitions. `None` → no genotype columns.
    /// - `samples`: sample names. `None` → no genotype columns.
    /// - `genotype_by`: layout mode. Defaults to `GenotypeBy::Sample`.
    /// - `samples_nested`: if `true`, genotype columns are wrapped in a
    ///   single `"samples"` struct column. If `false` (default), they are
    ///   top-level.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        fields: Select<String>,
        info_defs: Option<Vec<InfoDef>>,
        genotype_defs: Option<Vec<GenotypeDef>>,
        genotype_by: Option<GenotypeBy>,
        samples: Option<Vec<String>>,
        samples_nested: Option<bool>,
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

        let genotype_by = genotype_by.unwrap_or(GenotypeBy::Sample);
        let samples_nested = samples_nested.unwrap_or(false);

        let schema = Self::build_schema(
            &parsed_fields,
            info_defs.as_deref(),
            genotype_defs.as_deref(),
            samples.as_deref(),
            &genotype_by,
            samples_nested,
        );

        Ok(Self {
            fields: parsed_fields,
            info_defs,
            genotype_defs,
            genotype_by,
            samples,
            samples_nested,
            schema,
        })
    }

    /// Create a model by deriving INFO and FORMAT definitions from a VCF header.
    ///
    /// - `fields`: standard VCF field selection. `All` → all 7 standard fields.
    /// - `info_field_names`: `All` → all INFO from header. `Some(vec)` → filter
    ///   by name. `Omit` → no info column.
    /// - `genotype_field_names`: `All` → all FORMAT from header. `Some(vec)` →
    ///   filter by name. `Omit` → no genotype columns.
    /// - `samples`: `All` → all samples from header. `Some(vec)` → filter by
    ///   name. `Omit` → no sample columns.
    #[allow(clippy::too_many_arguments)]
    pub fn from_header(
        header: &noodles::vcf::Header,
        fields: Select<String>,
        info_field_names: Select<String>,
        genotype_field_names: Select<String>,
        genotype_by: Option<GenotypeBy>,
        samples: Select<String>,
        samples_nested: Option<bool>,
    ) -> crate::Result<Self> {
        // Derive info defs from header
        // Omit → no info column. All/Some → info column present (even if empty struct).
        let info_defs = match info_field_names {
            Select::Omit => None,
            sel => {
                let names: Vec<String> = match sel {
                    Select::All => header
                        .infos()
                        .iter()
                        .map(|(name, _)| name.to_string())
                        .collect(),
                    Select::Some(names) => names,
                    Select::Omit => unreachable!(),
                };
                let defs: Vec<InfoDef> = names
                    .into_iter()
                    .filter_map(|name| {
                        let info = header.infos().get(&name)?;
                        Some(InfoDef::new(name, &info.number(), &info.ty()))
                    })
                    .collect();
                Some(defs)
            }
        };

        // Derive genotype defs from header
        // Omit → deactivate genotype output. All/Some → genotype active (even if empty).
        let genotype_defs = match genotype_field_names {
            Select::Omit => None,
            sel => {
                let names: Vec<String> = match sel {
                    Select::All => header
                        .formats()
                        .iter()
                        .map(|(name, _)| name.to_string())
                        .collect(),
                    Select::Some(names) => names,
                    Select::Omit => unreachable!(),
                };
                let defs: Vec<GenotypeDef> = names
                    .into_iter()
                    .filter_map(|name| {
                        let format = header.formats().get(&name)?;
                        Some(GenotypeDef::new(name, &format.number(), &format.ty()))
                    })
                    .collect();
                Some(defs)
            }
        };

        // Derive sample names from header
        // Omit → no sample output. All (with samples in header) / Some → sample output active.
        // Select::All with no samples in header → None (nothing to include).
        // Select::Some([]) → Some(vec![]) (column present, empty).
        let samples = match samples {
            Select::Omit => None,
            Select::All => {
                let s: Vec<String> = header.sample_names().iter().cloned().collect();
                if s.is_empty() {
                    None
                } else {
                    Some(s)
                }
            }
            Select::Some(names) => Some(names),
        };

        Self::new(
            fields,
            info_defs,
            genotype_defs,
            genotype_by,
            samples,
            samples_nested,
        )
    }

    fn build_schema(
        fields: &[Field],
        info_defs: Option<&[InfoDef]>,
        genotype_defs: Option<&[GenotypeDef]>,
        samples: Option<&[String]>,
        genotype_by: &GenotypeBy,
        samples_nested: bool,
    ) -> SchemaRef {
        let mut arrow_fields: Vec<ArrowField> =
            fields.iter().map(|f| f.get_arrow_field()).collect();

        // INFO struct column
        if let Some(defs) = info_defs {
            let nested: Vec<ArrowField> = defs.iter().map(|d| d.get_arrow_field()).collect();
            arrow_fields.push(ArrowField::new(
                "info",
                DataType::Struct(arrow::datatypes::Fields::from(nested)),
                true,
            ));
        }

        // Genotype columns: both samples and genotype_defs must be Some to activate.
        // Some([]) produces empty struct content; None deactivates entirely.
        if let (Some(samples), Some(gt_defs)) = (samples, genotype_defs) {
            let genotype_columns: Vec<ArrowField> = match genotype_by {
                GenotypeBy::Sample => samples
                    .iter()
                    .map(|sample_name| {
                        let nested: Vec<ArrowField> = gt_defs
                            .iter()
                            .map(|def| def.get_arrow_field(&def.name))
                            .collect();
                        ArrowField::new(
                            sample_name,
                            DataType::Struct(arrow::datatypes::Fields::from(nested)),
                            true,
                        )
                    })
                    .collect(),
                GenotypeBy::Field => gt_defs
                    .iter()
                    .map(|def| {
                        let nested: Vec<ArrowField> = samples
                            .iter()
                            .map(|sample_name| def.get_arrow_field(sample_name))
                            .collect();
                        ArrowField::new(
                            &def.name,
                            DataType::Struct(arrow::datatypes::Fields::from(nested)),
                            true,
                        )
                    })
                    .collect(),
            };

            if !samples_nested {
                arrow_fields.extend(genotype_columns);
            } else {
                arrow_fields.push(ArrowField::new(
                    "samples",
                    DataType::Struct(arrow::datatypes::Fields::from(genotype_columns)),
                    true,
                ));
            }
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

    /// The validated INFO definitions, if included.
    pub fn info_defs(&self) -> Option<&[InfoDef]> {
        self.info_defs.as_deref()
    }

    /// Whether the INFO struct column is included.
    pub fn has_info(&self) -> bool {
        self.info_defs.is_some()
    }

    /// The validated genotype (FORMAT) definitions, if included.
    pub fn genotype_defs(&self) -> Option<&[GenotypeDef]> {
        self.genotype_defs.as_deref()
    }

    /// The sample names, if included.
    pub fn samples(&self) -> Option<&[String]> {
        self.samples.as_deref()
    }

    /// The genotype layout mode.
    pub fn genotype_by(&self) -> &GenotypeBy {
        &self.genotype_by
    }

    /// Whether genotype columns are nested (wrapped in a `"samples"` struct
    /// column) or top-level.
    pub fn samples_nested(&self) -> bool {
        self.samples_nested
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

        // Fixed fields
        let projected_fields: Vec<String> = self
            .fields
            .iter()
            .filter(|f| columns.iter().any(|c| c.eq_ignore_ascii_case(f.name())))
            .map(|f| f.name().to_string())
            .collect();

        // INFO
        let info_defs = if self.has_info() && columns.iter().any(|c| c.eq_ignore_ascii_case("info"))
        {
            self.info_defs.clone()
        } else {
            None
        };

        // Genotype columns
        let (genotype_defs, samples) = if self.samples_nested {
            // Nested mode: "samples" is an atomic top-level column
            if columns.iter().any(|c| c.eq_ignore_ascii_case("samples")) {
                (self.genotype_defs.clone(), self.samples.clone())
            } else {
                (None, None)
            }
        } else {
            match self.genotype_by {
                GenotypeBy::Sample => {
                    // Unnested: non-fixed, non-info columns are sample names
                    let projected_samples: Option<Vec<String>> = self.samples.as_ref().map(|s| {
                        s.iter()
                            .filter(|name| columns.iter().any(|c| c == *name))
                            .cloned()
                            .collect()
                    });
                    match &projected_samples {
                        Some(s) if s.is_empty() => (None, None),
                        _ => (self.genotype_defs.clone(), projected_samples),
                    }
                }
                GenotypeBy::Field => {
                    // Unnested: non-fixed, non-info columns are FORMAT field names
                    let projected_gt: Option<Vec<GenotypeDef>> =
                        self.genotype_defs.as_ref().map(|defs| {
                            defs.iter()
                                .filter(|def| columns.iter().any(|c| c == &def.name))
                                .cloned()
                                .collect()
                        });
                    match &projected_gt {
                        Some(defs) if defs.is_empty() => (None, None),
                        _ => (projected_gt, self.samples.clone()),
                    }
                }
            }
        };

        Self::new(
            Select::Some(projected_fields),
            info_defs,
            genotype_defs,
            Some(self.genotype_by.clone()),
            samples,
            Some(self.samples_nested),
        )
    }
}

impl PartialEq for Model {
    fn eq(&self, other: &Self) -> bool {
        self.fields == other.fields
            && self.info_defs == other.info_defs
            && self.genotype_defs == other.genotype_defs
            && self.samples == other.samples
            && self.genotype_by == other.genotype_by
            && self.samples_nested == other.samples_nested
    }
}

impl Eq for Model {}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::vcf::header::record::value::map::{Contig, Format, Info, Map};
    use noodles::vcf::Header;

    fn create_test_header() -> Header {
        Header::builder()
            .add_contig("sq0", Map::<Contig>::new())
            .add_info("DP", Map::<Info>::from("DP"))
            .add_format("GT", Map::<Format>::from("GT"))
            .add_sample_name("sample1")
            .add_sample_name("sample2")
            .build()
    }

    #[test]
    fn test_default_model() {
        let model = Model::new(Select::All, None, None, None, None, None).unwrap();
        assert_eq!(model.field_names().len(), 7);
        assert!(!model.has_info());
        assert!(model.genotype_defs().is_none());
        assert!(model.samples().is_none());
        assert_eq!(model.schema().fields().len(), 7);
    }

    #[test]
    fn test_from_header() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::All,
            Select::All,
            None,
            Select::All,
            None,
        )
        .unwrap();
        assert_eq!(model.field_names().len(), 7);
        assert!(model.has_info());
        assert_eq!(model.info_defs().unwrap().len(), 1);
        assert_eq!(model.genotype_defs().unwrap().len(), 1);
        assert_eq!(model.samples().unwrap().len(), 2);
        // 7 fields + info + 2 samples
        assert_eq!(model.schema().fields().len(), 10);
    }

    #[test]
    fn test_from_header_filtered() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::Some(vec!["chrom".into(), "pos".into()]),
            Select::Some(vec!["DP".into()]),
            Select::Some(vec!["GT".into()]),
            None,
            Select::Some(vec!["sample1".into()]),
            None,
        )
        .unwrap();
        assert_eq!(model.field_names(), vec!["chrom", "pos"]);
        assert_eq!(model.info_defs().unwrap().len(), 1);
        assert_eq!(model.genotype_defs().unwrap().len(), 1);
        assert_eq!(model.samples().unwrap(), &["sample1"]);
        // 2 fields + info + 1 sample
        assert_eq!(model.schema().fields().len(), 4);
    }

    #[test]
    fn test_from_header_omit() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::Omit,
            Select::Omit,
            None,
            Select::Omit,
            None,
        )
        .unwrap();
        assert_eq!(model.field_names().len(), 7);
        assert!(!model.has_info());
        assert!(model.genotype_defs().is_none());
        assert!(model.samples().is_none());
        assert_eq!(model.schema().fields().len(), 7);
    }

    #[test]
    fn test_no_info_no_genotype() {
        let model = Model::new(
            Select::Some(vec!["chrom".into(), "pos".into()]),
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();
        assert!(!model.has_info());
        assert!(model.genotype_defs().is_none());
        assert_eq!(model.schema().fields().len(), 2);
    }

    #[test]
    fn test_project_drops_info() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::All,
            Select::All,
            None,
            Select::All,
            None,
        )
        .unwrap();
        let projected = model.project(&["chrom".into(), "pos".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["chrom", "pos"]);
        assert!(!projected.has_info());
        assert!(projected.genotype_defs().is_none());
    }

    #[test]
    fn test_project_keeps_info() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::All,
            Select::All,
            None,
            Select::All,
            None,
        )
        .unwrap();
        let projected = model.project(&["chrom".into(), "info".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["chrom"]);
        assert!(projected.has_info());
    }

    #[test]
    fn test_project_samples() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::All,
            Select::All,
            None,
            Select::All,
            None,
        )
        .unwrap();
        let projected = model.project(&["chrom".into(), "sample1".into()]).unwrap();
        assert_eq!(projected.samples().unwrap(), &["sample1"]);
        assert!(projected.genotype_defs().is_some());
    }

    #[test]
    fn test_invalid_field() {
        let result = Model::new(
            Select::Some(vec!["invalid".into()]),
            None,
            None,
            None,
            None,
            None,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_nested_samples() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::All,
            Select::All,
            None,
            Select::All,
            Some(true),
        )
        .unwrap();
        assert!(model.samples_nested());
        // 7 fields + info + 1 "samples" struct
        assert_eq!(model.schema().fields().len(), 9);
        let samples_field = model.schema().field_with_name("samples").unwrap();
        match samples_field.data_type() {
            DataType::Struct(fields) => {
                assert_eq!(fields.len(), 2); // sample1, sample2
                assert_eq!(fields[0].name(), "sample1");
                assert_eq!(fields[1].name(), "sample2");
            }
            other => panic!("Expected Struct, got {:?}", other),
        }
    }

    #[test]
    fn test_nested_samples_projection() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::All,
            Select::All,
            None,
            Select::All,
            Some(true),
        )
        .unwrap();
        // "samples" is an atomic column in nested mode
        let projected = model.project(&["chrom".into(), "samples".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["chrom"]);
        assert!(projected.samples().is_some());
        assert_eq!(projected.schema().fields().len(), 2); // chrom + samples
        assert_eq!(projected.schema().field(1).name(), "samples");
    }

    #[test]
    fn test_nested_samples_excluded() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::All,
            Select::All,
            None,
            Select::All,
            Some(true),
        )
        .unwrap();
        let projected = model.project(&["chrom".into(), "info".into()]).unwrap();
        assert!(projected.samples().is_none());
        assert!(projected.genotype_defs().is_none());
        assert_eq!(projected.schema().fields().len(), 2); // chrom + info
    }

    #[test]
    fn test_samples_nested_default() {
        let header = create_test_header();
        // Default (samples_nested = false) produces top-level sample columns
        let model = Model::from_header(
            &header,
            Select::All,
            Select::All,
            Select::All,
            None,
            Select::All,
            None,
        )
        .unwrap();
        assert!(!model.samples_nested());
        // 7 fields + info + 2 sample columns
        assert_eq!(model.schema().fields().len(), 10);
    }

    // --- Empty-select (Some([])) edge cases ---

    #[test]
    fn test_info_empty_select_produces_empty_struct_column() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::Some(vec![]), // info present, empty
            Select::Omit,
            None,
            Select::Omit,
            None,
        )
        .unwrap();
        assert!(model.has_info());
        assert_eq!(model.info_defs().unwrap().len(), 0);
        // 7 fields + empty info struct
        assert_eq!(model.schema().fields().len(), 8);
        let info_field = model.schema().field_with_name("info").unwrap();
        match info_field.data_type() {
            DataType::Struct(fields) => assert!(fields.is_empty()),
            other => panic!("Expected empty Struct, got {:?}", other),
        }
    }

    #[test]
    fn test_info_omit_produces_no_column() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::Omit,
            Select::Omit,
            None,
            Select::Omit,
            None,
        )
        .unwrap();
        assert!(!model.has_info());
        assert_eq!(model.schema().fields().len(), 7);
    }

    #[test]
    fn test_genotype_empty_select_by_sample_produces_empty_struct_columns() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::Omit,
            Select::Some(vec![]), // genotype active, no fields
            Some(GenotypeBy::Sample),
            Select::Some(vec!["sample1".into()]),
            None,
        )
        .unwrap();
        assert!(model.genotype_defs().is_some());
        assert_eq!(model.genotype_defs().unwrap().len(), 0);
        // 7 fields + 1 sample column (empty struct)
        assert_eq!(model.schema().fields().len(), 8);
        let sample_field = model.schema().field_with_name("sample1").unwrap();
        match sample_field.data_type() {
            DataType::Struct(fields) => assert!(fields.is_empty()),
            other => panic!("Expected empty Struct, got {:?}", other),
        }
    }

    #[test]
    fn test_genotype_empty_select_by_field_produces_no_columns() {
        // by_field: columns are keyed by FORMAT field name; no fields → no columns
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::Omit,
            Select::Some(vec![]), // genotype active, no fields
            Some(GenotypeBy::Field),
            Select::Some(vec!["sample1".into()]),
            None,
        )
        .unwrap();
        assert_eq!(model.schema().fields().len(), 7);
    }

    #[test]
    fn test_genotype_omit_deactivates_output_despite_samples() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::Omit,
            Select::Omit, // genotype deactivated
            None,
            Select::Some(vec!["sample1".into()]),
            Some(true),
        )
        .unwrap();
        assert!(model.genotype_defs().is_none());
        assert_eq!(model.schema().fields().len(), 7);
    }

    #[test]
    fn test_samples_empty_select_nested_produces_empty_struct_column() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::Omit,
            Select::Some(vec!["GT".into()]),
            None,
            Select::Some(vec![]), // samples active, none selected
            Some(true),
        )
        .unwrap();
        assert!(model.samples().is_some());
        assert_eq!(model.samples().unwrap().len(), 0);
        // 7 fields + "samples" struct (empty: no samples → no sub-fields)
        assert_eq!(model.schema().fields().len(), 8);
        let samples_field = model.schema().field_with_name("samples").unwrap();
        match samples_field.data_type() {
            DataType::Struct(fields) => assert!(fields.is_empty()),
            other => panic!("Expected empty Struct, got {:?}", other),
        }
    }

    #[test]
    fn test_samples_empty_select_unnested_produces_no_columns() {
        // unnested: columns are keyed by sample name; no samples → no columns
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::Omit,
            Select::Some(vec!["GT".into()]),
            None,
            Select::Some(vec![]), // samples active, none selected
            Some(false),
        )
        .unwrap();
        assert_eq!(model.schema().fields().len(), 7);
    }

    #[test]
    fn test_samples_omit_deactivates_nested_column() {
        let header = create_test_header();
        let model = Model::from_header(
            &header,
            Select::All,
            Select::Omit,
            Select::Some(vec!["GT".into()]),
            None,
            Select::Omit, // samples deactivated
            Some(true),
        )
        .unwrap();
        assert!(model.samples().is_none());
        assert_eq!(model.schema().fields().len(), 7);
    }
}

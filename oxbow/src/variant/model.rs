pub mod batch_builder;
pub mod field;
pub mod genotype;
pub mod info;

pub use batch_builder::{BatchBuilder, GenotypeBy};

use std::sync::Arc;

use arrow::datatypes::{DataType, Field as ArrowField, Schema, SchemaRef};

use crate::OxbowError;
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
///
/// The model can produce an Arrow schema independently of any file header.
/// Use `from_header()` to derive definitions from a VCF header.
#[derive(Clone, Debug)]
pub struct Model {
    fields: Vec<Field>,
    info_defs: Option<Vec<InfoDef>>,
    genotype_defs: Option<Vec<GenotypeDef>>,
    samples: Option<Vec<String>>,
    genotype_by: GenotypeBy,
    schema: SchemaRef,
}

impl Model {
    /// Create a new variant model from validated definitions.
    ///
    /// - `fields`: standard VCF field names. `None` → all 7 standard fields.
    /// - `info_defs`: validated INFO definitions. `None` → no info column.
    /// - `genotype_defs`: validated FORMAT definitions. `None` → no genotype columns.
    /// - `samples`: sample names. `None` → no genotype columns.
    /// - `genotype_by`: layout mode. Defaults to `GenotypeBy::Sample`.
    pub fn new(
        fields: Option<Vec<String>>,
        info_defs: Option<Vec<InfoDef>>,
        genotype_defs: Option<Vec<GenotypeDef>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<GenotypeBy>,
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

        let genotype_by = genotype_by.unwrap_or(GenotypeBy::Sample);

        let schema = Self::build_schema(
            &parsed_fields,
            info_defs.as_deref(),
            genotype_defs.as_deref(),
            samples.as_deref(),
            &genotype_by,
        );

        Ok(Self {
            fields: parsed_fields,
            info_defs,
            genotype_defs,
            samples,
            genotype_by,
            schema,
        })
    }

    /// Create a model by deriving INFO and FORMAT definitions from a VCF header.
    ///
    /// Explicit name lists filter which definitions are included.
    /// `None` → include all definitions from the header.
    #[allow(clippy::too_many_arguments)]
    pub fn from_header(
        header: &noodles::vcf::Header,
        fields: Option<Vec<String>>,
        info_field_names: Option<Vec<String>>,
        genotype_field_names: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<GenotypeBy>,
    ) -> crate::Result<Self> {
        // Derive info defs from header
        let info_names: Vec<String> = info_field_names.unwrap_or_else(|| {
            header
                .infos()
                .iter()
                .map(|(name, _)| name.to_string())
                .collect()
        });
        let info_defs: Vec<InfoDef> = info_names
            .into_iter()
            .filter_map(|name| {
                let info = header.infos().get(&name)?;
                Some(InfoDef::new(name, &info.number(), &info.ty()))
            })
            .collect();

        // Derive genotype defs from header
        let gt_names: Vec<String> = genotype_field_names.unwrap_or_else(|| {
            header
                .formats()
                .iter()
                .map(|(name, _)| name.to_string())
                .collect()
        });
        let genotype_defs: Vec<GenotypeDef> = gt_names
            .into_iter()
            .filter_map(|name| {
                let format = header.formats().get(&name)?;
                Some(GenotypeDef::new(name, &format.number(), &format.ty()))
            })
            .collect();

        // Derive sample names from header
        let samples = samples.unwrap_or_else(|| header.sample_names().iter().cloned().collect());

        // Wrap in Option: empty → None for info/genotype, Some for samples
        let info_defs = if info_defs.is_empty() {
            None
        } else {
            Some(info_defs)
        };
        let genotype_defs = if genotype_defs.is_empty() {
            None
        } else {
            Some(genotype_defs)
        };
        let samples = if samples.is_empty() {
            None
        } else {
            Some(samples)
        };

        Self::new(fields, info_defs, genotype_defs, samples, genotype_by)
    }

    fn build_schema(
        fields: &[Field],
        info_defs: Option<&[InfoDef]>,
        genotype_defs: Option<&[GenotypeDef]>,
        samples: Option<&[String]>,
        genotype_by: &GenotypeBy,
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

        // Genotype columns (require both samples and genotype_defs)
        let samples = samples.unwrap_or(&[]);
        let gt_defs = genotype_defs.unwrap_or(&[]);
        if !samples.is_empty() && !gt_defs.is_empty() {
            match genotype_by {
                GenotypeBy::Sample => {
                    // One struct column per sample, with FORMAT fields as sub-fields
                    for sample_name in samples {
                        let nested: Vec<ArrowField> = gt_defs
                            .iter()
                            .map(|def| def.get_arrow_field(&def.name))
                            .collect();
                        arrow_fields.push(ArrowField::new(
                            sample_name,
                            DataType::Struct(arrow::datatypes::Fields::from(nested)),
                            true,
                        ));
                    }
                }
                GenotypeBy::Field => {
                    // One struct column per FORMAT field, with sample names as sub-fields
                    for def in gt_defs {
                        let nested: Vec<ArrowField> = samples
                            .iter()
                            .map(|sample_name| def.get_arrow_field(sample_name))
                            .collect();
                        arrow_fields.push(ArrowField::new(
                            &def.name,
                            DataType::Struct(arrow::datatypes::Fields::from(nested)),
                            true,
                        ));
                    }
                }
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
        let (genotype_defs, samples) = match self.genotype_by {
            GenotypeBy::Sample => {
                // Non-fixed, non-info columns are sample names
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
                // Non-fixed, non-info columns are FORMAT field names
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
        };

        Self::new(
            Some(projected_fields),
            info_defs,
            genotype_defs,
            samples,
            Some(self.genotype_by.clone()),
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
        let model = Model::new(None, None, None, None, None).unwrap();
        assert_eq!(model.field_names().len(), 7);
        assert!(!model.has_info());
        assert!(model.genotype_defs().is_none());
        assert!(model.samples().is_none());
        assert_eq!(model.schema().fields().len(), 7);
    }

    #[test]
    fn test_from_header() {
        let header = create_test_header();
        let model = Model::from_header(&header, None, None, None, None, None).unwrap();
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
            Some(vec!["chrom".into(), "pos".into()]),
            Some(vec!["DP".into()]),
            Some(vec!["GT".into()]),
            Some(vec!["sample1".into()]),
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
    fn test_no_info_no_genotype() {
        let model = Model::new(
            Some(vec!["chrom".into(), "pos".into()]),
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
        let model = Model::from_header(&header, None, None, None, None, None).unwrap();
        let projected = model.project(&["chrom".into(), "pos".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["chrom", "pos"]);
        assert!(!projected.has_info());
        assert!(projected.genotype_defs().is_none());
    }

    #[test]
    fn test_project_keeps_info() {
        let header = create_test_header();
        let model = Model::from_header(&header, None, None, None, None, None).unwrap();
        let projected = model.project(&["chrom".into(), "info".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["chrom"]);
        assert!(projected.has_info());
    }

    #[test]
    fn test_project_samples() {
        let header = create_test_header();
        let model = Model::from_header(&header, None, None, None, None, None).unwrap();
        let projected = model.project(&["chrom".into(), "sample1".into()]).unwrap();
        assert_eq!(projected.samples().unwrap(), &["sample1"]);
        assert!(projected.genotype_defs().is_some());
    }

    #[test]
    fn test_invalid_field() {
        let result = Model::new(Some(vec!["invalid".into()]), None, None, None, None);
        assert!(result.is_err());
    }
}

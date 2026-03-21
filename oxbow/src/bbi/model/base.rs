pub mod batch;
pub mod field;

use std::sync::Arc;

use arrow::datatypes::{Field as ArrowField, Schema, SchemaRef};

pub use crate::bed::model::schema::BedSchema;
use crate::{CoordSystem, OxbowError, Select};
pub use batch::BatchBuilder;
use field::{bed_standard_fields, FieldDef};

pub struct BigBedRecord<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
    pub rest: &'a String,
}

impl<'a> BigBedRecord<'a> {
    pub fn new(chrom: &'a str, entry: &'a bigtools::BedEntry) -> Self {
        Self {
            chrom,
            start: entry.start,
            end: entry.end,
            rest: &entry.rest,
        }
    }
}

pub struct BigWigRecord<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

impl<'a> BigWigRecord<'a> {
    pub fn new(chrom: &'a str, entry: &'a bigtools::Value) -> Self {
        Self {
            chrom,
            start: entry.start,
            end: entry.end,
            value: entry.value,
        }
    }
}

/// Build the full list of FieldDefs for this BedSchema using BBI types.
fn bed_schema_field_defs(bed_schema: &BedSchema) -> Vec<FieldDef> {
    let standard = bed_standard_fields();
    let mut defs: Vec<FieldDef> = standard
        .into_iter()
        .take(bed_schema.standard_field_count())
        .collect();
    defs.extend(bed_schema.custom_fields().iter().cloned());
    defs
}

/// A data model for BBI base (BigBed/BigWig) records.
///
/// Wraps a [`BedSchema`] with field projection, using AutoSql-based Arrow
/// types (e.g., UInt32 for positions) rather than BED's noodles-based types.
///
/// `coord_system` controls the coordinate system to return the positions in.
/// The default is 0-based half-open (BED/BBI convention).
#[derive(Clone, Debug)]
pub struct Model {
    bed_schema: BedSchema,
    fields: Vec<FieldDef>,
    coord_system: CoordSystem,
    schema: SchemaRef,
}

impl Model {
    /// Create a new BBI base model.
    ///
    /// - `bed_schema`: the parsing interpretation.
    /// - `fields`: column names to project. `None` → all fields from the schema.
    pub fn new(
        bed_schema: BedSchema,
        fields: Select<String>,
        coord_system: CoordSystem,
    ) -> crate::Result<Self> {
        let all_defs = bed_schema_field_defs(&bed_schema);
        let projected = match fields {
            Select::All => all_defs,
            Select::Omit => Vec::new(),
            Select::Some(names) => {
                let mut projected = Vec::new();
                for name in &names {
                    let def = all_defs
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
            coord_system,
            schema,
        })
    }

    pub fn bed_schema(&self) -> &BedSchema {
        &self.bed_schema
    }

    pub fn bed_schema_field_defs(&self) -> Vec<FieldDef> {
        bed_schema_field_defs(&self.bed_schema)
    }

    pub fn fields(&self) -> &[FieldDef] {
        &self.fields
    }

    pub fn field_names(&self) -> Vec<String> {
        self.fields.iter().map(|d| d.name.clone()).collect()
    }

    pub fn coord_system(&self) -> CoordSystem {
        self.coord_system
    }

    pub fn schema(&self) -> &SchemaRef {
        &self.schema
    }

    pub fn column_names(&self) -> Vec<String> {
        self.schema
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect()
    }

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

        Self::new(
            self.bed_schema.clone(),
            Select::Some(projected),
            self.coord_system,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bed::FieldType;
    use arrow::datatypes::DataType;

    #[test]
    fn test_bedgraph_model() {
        let bed_schema = BedSchema::new_bedgraph().unwrap();
        let model = Model::new(bed_schema, Select::All, CoordSystem::ZeroHalfOpen).unwrap();
        assert_eq!(model.field_names(), vec!["chrom", "start", "end", "value"]);
        // BBI uses UInt32 for positions (AutoSql types)
        assert_eq!(model.schema().field(1).data_type(), &DataType::UInt32);
        assert_eq!(model.schema().field(3).data_type(), &DataType::Float32);
    }

    #[test]
    fn test_bed6_model() {
        let bed_schema: BedSchema = "bed6".parse().unwrap();
        let model = Model::new(bed_schema, Select::All, CoordSystem::ZeroHalfOpen).unwrap();
        assert_eq!(model.field_names().len(), 6);
        assert_eq!(model.schema().field(1).data_type(), &DataType::UInt32);
    }

    #[test]
    fn test_bed6_projection() {
        let bed_schema: BedSchema = "bed6".parse().unwrap();
        let model = Model::new(bed_schema, Select::All, CoordSystem::ZeroHalfOpen).unwrap();
        let projected = model.project(&["chrom".into(), "strand".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["chrom", "strand"]);
    }

    #[test]
    fn test_custom_fields() {
        let defs = vec![FieldDef::new("extra".into(), FieldType::Float)];
        let bed_schema = BedSchema::new(3, Some(defs)).unwrap();
        let model = Model::new(bed_schema, Select::All, CoordSystem::ZeroHalfOpen).unwrap();
        assert_eq!(model.field_names(), vec!["chrom", "start", "end", "extra"]);
        assert_eq!(model.schema().field(3).data_type(), &DataType::Float32);
    }

    #[test]
    fn test_project_noncontiguous_custom() {
        let defs = vec![
            FieldDef::new("a".into(), FieldType::String),
            FieldDef::new("b".into(), FieldType::Float),
            FieldDef::new("c".into(), FieldType::String),
        ];
        let bed_schema = BedSchema::new(3, Some(defs)).unwrap();
        let model = Model::new(bed_schema, Select::All, CoordSystem::ZeroHalfOpen).unwrap();
        let projected = model.project(&["chrom".into(), "c".into()]).unwrap();
        assert_eq!(projected.field_names(), vec!["chrom", "c"]);
    }
}

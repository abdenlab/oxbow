use std::io;

use bigtools::bed::autosql::parse::Field as AutosqlField;
use bigtools::bed::autosql::parse::FieldType as AutosqlType;

// Re-export shared types from the bed module.
pub use crate::bed::model::field_def::{
    bed_standard_fields, FieldBuilder, FieldDef, FieldType, Push,
};

/// Convert an AutoSql field to a [`FieldDef`].
impl TryFrom<&AutosqlField> for FieldDef {
    type Error = io::Error;

    fn try_from(field: &AutosqlField) -> Result<Self, Self::Error> {
        let ty = FieldType::try_from(field)?;
        Ok(Self {
            name: field.name.clone(),
            ty,
        })
    }
}

/// Convert an AutoSql field type to a [`FieldType`].
impl TryFrom<&AutosqlField> for FieldType {
    type Error = io::Error;

    fn try_from(f: &AutosqlField) -> Result<Self, Self::Error> {
        let ty = match &f.field_type {
            AutosqlType::Byte => match &f.field_size {
                None => Self::Byte,
                Some(size_str) => match size_str.parse::<usize>() {
                    Ok(size) => Self::ByteFixedSizeList(size),
                    Err(_) => Self::ByteList,
                },
            },
            AutosqlType::Ubyte => match &f.field_size {
                None => Self::Ubyte,
                Some(size_str) => match size_str.parse::<usize>() {
                    Ok(size) => Self::UbyteFixedSizeList(size),
                    Err(_) => Self::UbyteList,
                },
            },
            AutosqlType::Short => match &f.field_size {
                None => Self::Short,
                Some(size_str) => match size_str.parse::<usize>() {
                    Ok(size) => Self::ShortFixedSizeList(size),
                    Err(_) => Self::ShortList,
                },
            },
            AutosqlType::Ushort => match &f.field_size {
                None => Self::Ushort,
                Some(size_str) => match size_str.parse::<usize>() {
                    Ok(size) => Self::UshortFixedSizeList(size),
                    Err(_) => Self::UshortList,
                },
            },
            AutosqlType::Int => match &f.field_size {
                None => Self::Int,
                Some(size_str) => match size_str.parse::<usize>() {
                    Ok(size) => Self::IntFixedSizeList(size),
                    Err(_) => Self::IntList,
                },
            },
            AutosqlType::Uint => match &f.field_size {
                None => Self::Uint,
                Some(size_str) => match size_str.parse::<usize>() {
                    Ok(size) => Self::UintFixedSizeList(size),
                    Err(_) => Self::UintList,
                },
            },
            AutosqlType::Bigint => match &f.field_size {
                None => Self::Bigint,
                Some(size_str) => match size_str.parse::<usize>() {
                    Ok(size) => Self::BigintFixedSizeList(size),
                    Err(_) => Self::BigintList,
                },
            },
            AutosqlType::Float => match &f.field_size {
                None => Self::Float,
                Some(size_str) => match size_str.parse::<usize>() {
                    Ok(size) => Self::FloatFixedSizeList(size),
                    Err(_) => Self::FloatList,
                },
            },
            AutosqlType::Double => match &f.field_size {
                None => Self::Double,
                Some(size_str) => match size_str.parse::<usize>() {
                    Ok(size) => Self::DoubleFixedSizeList(size),
                    Err(_) => Self::DoubleList,
                },
            },
            AutosqlType::Char => Self::Char,
            AutosqlType::String => Self::String,
            AutosqlType::Lstring => Self::Lstring,
            AutosqlType::Enum(items) => Self::Enum(items.clone()),
            AutosqlType::Set(items) => Self::Set(items.clone()),
            AutosqlType::Declaration(_, _) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "Unsupported AutoSql type: Declaration.",
                ));
            }
        };
        Ok(ty)
    }
}

/// The 12 standard BED field definitions as AutoSql fields.
///
/// This is the bigtools-dependent version, kept for compatibility with BBI code
/// that receives AutoSql fields from bigtools parsers.
pub fn bed_standard_autosql_fields() -> [AutosqlField; 12] {
    [
        AutosqlField {
            name: "chrom".to_string(),
            field_type: AutosqlType::String,
            field_size: None,
            index_type: None,
            auto: false,
            comment: "Reference sequence chromosome or scaffold.".to_string(),
        },
        AutosqlField {
            name: "start".to_string(),
            field_type: AutosqlType::Uint,
            field_size: None,
            index_type: None,
            auto: false,
            comment: "Start position in chromosome.".to_string(),
        },
        AutosqlField {
            name: "end".to_string(),
            field_type: AutosqlType::Uint,
            field_size: None,
            index_type: None,
            auto: false,
            comment: "End position in chromosome.".to_string(),
        },
        AutosqlField {
            name: "name".to_string(),
            field_type: AutosqlType::String,
            field_size: None,
            index_type: None,
            auto: false,
            comment: "Name or ID of item.".to_string(),
        },
        AutosqlField {
            name: "score".to_string(),
            field_type: AutosqlType::Ushort,
            field_size: None,
            index_type: None,
            auto: false,
            comment: "Score (0-1000).".to_string(),
        },
        AutosqlField {
            name: "strand".to_string(),
            field_type: AutosqlType::Char,
            field_size: None,
            index_type: None,
            auto: false,
            comment: "+ or - for strand.".to_string(),
        },
        AutosqlField {
            name: "thickStart".to_string(),
            field_type: AutosqlType::Uint,
            field_size: None,
            index_type: None,
            auto: false,
            comment: "Start of where display should be thick (start codon).".to_string(),
        },
        AutosqlField {
            name: "thickEnd".to_string(),
            field_type: AutosqlType::Uint,
            field_size: None,
            index_type: None,
            auto: false,
            comment: "End of where display should be thick (stop codon).".to_string(),
        },
        AutosqlField {
            name: "itemRgb".to_string(),
            field_type: AutosqlType::Ubyte,
            field_size: Some("3".to_string()),
            index_type: None,
            auto: false,
            comment: "RGB value (use R,G,B string in input file).".to_string(),
        },
        AutosqlField {
            name: "blockCount".to_string(),
            field_type: AutosqlType::Uint,
            field_size: None,
            index_type: None,
            auto: false,
            comment: "Number of blocks.".to_string(),
        },
        AutosqlField {
            name: "blockSizes".to_string(),
            field_type: AutosqlType::Uint,
            field_size: Some("blockCount".to_string()),
            index_type: None,
            auto: false,
            comment: "Comma separated list of block sizes.".to_string(),
        },
        AutosqlField {
            name: "blockStarts".to_string(),
            field_type: AutosqlType::Uint,
            field_size: Some("blockCount".to_string()),
            index_type: None,
            auto: false,
            comment: "Start positions relative to chromStart.".to_string(),
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fielddef_try_from_autosql_field() {
        let autosql_field = AutosqlField {
            name: "foo".to_string(),
            field_type: AutosqlType::Int,
            field_size: None,
            index_type: None,
            auto: false,
            comment: "Just a field".to_string(),
        };
        let field_def = FieldDef::try_from(&autosql_field);
        assert!(field_def.is_ok());
        let field_def = field_def.unwrap();
        assert_eq!(field_def.name, "foo");
        assert_eq!(field_def.ty, FieldType::Int);
    }
}

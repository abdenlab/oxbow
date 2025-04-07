use std::io;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::ArrayRef;
use arrow::array::{
    FixedSizeListBuilder, Float32Builder, Float64Builder, GenericStringBuilder, Int16Builder,
    Int32Builder, Int64Builder, Int8Builder, LargeStringBuilder, ListBuilder,
    StringDictionaryBuilder, UInt16Builder, UInt32Builder, UInt64Builder, UInt8Builder,
};
use arrow::datatypes::Field as ArrowField;
use arrow::datatypes::{DataType, Int32Type};

use bigtools::bed::autosql::parse::Field as AutosqlField;
use bigtools::bed::autosql::parse::FieldType as AutosqlType;

pub fn bed_standard_fields() -> [AutosqlField; 12] {
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

/// An AutoSql-compatible field definition.
///
/// We flatten AutoSql's type/container hierarchy to a [`FieldType`] representation that maps to
/// Arrow structures.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct FieldDef {
    pub name: String,
    pub ty: FieldType,
}

impl FieldDef {
    pub fn new(name: String, ty: FieldType) -> Self {
        Self { name, ty }
    }

    pub fn get_arrow_field(&self) -> ArrowField {
        ArrowField::new(&self.name, self.ty.arrow_type(), true)
    }
}

impl TryFrom<(String, String)> for FieldDef {
    type Error = io::Error;

    fn try_from(def: (String, String)) -> Result<Self, Self::Error> {
        let (name, type_str) = def;
        let ty: FieldType = type_str.parse().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid tag type: {}", type_str),
            )
        })?;
        Ok(Self { name, ty })
    }
}

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

/// A mapping of AutoSql types to Arrow data types.
///
/// Scalar types and fixed-size array types map naturally. Parameterized-length array types in
/// AutoSql are mapped to variable-sized lists in Arrow with no parameter. Enum and set types in
/// AutoSql are mapped to Dictionary and list-of-Dictionary data types in Arrow, respectively.
///
/// # References
/// - [AutoSql and AutoXml Linux Journal article](https://www.linuxjournal.com/article/5949)
/// - [History of AutoSql](https://genomewiki.ucsc.edu/index.php/AutoSql)
/// - [Arrow DataTypes](https://docs.rs/arrow/latest/arrow/datatypes/enum.DataType.html)
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum FieldType {
    Byte,
    ByteList,
    ByteFixedSizeList(usize),
    Ubyte,
    UbyteList,
    UbyteFixedSizeList(usize),
    Short,
    ShortList,
    ShortFixedSizeList(usize),
    Ushort,
    UshortList,
    UshortFixedSizeList(usize),
    Int,
    IntList,
    IntFixedSizeList(usize),
    Uint,
    UintList,
    UintFixedSizeList(usize),
    Bigint,
    BigintList,
    BigintFixedSizeList(usize),
    Ubigint,
    UbigintList,
    UbigintFixedSizeList(usize),
    Float,
    FloatList,
    FloatFixedSizeList(usize),
    Double,
    DoubleList,
    DoubleFixedSizeList(usize),
    Char,
    String,
    Lstring,
    Enum(Vec<String>),
    Set(Vec<String>),
}

impl FieldType {
    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::Byte => DataType::Int8,
            Self::ByteList => {
                let item = ArrowField::new("item", DataType::Int8, true);
                DataType::List(Arc::new(item))
            }
            Self::ByteFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::Int8, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }
            Self::Ubyte => DataType::UInt8,
            Self::UbyteList => {
                let item = ArrowField::new("item", DataType::UInt8, true);
                DataType::List(Arc::new(item))
            }
            Self::UbyteFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::UInt8, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }
            Self::Short => DataType::Int16,
            Self::ShortList => {
                let item = ArrowField::new("item", DataType::Int16, true);
                DataType::List(Arc::new(item))
            }
            Self::ShortFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::Int16, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }
            Self::Ushort => DataType::UInt16,
            Self::UshortList => {
                let item = ArrowField::new("item", DataType::UInt16, true);
                DataType::List(Arc::new(item))
            }
            Self::UshortFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::UInt16, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }
            Self::Int => DataType::Int32,
            Self::IntList => {
                let item = ArrowField::new("item", DataType::Int32, true);
                DataType::List(Arc::new(item))
            }
            Self::IntFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::Int32, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }
            Self::Uint => DataType::UInt32,
            Self::UintList => {
                let item = ArrowField::new("item", DataType::UInt32, true);
                DataType::List(Arc::new(item))
            }
            Self::UintFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::UInt32, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }
            Self::Bigint => DataType::Int64,
            Self::BigintList => {
                let item = ArrowField::new("item", DataType::Int64, true);
                DataType::List(Arc::new(item))
            }
            Self::BigintFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::Int64, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }
            Self::Ubigint => DataType::UInt64,
            Self::UbigintList => {
                let item = ArrowField::new("item", DataType::UInt64, true);
                DataType::List(Arc::new(item))
            }
            Self::UbigintFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::UInt64, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }
            Self::Float => DataType::Float32,
            Self::FloatList => {
                let item = ArrowField::new("item", DataType::Float32, true);
                DataType::List(Arc::new(item))
            }
            Self::FloatFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::Float32, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }
            Self::Double => DataType::Float64,
            Self::DoubleList => {
                let item = ArrowField::new("item", DataType::Float64, true);
                DataType::List(Arc::new(item))
            }
            Self::DoubleFixedSizeList(n) => {
                let item = ArrowField::new("item", DataType::Float64, true);
                DataType::FixedSizeList(Arc::new(item), *n as i32)
            }
            Self::Char => DataType::Utf8,
            Self::String => DataType::Utf8,
            Self::Lstring => DataType::LargeUtf8,
            Self::Enum(_) => {
                DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8))
            }
            Self::Set(_) => {
                let item = ArrowField::new(
                    "item",
                    DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8)),
                    true,
                );
                DataType::List(Arc::new(item))
            }
        }
    }
}

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

impl FromStr for FieldType {
    type Err = std::io::Error;

    /// Parse a string representation into a [`FieldType`].
    ///
    /// A field type can be:
    /// - A numerical scalar type
    /// - A variable-length list of numbers
    /// - A fixed-size list of numbers
    /// - A string type (`char`, `string`, `lstring`)
    /// - An `enum` or `set` (list of enum) type
    ///
    /// We support two different notations (case-insensitive):
    ///
    /// AutoSql notation:
    /// - `byte`, `byte[]`, `byte[10]`
    /// - `ubyte`, `ubyte[]`, `ubyte[10]`
    /// - `short`, `short[]`, `short[10]`
    /// - `ushort`, `ushort[]`, `ushort[10]`
    /// - `int`, `int[]`, `int[10]`
    /// - `uint`, `uint[]`, `uint[10]`
    /// - `bigint`, `bigint[]`, `bigint[10]`
    /// - `ubigint`, `ubigint[]`, `ubigint[10]`
    /// - `float`, `float[]`, `float[10]`
    /// - `double`, `double[]`, `double[10]`
    /// - `char`, `string`, `lstring`
    /// - `enum(item1,item2,item3)`, `set(item1,item2,item3)`
    ///
    /// Rust notation for numerical types and arrays:
    /// - `i8`, `[i8]`, `[i8; 10]`
    /// - `u8`, `[u8]`, `[u8; 10]`
    /// - `i16`, `[i16]`, `[i16; 10]`
    /// - `u16`, `[u16]`, `[u16; 10]`
    /// - `i32`, `[i32]`, `[i32; 10]`
    /// - `u32`, `[u32]`, `[u32; 10]`
    /// - `i64`, `[i64]`, `[i64; 10]`
    /// - `u64`, `[u64]`, `[u64; 10]`
    /// - `f32`, `[f32]`, `[f32; 10]`
    /// - `f64`, `[f64]`, `[f64; 10]`
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_lowercase();

        match s.as_str() {
            "byte" | "i8" => Ok(Self::Byte),
            "byte[]" | "[i8]" => Ok(Self::ByteList),
            s if s.starts_with("byte[") => match parse_size_param(s, "byte[", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            s if s.starts_with("[i8;") => match parse_size_param(s, "[i8;", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            "ubyte" | "u8" => Ok(Self::Ubyte),
            "ubyte[]" | "[u8]" => Ok(Self::UbyteList),
            s if s.starts_with("ubyte[") => match parse_size_param(s, "ubyte[", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            s if s.starts_with("[u8;") => match parse_size_param(s, "[u8;", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            "short" | "i16" => Ok(Self::Short),
            "short[]" | "[i16]" => Ok(Self::ShortList),
            s if s.starts_with("short[") => match parse_size_param(s, "short[", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            s if s.starts_with("[i16;") => match parse_size_param(s, "[i16;", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            "ushort" | "u16" => Ok(Self::Ushort),
            "ushort[]" | "[u16]" => Ok(Self::UshortList),
            s if s.starts_with("ushort[") => match parse_size_param(s, "ushort[", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            s if s.starts_with("[u16;") => match parse_size_param(s, "[u16;", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            "int" | "i32" => Ok(Self::Int),
            "int[]" | "[i32]" => Ok(Self::IntList),
            s if s.starts_with("int[") => match parse_size_param(s, "int[", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            s if s.starts_with("[i32;") => match parse_size_param(s, "[i32;", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            "uint" | "u32" => Ok(Self::Uint),
            "uint[]" | "[u32]" => Ok(Self::UintList),
            s if s.starts_with("uint[") => match parse_size_param(s, "uint[", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            s if s.starts_with("[u32;") => match parse_size_param(s, "[u32;", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            "bigint" | "i64" => Ok(Self::Bigint),
            "bigint[]" | "[i64]" => Ok(Self::BigintList),
            s if s.starts_with("bigint[") => match parse_size_param(s, "bigint[", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            s if s.starts_with("[i64;") => match parse_size_param(s, "[i64;", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            "ubigint" | "u64" => Ok(Self::Ubigint),
            "ubigint[]" | "[u64]" => Ok(Self::UbigintList),
            s if s.starts_with("ubigint[") => match parse_size_param(s, "ubigint[", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            s if s.starts_with("[u64;") => match parse_size_param(s, "[u64;", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            "float" | "f32" => Ok(Self::Float),
            "float[]" | "[f32]" => Ok(Self::FloatList),
            s if s.starts_with("float[") => match parse_size_param(s, "float[", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            s if s.starts_with("[f32;") => match parse_size_param(s, "[f32;", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            "double" | "f64" => Ok(Self::Double),
            "double[]" | "[f64]" => Ok(Self::DoubleList),
            s if s.starts_with("double[") => match parse_size_param(s, "double[", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            s if s.starts_with("[f64;") => match parse_size_param(s, "[f64;", "]") {
                Ok(size) => Ok(Self::ByteFixedSizeList(size)),
                Err(e) => Err(e),
            },
            "char" => Ok(Self::Char),
            "string" => Ok(Self::String),
            "lstring" => Ok(Self::Lstring),
            s if s.starts_with("enum(") => match parse_enum_param(s, "enum(", ")") {
                Ok(items) => Ok(Self::Enum(items)),
                Err(e) => Err(e),
            },
            s if s.starts_with("set(") => match parse_enum_param(s, "set(", ")") {
                Ok(items) => Ok(Self::Set(items)),
                Err(e) => Err(e),
            },
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid AutoSql type: {}", s),
            )),
        }
    }
}

fn parse_size_param(s: &str, prefix: &str, suffix: &str) -> io::Result<usize> {
    if s.starts_with(prefix) && s.ends_with(suffix) {
        let size_str = &s[prefix.len()..s.len() - suffix.len()];
        if let Ok(size) = size_str.parse::<usize>() {
            if size > 0 {
                return Ok(size);
            }
        }
    };
    Err(std::io::Error::new(
        std::io::ErrorKind::InvalidInput,
        format!("Invalid fixed-size list format: {}", s),
    ))
}

fn parse_enum_param(s: &str, prefix: &str, suffix: &str) -> io::Result<Vec<String>> {
    if s.starts_with(prefix) && s.ends_with(suffix) {
        let items: Vec<String> = s[prefix.len()..s.len() - suffix.len()]
            .split(',')
            .map(|item| item.trim().to_string())
            .filter(|item| !item.is_empty())
            .collect();
        return Ok(items);
    };
    Err(std::io::Error::new(
        std::io::ErrorKind::InvalidInput,
        format!("Invalid enum/set format: {}", s),
    ))
}

/// A builder for an Arrow array (column) corresponding to an AutoSql field.
pub enum FieldBuilder {
    Byte(Int8Builder),
    ByteList(ListBuilder<Int8Builder>),
    ByteFixedSizeList(FixedSizeListBuilder<Int8Builder>),

    Ubyte(UInt8Builder),
    UbyteList(ListBuilder<UInt8Builder>),
    UbyteFixedSizeList(FixedSizeListBuilder<UInt8Builder>),

    Short(Int16Builder),
    ShortList(ListBuilder<Int16Builder>),
    ShortFixedSizeList(FixedSizeListBuilder<Int16Builder>),

    Ushort(UInt16Builder),
    UshortList(ListBuilder<UInt16Builder>),
    UshortFixedSizeList(FixedSizeListBuilder<UInt16Builder>),

    Int(Int32Builder),
    IntList(ListBuilder<Int32Builder>),
    IntFixedSizeList(FixedSizeListBuilder<Int32Builder>),

    Uint(UInt32Builder),
    UintList(ListBuilder<UInt32Builder>),
    UintFixedSizeList(FixedSizeListBuilder<UInt32Builder>),

    Bigint(Int64Builder),
    BigintList(ListBuilder<Int64Builder>),
    BigintFixedSizeList(FixedSizeListBuilder<Int64Builder>),

    Ubigint(UInt64Builder),
    UbigintList(ListBuilder<UInt64Builder>),
    UbigintFixedSizeList(FixedSizeListBuilder<UInt64Builder>),

    Float(Float32Builder),
    FloatList(ListBuilder<Float32Builder>),
    FloatFixedSizeList(FixedSizeListBuilder<Float32Builder>),

    Double(Float64Builder),
    DoubleList(ListBuilder<Float64Builder>),
    DoubleFixedSizeList(FixedSizeListBuilder<Float64Builder>),

    Char(GenericStringBuilder<i32>),
    String(GenericStringBuilder<i32>),
    Lstring(LargeStringBuilder),

    Enum(StringDictionaryBuilder<Int32Type>),
    Set(ListBuilder<StringDictionaryBuilder<Int32Type>>),
}

impl FieldBuilder {
    pub fn new(field_type: &FieldType, capacity: usize) -> io::Result<Self> {
        let builder = match field_type {
            FieldType::Byte => Self::Byte(Int8Builder::with_capacity(capacity)),
            FieldType::ByteList => {
                Self::ByteList(ListBuilder::with_capacity(Int8Builder::new(), capacity))
            }
            FieldType::ByteFixedSizeList(size) => Self::ByteFixedSizeList(
                FixedSizeListBuilder::with_capacity(Int8Builder::new(), *size as i32, capacity),
            ),
            FieldType::Ubyte => Self::Ubyte(UInt8Builder::with_capacity(capacity)),
            FieldType::UbyteList => {
                Self::UbyteList(ListBuilder::with_capacity(UInt8Builder::new(), capacity))
            }
            FieldType::UbyteFixedSizeList(size) => Self::UbyteFixedSizeList(
                FixedSizeListBuilder::with_capacity(UInt8Builder::new(), *size as i32, capacity),
            ),
            FieldType::Short => Self::Short(Int16Builder::with_capacity(capacity)),
            FieldType::ShortList => {
                Self::ShortList(ListBuilder::with_capacity(Int16Builder::new(), capacity))
            }
            FieldType::ShortFixedSizeList(size) => Self::ShortFixedSizeList(
                FixedSizeListBuilder::with_capacity(Int16Builder::new(), *size as i32, capacity),
            ),
            FieldType::Ushort => Self::Ushort(UInt16Builder::with_capacity(capacity)),
            FieldType::UshortList => {
                Self::UshortList(ListBuilder::with_capacity(UInt16Builder::new(), capacity))
            }
            FieldType::UshortFixedSizeList(size) => Self::UshortFixedSizeList(
                FixedSizeListBuilder::with_capacity(UInt16Builder::new(), *size as i32, capacity),
            ),
            FieldType::Int => Self::Int(Int32Builder::with_capacity(capacity)),
            FieldType::IntList => {
                Self::IntList(ListBuilder::with_capacity(Int32Builder::new(), capacity))
            }
            FieldType::IntFixedSizeList(size) => Self::IntFixedSizeList(
                FixedSizeListBuilder::with_capacity(Int32Builder::new(), *size as i32, capacity),
            ),
            FieldType::Uint => Self::Uint(UInt32Builder::with_capacity(capacity)),
            FieldType::UintList => {
                Self::UintList(ListBuilder::with_capacity(UInt32Builder::new(), capacity))
            }
            FieldType::UintFixedSizeList(size) => Self::UintFixedSizeList(
                FixedSizeListBuilder::with_capacity(UInt32Builder::new(), *size as i32, capacity),
            ),
            FieldType::Bigint => Self::Bigint(Int64Builder::with_capacity(capacity)),
            FieldType::BigintList => {
                Self::BigintList(ListBuilder::with_capacity(Int64Builder::new(), capacity))
            }
            FieldType::BigintFixedSizeList(size) => Self::BigintFixedSizeList(
                FixedSizeListBuilder::with_capacity(Int64Builder::new(), *size as i32, capacity),
            ),
            FieldType::Ubigint => Self::Ubigint(UInt64Builder::with_capacity(capacity)),
            FieldType::UbigintList => {
                Self::UbigintList(ListBuilder::with_capacity(UInt64Builder::new(), capacity))
            }
            FieldType::UbigintFixedSizeList(size) => Self::UbigintFixedSizeList(
                FixedSizeListBuilder::with_capacity(UInt64Builder::new(), *size as i32, capacity),
            ),
            FieldType::Float => Self::Float(Float32Builder::with_capacity(capacity)),
            FieldType::FloatList => {
                Self::FloatList(ListBuilder::with_capacity(Float32Builder::new(), capacity))
            }
            FieldType::FloatFixedSizeList(size) => Self::FloatFixedSizeList(
                FixedSizeListBuilder::with_capacity(Float32Builder::new(), *size as i32, capacity),
            ),
            FieldType::Double => Self::Double(Float64Builder::with_capacity(capacity)),
            FieldType::DoubleList => {
                Self::DoubleList(ListBuilder::with_capacity(Float64Builder::new(), capacity))
            }
            FieldType::DoubleFixedSizeList(size) => Self::DoubleFixedSizeList(
                FixedSizeListBuilder::with_capacity(Float64Builder::new(), *size as i32, capacity),
            ),
            FieldType::Char => {
                Self::Char(GenericStringBuilder::<i32>::with_capacity(capacity, 1024))
            }
            FieldType::String => {
                Self::String(GenericStringBuilder::<i32>::with_capacity(capacity, 1024))
            }
            FieldType::Lstring => Self::Lstring(LargeStringBuilder::with_capacity(capacity, 1024)),
            FieldType::Enum(items) => {
                let n = items.len();
                let values = arrow::array::StringArray::from(items.to_owned());
                let builder = StringDictionaryBuilder::<Int32Type>::new_with_dictionary(n, &values)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;
                Self::Enum(builder)
            }
            FieldType::Set(items) => {
                let n = items.len();
                let values = arrow::array::StringArray::from(items.to_owned());
                let item_builder =
                    StringDictionaryBuilder::<Int32Type>::new_with_dictionary(n, &values)
                        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e.to_string()))?;
                Self::Set(ListBuilder::<StringDictionaryBuilder<Int32Type>>::new(
                    item_builder,
                ))
            }
        };
        Ok(builder)
    }

    pub fn finish(&mut self) -> ArrayRef {
        match self {
            Self::Byte(builder) => Arc::new(builder.finish()),
            Self::ByteList(builder) => Arc::new(builder.finish()),
            Self::ByteFixedSizeList(builder) => Arc::new(builder.finish()),

            Self::Ubyte(builder) => Arc::new(builder.finish()),
            Self::UbyteList(builder) => Arc::new(builder.finish()),
            Self::UbyteFixedSizeList(builder) => Arc::new(builder.finish()),

            Self::Short(builder) => Arc::new(builder.finish()),
            Self::ShortList(builder) => Arc::new(builder.finish()),
            Self::ShortFixedSizeList(builder) => Arc::new(builder.finish()),

            Self::Ushort(builder) => Arc::new(builder.finish()),
            Self::UshortList(builder) => Arc::new(builder.finish()),
            Self::UshortFixedSizeList(builder) => Arc::new(builder.finish()),

            Self::Int(builder) => Arc::new(builder.finish()),
            Self::IntList(builder) => Arc::new(builder.finish()),
            Self::IntFixedSizeList(builder) => Arc::new(builder.finish()),

            Self::Uint(builder) => Arc::new(builder.finish()),
            Self::UintList(builder) => Arc::new(builder.finish()),
            Self::UintFixedSizeList(builder) => Arc::new(builder.finish()),

            Self::Bigint(builder) => Arc::new(builder.finish()),
            Self::BigintList(builder) => Arc::new(builder.finish()),
            Self::BigintFixedSizeList(builder) => Arc::new(builder.finish()),

            Self::Ubigint(builder) => Arc::new(builder.finish()),
            Self::UbigintList(builder) => Arc::new(builder.finish()),
            Self::UbigintFixedSizeList(builder) => Arc::new(builder.finish()),

            Self::Float(builder) => Arc::new(builder.finish()),
            Self::FloatList(builder) => Arc::new(builder.finish()),
            Self::FloatFixedSizeList(builder) => Arc::new(builder.finish()),

            Self::Double(builder) => Arc::new(builder.finish()),
            Self::DoubleList(builder) => Arc::new(builder.finish()),
            Self::DoubleFixedSizeList(builder) => Arc::new(builder.finish()),

            Self::Char(builder) => Arc::new(builder.finish()),
            Self::String(builder) => Arc::new(builder.finish()),
            Self::Lstring(builder) => Arc::new(builder.finish()),

            Self::Enum(builder) => Arc::new(builder.finish()),
            Self::Set(builder) => Arc::new(builder.finish()),
        }
    }
}

pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

/// Append a field value from a string to the column.
impl Push<&str> for FieldBuilder {
    fn push(&mut self, value: &str) -> io::Result<()> {
        match self {
            Self::Byte(builder) => {
                match value.parse::<i8>() {
                    Ok(val) => {
                        builder.append_value(val);
                    }
                    Err(_) => {
                        builder.append_null();
                    }
                };
            }
            Self::ByteList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<i8>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }
            Self::ByteFixedSizeList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<i8>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }

            Self::Ubyte(builder) => {
                match value.parse::<u8>() {
                    Ok(val) => {
                        builder.append_value(val);
                    }
                    Err(_) => {
                        builder.append_null();
                    }
                };
            }
            Self::UbyteList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<u8>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }
            Self::UbyteFixedSizeList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<u8>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }

            Self::Short(builder) => {
                match value.parse::<i16>() {
                    Ok(val) => {
                        builder.append_value(val);
                    }
                    Err(_) => {
                        builder.append_null();
                    }
                };
            }
            Self::ShortList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<i16>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }
            Self::ShortFixedSizeList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<i16>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }

            Self::Ushort(builder) => {
                match value.parse::<u16>() {
                    Ok(val) => {
                        builder.append_value(val);
                    }
                    Err(_) => {
                        builder.append_null();
                    }
                };
            }
            Self::UshortList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<u16>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }
            Self::UshortFixedSizeList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<u16>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }

            Self::Int(builder) => {
                match value.parse::<i32>() {
                    Ok(val) => {
                        builder.append_value(val);
                    }
                    Err(_) => {
                        builder.append_null();
                    }
                };
            }
            Self::IntList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<i32>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }
            Self::IntFixedSizeList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<i32>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }

            Self::Uint(builder) => {
                match value.parse::<u32>() {
                    Ok(val) => {
                        builder.append_value(val);
                    }
                    Err(_) => {
                        builder.append_null();
                    }
                };
            }
            Self::UintList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<u32>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }
            Self::UintFixedSizeList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<u32>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }

            Self::Bigint(builder) => {
                match value.parse::<i64>() {
                    Ok(val) => {
                        builder.append_value(val);
                    }
                    Err(_) => {
                        builder.append_null();
                    }
                };
            }
            Self::BigintList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<i64>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }
            Self::BigintFixedSizeList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<i64>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }

            Self::Ubigint(builder) => {
                match value.parse::<u64>() {
                    Ok(val) => {
                        builder.append_value(val);
                    }
                    Err(_) => {
                        builder.append_null();
                    }
                };
            }
            Self::UbigintList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<u64>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }
            Self::UbigintFixedSizeList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<u64>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }

            Self::Float(builder) => {
                match value.parse::<f32>() {
                    Ok(val) => {
                        builder.append_value(val);
                    }
                    Err(_) => {
                        builder.append_null();
                    }
                };
            }
            Self::FloatList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<f32>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }
            Self::FloatFixedSizeList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<f32>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }

            Self::Double(builder) => {
                match value.parse::<f64>() {
                    Ok(val) => {
                        builder.append_value(val);
                    }
                    Err(_) => {
                        builder.append_null();
                    }
                };
            }
            Self::DoubleList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<f64>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }
            Self::DoubleFixedSizeList(builder) => {
                let mut is_valid = false;
                let value = value.trim_end_matches(',');
                for elem in value.split(",") {
                    match elem.parse::<f64>() {
                        Ok(val) => {
                            builder.values().append_value(val);
                            is_valid = true;
                        }
                        Err(_) => {
                            builder.values().append_null();
                        }
                    };
                }
                builder.append(is_valid);
            }

            Self::Char(builder) => {
                if value.is_empty() {
                    builder.append_null();
                } else {
                    builder.append_value(value);
                }
            }
            Self::String(builder) => {
                if value.is_empty() {
                    builder.append_null();
                } else {
                    builder.append_value(value);
                }
            }
            Self::Lstring(builder) => {
                if value.is_empty() {
                    builder.append_null();
                } else {
                    builder.append_value(value);
                }
            }

            Self::Enum(builder) => {
                if value.is_empty() {
                    builder.append_null();
                } else {
                    builder.append_value(value);
                }
            }

            Self::Set(builder) => {
                if value.is_empty() {
                    builder.append(false);
                } else {
                    let value = value.trim_end_matches(',');
                    for elem in value.split(",") {
                        builder.values().append_value(elem);
                    }
                    builder.append(true);
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fielddef_new() {
        let field_def = FieldDef::new("foo".to_string(), FieldType::Int);
        assert_eq!(field_def.name, "foo");
        assert_eq!(field_def.ty, FieldType::Int);
    }

    #[test]
    fn test_fielddef_try_from_tuple() {
        let field_def = FieldDef::try_from(("foo".to_string(), "int".to_string()));
        assert!(field_def.is_ok());
        let field_def = field_def.unwrap();
        assert_eq!(field_def.name, "foo");
        assert_eq!(field_def.ty, FieldType::Int);
    }

    #[test]
    fn test_fielddef_try_from_tuple_invalid_type() {
        let field_def = FieldDef::try_from(("foo".to_string(), "invalid".to_string()));
        assert!(field_def.is_err());
    }

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

    #[test]
    fn test_field_arrow_type() {
        for ty in [
            FieldType::Byte,
            FieldType::ByteList,
            FieldType::ByteFixedSizeList(10),
            FieldType::Ubyte,
            FieldType::UbyteList,
            FieldType::UbyteFixedSizeList(10),
            FieldType::Short,
            FieldType::ShortList,
            FieldType::ShortFixedSizeList(10),
            FieldType::Ushort,
            FieldType::UshortList,
            FieldType::UshortFixedSizeList(10),
            FieldType::Int,
            FieldType::IntList,
            FieldType::IntFixedSizeList(10),
            FieldType::Uint,
            FieldType::UintList,
            FieldType::UintFixedSizeList(10),
            FieldType::Bigint,
            FieldType::BigintList,
            FieldType::BigintFixedSizeList(10),
            FieldType::Ubigint,
            FieldType::UbigintList,
            FieldType::UbigintFixedSizeList(10),
            FieldType::Float,
            FieldType::FloatList,
            FieldType::FloatFixedSizeList(10),
            FieldType::Double,
            FieldType::DoubleList,
            FieldType::DoubleFixedSizeList(10),
            FieldType::Char,
            FieldType::String,
            FieldType::Lstring,
            FieldType::Enum(vec!["item1".to_string(), "item2".to_string()]),
            FieldType::Set(vec!["item1".to_string(), "item2".to_string()]),
        ] {
            let mut builder = FieldBuilder::new(&ty, 1024).unwrap();
            let data_type = builder.finish().data_type().clone();
            assert_eq!(ty.arrow_type(), data_type);
        }
    }
}

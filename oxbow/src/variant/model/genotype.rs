use indexmap::IndexMap;
use std::io;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, BooleanBuilder, FixedSizeListBuilder, Float32Builder, GenericStringBuilder,
    Int32Builder, ListBuilder, StructArray, StructBuilder,
};
use arrow::datatypes::{DataType, Field as ArrowField};

use noodles::vcf::header::record::value::map::format::Number;
use noodles::vcf::header::record::value::map::format::Type;
use noodles::vcf::variant::record::samples::series::value::array::Array as Values;
use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
use noodles::vcf::variant::record::samples::series::Value;

pub const GT_ALLELE_TYPE: DataType = DataType::Int32;
pub const GT_PHASED_TYPE: DataType = DataType::Boolean;

/// A variant (VCF/BCF) sample-specific genotype field definition.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct GenotypeDef {
    pub name: String,
    pub ty: GenotypeType,
}

impl GenotypeDef {
    pub fn new(name: String, number: &Number, ty: &Type) -> Self {
        let ty = GenotypeType::from_noodles(number, ty);
        Self { name, ty }
    }

    pub fn try_from_strings(def: (String, String, String)) -> Result<Self, io::Error> {
        let (name, number, ty) = def;
        let number = match number.as_str() {
            "A" => Number::AlternateBases,
            "R" => Number::ReferenceAlternateBases,
            "G" => Number::Samples,
            "LA" => Number::LocalAlternateBases,
            "LR" => Number::LocalReferenceAlternateBases,
            "LG" => Number::LocalSamples,
            "P" => Number::Ploidy,
            "M" => Number::BaseModifications,
            "." => Number::Unknown,
            _ => {
                if let Ok(n) = number.parse::<usize>() {
                    Number::Count(n)
                } else {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "Invalid number parameter for FORMAT field '{}': {}",
                            name, number
                        ),
                    ));
                }
            }
        };
        let ty: Type = ty.parse().map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid type for FORMAT field '{}': {}", name, ty),
            )
        })?;
        Ok(Self::new(name, &number, &ty))
    }
}

/// A mapping of native sample-specific genotype field types to Arrow data types.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum GenotypeType {
    Character,
    CharacterList,
    CharacterFixedSizeList(usize),
    String,
    StringList,
    StringFixedSizeList(usize),
    Integer,
    IntegerList,
    IntegerFixedSizeList(usize),
    Float,
    FloatList,
    FloatFixedSizeList(usize),
}

impl GenotypeType {
    pub fn from_noodles(number: &Number, ty: &Type) -> Self {
        match number {
            Number::Count(1) => match ty {
                Type::Character => GenotypeType::Character,
                Type::String => GenotypeType::String,
                Type::Integer => GenotypeType::Integer,
                Type::Float => GenotypeType::Float,
            },
            Number::Count(n) => {
                let n = *n;
                match ty {
                    Type::Character => GenotypeType::CharacterFixedSizeList(n),
                    Type::String => GenotypeType::StringFixedSizeList(n),
                    Type::Integer => GenotypeType::IntegerFixedSizeList(n),
                    Type::Float => GenotypeType::FloatFixedSizeList(n),
                }
            }
            _ => {
                // "A", "R", "G", "LA", "LR", "LG", "P", "M", "."
                match ty {
                    Type::Character => GenotypeType::CharacterList,
                    Type::String => GenotypeType::StringList,
                    Type::Integer => GenotypeType::IntegerList,
                    Type::Float => GenotypeType::FloatList,
                }
            }
        }
    }

    pub fn number(&self) -> Option<usize> {
        match self {
            Self::Character => Some(1),
            Self::String => Some(1),
            Self::Integer => Some(1),
            Self::Float => Some(1),
            Self::CharacterFixedSizeList(n) => Some(*n),
            Self::StringFixedSizeList(n) => Some(*n),
            Self::IntegerFixedSizeList(n) => Some(*n),
            Self::FloatFixedSizeList(n) => Some(*n),
            Self::CharacterList => None,
            Self::StringList => None,
            Self::IntegerList => None,
            Self::FloatList => None,
        }
    }
}

/// A sample-specific genotype data builder.
pub enum GenotypeBuilder {
    Genotype(StructBuilder),
    Character(GenericStringBuilder<i32>),
    CharacterList(ListBuilder<GenericStringBuilder<i32>>),
    CharacterFixedSizeList(FixedSizeListBuilder<GenericStringBuilder<i32>>),
    String(GenericStringBuilder<i32>),
    StringList(ListBuilder<GenericStringBuilder<i32>>),
    StringFixedSizeList(FixedSizeListBuilder<GenericStringBuilder<i32>>),
    Integer(Int32Builder),
    IntegerList(ListBuilder<Int32Builder>),
    IntegerFixedSizeList(FixedSizeListBuilder<Int32Builder>),
    Float(Float32Builder),
    FloatList(ListBuilder<Float32Builder>),
    FloatFixedSizeList(FixedSizeListBuilder<Float32Builder>),
}

impl GenotypeBuilder {
    pub fn new(def: &GenotypeDef) -> Self {
        match def.ty {
            GenotypeType::Character => {
                GenotypeBuilder::Character(GenericStringBuilder::<i32>::new())
            }
            GenotypeType::CharacterList => {
                GenotypeBuilder::CharacterList(ListBuilder::new(GenericStringBuilder::<i32>::new()))
            }
            GenotypeType::CharacterFixedSizeList(n) => GenotypeBuilder::CharacterFixedSizeList(
                FixedSizeListBuilder::new(GenericStringBuilder::<i32>::new(), n as i32),
            ),
            GenotypeType::String => GenotypeBuilder::String(GenericStringBuilder::<i32>::new()),
            GenotypeType::StringList => {
                GenotypeBuilder::StringList(ListBuilder::new(GenericStringBuilder::<i32>::new()))
            }
            GenotypeType::StringFixedSizeList(n) => GenotypeBuilder::StringFixedSizeList(
                FixedSizeListBuilder::new(GenericStringBuilder::<i32>::new(), n as i32),
            ),
            GenotypeType::Integer => GenotypeBuilder::Integer(Int32Builder::new()),
            GenotypeType::IntegerList => {
                GenotypeBuilder::IntegerList(ListBuilder::new(Int32Builder::new()))
            }
            GenotypeType::IntegerFixedSizeList(n) => GenotypeBuilder::IntegerFixedSizeList(
                FixedSizeListBuilder::new(Int32Builder::new(), n as i32),
            ),
            GenotypeType::Float => GenotypeBuilder::Float(Float32Builder::new()),
            GenotypeType::FloatList => {
                GenotypeBuilder::FloatList(ListBuilder::new(Float32Builder::new()))
            }
            GenotypeType::FloatFixedSizeList(n) => GenotypeBuilder::FloatFixedSizeList(
                FixedSizeListBuilder::new(Float32Builder::new(), n as i32),
            ),
        }
    }

    pub fn new_gt() -> Self {
        let gt_allele_list_type =
            DataType::List(Arc::new(ArrowField::new("item", GT_ALLELE_TYPE, true)));
        let gt_phased_list_type =
            DataType::List(Arc::new(ArrowField::new("item", GT_PHASED_TYPE, true)));
        let gt_allele_builder = ListBuilder::new(Int32Builder::new());
        let gt_phased_builder = ListBuilder::new(BooleanBuilder::new());
        Self::Genotype(StructBuilder::new(
            vec![
                ArrowField::new("allele", gt_allele_list_type, true),
                ArrowField::new("phased", gt_phased_list_type, true),
            ],
            vec![Box::new(gt_allele_builder), Box::new(gt_phased_builder)],
        ))
    }

    pub fn arrow_type(&self) -> DataType {
        match self {
            Self::Genotype(_) => {
                let list_of_i32 =
                    DataType::List(Arc::new(ArrowField::new("item", DataType::Int32, true)));
                let list_of_bool =
                    DataType::List(Arc::new(ArrowField::new("item", DataType::Boolean, true)));
                let nested_fields = arrow::datatypes::Fields::from(vec![
                    ArrowField::new("allele", list_of_i32, true),
                    ArrowField::new("phased", list_of_bool, true),
                ]);
                DataType::Struct(nested_fields)
            }
            Self::Character(_) => DataType::Utf8,
            Self::CharacterList(_) => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                DataType::List(Arc::new(item))
            }
            Self::CharacterFixedSizeList(builder) => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                let n = builder.value_length();
                DataType::FixedSizeList(Arc::new(item), n)
            }
            Self::String(_) => DataType::Utf8,
            Self::StringList(_) => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                DataType::List(Arc::new(item))
            }
            Self::StringFixedSizeList(builder) => {
                let item = ArrowField::new("item", DataType::Utf8, true);
                let n = builder.value_length();
                DataType::FixedSizeList(Arc::new(item), n)
            }
            Self::Integer(_) => DataType::Int32,
            Self::IntegerList(_) => {
                let item = ArrowField::new("item", DataType::Int32, true);
                DataType::List(Arc::new(item))
            }
            Self::IntegerFixedSizeList(builder) => {
                let item = ArrowField::new("item", DataType::Int32, true);
                let n = builder.value_length();
                DataType::FixedSizeList(Arc::new(item), n)
            }
            Self::Float(_) => DataType::Float32,
            Self::FloatList(_) => {
                let item = ArrowField::new("item", DataType::Float32, true);
                DataType::List(Arc::new(item))
            }
            Self::FloatFixedSizeList(builder) => {
                let item = ArrowField::new("item", DataType::Float32, true);
                let n = builder.value_length();
                DataType::FixedSizeList(Arc::new(item), n)
            }
        }
    }

    pub fn get_arrow_field(&self, name: &str) -> ArrowField {
        ArrowField::new(name, self.arrow_type(), true)
    }

    pub fn append_null(&mut self) {
        match self {
            Self::Genotype(builder) => {
                builder
                    .field_builder::<ListBuilder<Int32Builder>>(0)
                    .unwrap()
                    .append(false);
                builder
                    .field_builder::<ListBuilder<BooleanBuilder>>(1)
                    .unwrap()
                    .append(false);
                builder.append(false);
            }
            Self::Character(builder) => builder.append_null(),
            Self::CharacterList(builder) => builder.append(false),
            Self::CharacterFixedSizeList(builder) => {
                for _ in 0..builder.value_length() {
                    builder.values().append_null();
                }
                builder.append(false);
            }
            Self::String(builder) => builder.append_null(),
            Self::StringList(builder) => builder.append(false),
            Self::StringFixedSizeList(builder) => {
                for _ in 0..builder.value_length() {
                    builder.values().append_null();
                }
                builder.append(false);
            }
            Self::Integer(builder) => builder.append_null(),
            Self::IntegerList(builder) => builder.append(false),
            Self::IntegerFixedSizeList(builder) => {
                for _ in 0..builder.value_length() {
                    builder.values().append_null();
                }
                builder.append(false);
            }
            Self::Float(builder) => builder.append_null(),
            Self::FloatList(builder) => builder.append(false),
            Self::FloatFixedSizeList(builder) => {
                for _ in 0..builder.value_length() {
                    builder.values().append_null();
                }
                builder.append(false);
            }
        };
    }

    pub fn append_value(&mut self, value: &Value) -> io::Result<()> {
        match value {
            Value::Genotype(gt) => {
                let (alleles, phasings): (Vec<Option<i32>>, Vec<bool>) = gt
                    .iter()
                    .filter_map(|result| match result {
                        Ok((allele, phasing)) => {
                            let allele = allele.map(|a| a as i32);
                            let phased = match phasing {
                                Phasing::Phased => true,
                                Phasing::Unphased => false,
                            };
                            Some((allele, phased))
                        }
                        Err(_) => None,
                    })
                    .unzip();
                match self {
                    GenotypeBuilder::Genotype(builder) => {
                        let allele_builder = builder
                            .field_builder::<ListBuilder<Int32Builder>>(0)
                            .unwrap();
                        for allele in alleles {
                            allele_builder.values().append_option(allele);
                        }
                        allele_builder.append(true);

                        let phased_builder = builder
                            .field_builder::<ListBuilder<BooleanBuilder>>(1)
                            .unwrap();
                        for phased in phasings {
                            phased_builder.values().append_value(phased);
                        }
                        phased_builder.append(true);

                        builder.append(true);
                    }
                    _ => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "Error parsing FORMAT field: type mismatch",
                        ));
                    }
                }
            }
            Value::Character(c) => match self {
                GenotypeBuilder::Character(builder) => {
                    builder.append_value(c.to_string());
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing FORMAT field: type mismatch",
                    ));
                }
            },
            Value::String(s) => match self {
                GenotypeBuilder::String(builder) => {
                    builder.append_value(s);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing FORMAT field: type mismatch",
                    ));
                }
            },
            Value::Integer(n) => match self {
                GenotypeBuilder::Integer(builder) => {
                    builder.append_value(*n);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing FORMAT field: type mismatch",
                    ));
                }
            },
            Value::Float(f) => match self {
                GenotypeBuilder::Float(builder) => {
                    builder.append_value(*f);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing FORMAT field: type mismatch",
                    ));
                }
            },
            Value::Array(array) => {
                self.append_values(array)?;
            }
        };
        Ok(())
    }

    fn append_values(&mut self, array: &Values) -> io::Result<()> {
        match array {
            Values::Character(values) => match self {
                GenotypeBuilder::CharacterList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(c) => {
                                non_null = true;
                                builder.values().append_value(c.to_string())
                            }
                            None => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                GenotypeBuilder::CharacterFixedSizeList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(c) => {
                                non_null = true;
                                builder.values().append_value(c.to_string())
                            }
                            None => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing FORMAT field: type mismatch",
                    ));
                }
            },
            Values::String(values) => match self {
                GenotypeBuilder::StringList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(s) => {
                                non_null = true;
                                builder.values().append_value(&s)
                            }
                            None => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                GenotypeBuilder::StringFixedSizeList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(s) => {
                                non_null = true;
                                builder.values().append_value(&s)
                            }
                            None => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing FORMAT field: type mismatch",
                    ));
                }
            },
            Values::Integer(values) => match self {
                GenotypeBuilder::IntegerList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(n) => {
                                non_null = true;
                                builder.values().append_value(n);
                            }
                            _ => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                GenotypeBuilder::IntegerFixedSizeList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(n) => {
                                non_null = true;
                                builder.values().append_value(n);
                            }
                            _ => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing FORMAT field: type mismatch",
                    ));
                }
            },
            Values::Float(values) => match self {
                GenotypeBuilder::FloatList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(f) => {
                                non_null = true;
                                builder.values().append_value(f);
                            }
                            _ => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                GenotypeBuilder::FloatFixedSizeList(builder) => {
                    let mut non_null = false;
                    for value in values.iter() {
                        match value? {
                            Some(f) => {
                                non_null = true;
                                builder.values().append_value(f);
                            }
                            _ => builder.values().append_null(),
                        }
                    }
                    builder.append(non_null);
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Error parsing FORMAT field: type mismatch",
                    ));
                }
            },
        }
        Ok(())
    }

    pub fn finish(&mut self) -> ArrayRef {
        match self {
            Self::Genotype(builder) => Arc::new(builder.finish()),
            Self::Character(builder) => Arc::new(builder.finish()),
            Self::CharacterList(builder) => Arc::new(builder.finish()),
            Self::CharacterFixedSizeList(builder) => Arc::new(builder.finish()),
            Self::String(builder) => Arc::new(builder.finish()),
            Self::StringList(builder) => Arc::new(builder.finish()),
            Self::StringFixedSizeList(builder) => Arc::new(builder.finish()),
            Self::Integer(builder) => Arc::new(builder.finish()),
            Self::IntegerList(builder) => Arc::new(builder.finish()),
            Self::IntegerFixedSizeList(builder) => Arc::new(builder.finish()),
            Self::Float(builder) => Arc::new(builder.finish()),
            Self::FloatList(builder) => Arc::new(builder.finish()),
            Self::FloatFixedSizeList(builder) => Arc::new(builder.finish()),
        }
    }
}

/// A single-sample genotype data builder.
///
/// Builds a StructArray, one variant per row, of a single sample's genotype
/// data, covering a defined list of genotype fields.
pub struct SampleStructBuilder {
    genotype_defs: Vec<GenotypeDef>,
    builders: IndexMap<GenotypeDef, GenotypeBuilder>,
}

impl SampleStructBuilder {
    pub fn new(genotype_defs: Vec<GenotypeDef>) -> Self {
        let mut builders = IndexMap::new();
        genotype_defs.iter().for_each(|def| {
            let builder = if def.name == "GT" {
                GenotypeBuilder::new_gt()
            } else {
                GenotypeBuilder::new(def)
            };
            builders.insert(def.clone(), builder);
        });

        Self {
            genotype_defs,
            builders,
        }
    }

    pub fn get_arrow_fields(&self) -> Vec<ArrowField> {
        self.builders
            .iter()
            .map(|(def, builder)| builder.get_arrow_field(def.name.as_str()))
            .collect()
    }

    pub fn push(&mut self, sample_data: IndexMap<String, Option<Value>>) -> io::Result<()> {
        for (def, builder) in self.builders.iter_mut() {
            let value = match sample_data.get(def.name.as_str()) {
                Some(Some(value)) => value,
                None | Some(None) => {
                    // Field is invalid or is not reported for this sample
                    builder.append_null();
                    continue;
                }
            };
            builder.append_value(value).map_err(|e| {
                io::Error::new(
                    e.kind(),
                    format!("Error processing field '{}': {}", def.name, e),
                )
            })?;
        }
        Ok(())
    }

    pub fn finish(&mut self) -> StructArray {
        let fields = self.get_arrow_fields().into_iter().map(Arc::new);
        let arrays: Vec<ArrayRef> = self
            .genotype_defs
            .iter()
            .map(|def| {
                let builder = self.builders.get_mut(def).unwrap();
                builder.finish()
            })
            .collect();
        StructArray::from(fields.zip(arrays).collect::<Vec<_>>())
    }
}

/// A single-field genotype data builder.
///
/// Builds a StructArray, one variant per row, of a single genotype field's
/// data series, covering a defined series of biological samples.
pub struct SeriesStructBuilder {
    sample_names: Vec<String>,
    builders: IndexMap<String, GenotypeBuilder>,
}

impl SeriesStructBuilder {
    pub fn new(genotype_def: GenotypeDef, sample_names: Vec<String>) -> Self {
        let mut builders = IndexMap::new();
        if genotype_def.name == "GT" {
            for sample_name in sample_names.iter() {
                let builder = GenotypeBuilder::new_gt();
                builders.insert(sample_name.clone(), builder);
            }
        } else {
            for sample_name in sample_names.iter() {
                let builder = GenotypeBuilder::new(&genotype_def);
                builders.insert(sample_name.clone(), builder);
            }
        }
        Self {
            sample_names,
            builders,
        }
    }

    pub fn get_arrow_fields(&self) -> Vec<ArrowField> {
        self.builders
            .iter()
            .map(|(sample_name, builder)| builder.get_arrow_field(sample_name.as_str()))
            .collect()
    }

    pub fn push(&mut self, series_data: IndexMap<String, Option<Value>>) -> io::Result<()> {
        for (sample_name, builder) in self.builders.iter_mut() {
            let value = match series_data.get(sample_name) {
                Some(Some(value)) => value,
                None | Some(None) => {
                    // Field is invalid or is not reported for this sample
                    builder.append_null();
                    continue;
                }
            };
            builder.append_value(value).map_err(|e| {
                io::Error::new(
                    e.kind(),
                    format!("Error processing field for sample '{}': {}", sample_name, e),
                )
            })?;
        }
        Ok(())
    }

    pub fn finish(&mut self) -> StructArray {
        let fields = self.get_arrow_fields().into_iter().map(Arc::new);
        let arrays: Vec<ArrayRef> = self
            .sample_names
            .iter()
            .map(|sample_name| {
                let builder = self.builders.get_mut(sample_name).unwrap();
                builder.finish()
            })
            .collect();
        StructArray::from(fields.zip(arrays).collect::<Vec<_>>())
    }
}

use std::str::FromStr;

use super::field_def::{FieldDef, FieldType};
use crate::OxbowError;

/// The 12 standard BED field names.
pub const STANDARD_FIELD_NAMES: [&str; 12] = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
    "blockCount",
    "blockSizes",
    "blockStarts",
];

/// Represents a BED schema defining the logical field structure.
///
/// A BED schema specifies `n` standard fields (3–12) and an optional set of
/// custom extension fields. Standard field names are fixed by the spec; their
/// Arrow types are determined by the format-specific Model (BED vs BBI).
/// Custom fields carry explicit types via [`FieldDef`].
///
/// # Examples
///
/// ```
/// use oxbow::bed::model::BedSchema;
///
/// let schema: BedSchema = "bed".parse().unwrap();
/// assert_eq!(schema.standard_field_count(), 6);
/// assert_eq!(schema.custom_field_count(), Some(0));
///
/// let schema: BedSchema = "bed12".parse().unwrap();
/// assert_eq!(schema.standard_field_count(), 12);
/// assert_eq!(schema.custom_field_count(), Some(0));
///
/// let schema: BedSchema = "bed6+3".parse().unwrap();
/// assert_eq!(schema.standard_field_count(), 6);
/// assert_eq!(schema.custom_field_count(), Some(3));
///
/// let schema: BedSchema = "bed6+".parse().unwrap();
/// assert_eq!(schema.standard_field_count(), 6);
/// assert_eq!(schema.custom_field_count(), None);
///
/// let schema: BedSchema = "bedgraph".parse().unwrap();
/// assert_eq!(schema.standard_field_count(), 3);
/// assert_eq!(schema.custom_field_count(), Some(1));
/// ```
///
/// # References
///
/// - [UCSC Genome Browser format documentation](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
/// - [GA4GH BED specification](https://samtools.github.io/hts-specs/BEDv1.pdf)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BedSchema {
    n: usize,
    m: Option<usize>,
    custom: Vec<FieldDef>,
}

impl BedSchema {
    /// Define a BED schema with `n` standard fields and an optional list of custom field definitions.
    pub fn new(n: usize, custom: Option<Vec<FieldDef>>) -> crate::Result<Self> {
        if n < 3 {
            return Err(OxbowError::invalid_input(format!(
                "Invalid BED schema: n < 3 (n={})",
                n
            )));
        } else if n > 12 {
            return Err(OxbowError::invalid_input(format!(
                "Invalid BED schema: n > 12 (n={})",
                n
            )));
        }

        let (m, custom) = match custom {
            Some(custom) => (Some(custom.len()), custom),
            // extension m is undefined: lump into "rest"
            None => (
                None,
                vec![FieldDef::new("rest".to_string(), FieldType::String)],
            ),
        };

        Ok(Self { n, m, custom })
    }

    /// Define a BED schema by specifying *n* and *m* parameters.
    ///
    /// `n` is the number of standard fields (3 to 12), and `m` is an optional number of custom
    /// fields. If `m` is 0, then the schema will only contain the standard fields. If `m` is `None`,
    /// then the schema will contain `n` standard fields followed by a single custom field named
    /// `rest` that contains the remainder of the data line. If `m` is a number greater than 0, then
    /// the schema will contain `n` standard fields followed by `m` custom fields named `BED{n}+1`,
    /// `BED{n}+2`, ..., `BED{n}+m`.
    ///
    /// All custom fields are assigned type `String`.
    pub fn new_from_nm(n: usize, m: Option<usize>) -> crate::Result<Self> {
        let custom_fields = m.map(|m| {
            (1..=m)
                .map(|i| FieldDef::new(format!("BED{}+{}", n, i), FieldType::String))
                .collect()
        });
        Self::new(n, custom_fields)
    }

    /// Define a BED schema for the bedGraph format.
    ///
    /// BedGraph can be considered a BED3+1 format, where the first three fields are the standard
    /// BED fields (chrom, start, end) and the fourth field is a floating point value representing
    /// a unique score assigned to the bases in the corresponding interval.
    pub fn new_bedgraph() -> crate::Result<Self> {
        Self::new(
            3,
            Some(vec![FieldDef::new("value".to_string(), FieldType::Float)]),
        )
    }

    /// All field names (standard + custom).
    pub fn field_names(&self) -> Vec<String> {
        let mut names: Vec<String> = STANDARD_FIELD_NAMES
            .iter()
            .take(self.n)
            .map(|&s| s.to_string())
            .collect();
        names.extend(self.custom.iter().map(|d| d.name.clone()));
        names
    }

    /// The number of standard fields (3–12).
    pub fn standard_field_count(&self) -> usize {
        self.n
    }

    /// The standard field names.
    pub fn standard_field_names(&self) -> Vec<String> {
        STANDARD_FIELD_NAMES
            .iter()
            .take(self.n)
            .map(|&s| s.to_string())
            .collect()
    }

    /// The number of custom fields, if defined.
    ///
    /// `None` means the extension count is undefined (BEDn+ mode).
    pub fn custom_field_count(&self) -> Option<usize> {
        self.m
    }

    /// The custom field names.
    pub fn custom_field_names(&self) -> Vec<String> {
        self.custom.iter().map(|d| d.name.clone()).collect()
    }

    /// The custom field definitions (with types).
    pub fn custom_fields(&self) -> &[FieldDef] {
        &self.custom
    }
}

impl std::fmt::Display for BedSchema {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(m) = self.m {
            if m == 0 {
                write!(f, "bed{}", self.n)
            } else {
                write!(f, "bed{}+{}", self.n, m)
            }
        } else {
            write!(f, "bed{}+", self.n)
        }
    }
}

impl FromStr for BedSchema {
    type Err = OxbowError;
    /// Define a BED schema using the shorthand BED*n+m* notation.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_ascii_lowercase();

        if s == "bed" {
            return Self::new_from_nm(6, Some(0));
        }

        fn parse_error(s: &str) -> OxbowError {
            OxbowError::invalid_input(format!("Invalid BED format specifier: {}", s))
        }

        if let Some(rest) = s.strip_prefix("bed") {
            if rest == "graph" {
                Self::new_bedgraph()
            } else if let Some(n) = rest.strip_suffix('+') {
                let n = n.parse::<usize>().map_err(|_| parse_error(&s))?;
                Self::new_from_nm(n, None)
            } else if let Some(pos) = rest.find('+') {
                let n = rest[..pos].parse::<usize>().map_err(|_| parse_error(&s))?;
                let m = rest[pos + 1..]
                    .parse::<usize>()
                    .map_err(|_| parse_error(&s))?;
                Self::new_from_nm(n, Some(m))
            } else {
                let n = rest.parse::<usize>().map_err(|_| parse_error(&s))?;
                Self::new_from_nm(n, Some(0))
            }
        } else {
            Err(parse_error(&s))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed_schema_bedn() {
        let spec: BedSchema = "bed".parse().unwrap();
        assert_eq!(spec.standard_field_count(), 6);
        assert_eq!(spec.custom_field_count(), Some(0));

        let spec: BedSchema = "bed6".parse().unwrap();
        assert_eq!(spec.standard_field_count(), 6);
        assert_eq!(spec.custom_field_count(), Some(0));

        let spec: BedSchema = "bed12".parse().unwrap();
        assert_eq!(spec.standard_field_count(), 12);
        assert_eq!(spec.custom_field_count(), Some(0));
    }

    #[test]
    fn test_bed_schema_bedn_plus_zero() {
        assert_eq!(
            "bed6".parse::<BedSchema>().unwrap(),
            "bed6+0".parse::<BedSchema>().unwrap()
        );
    }

    #[test]
    fn test_bed_schema_bedn_plus_m() {
        let spec: BedSchema = "bed6+3".parse().unwrap();
        assert_eq!(spec.standard_field_count(), 6);
        assert_eq!(spec.custom_field_count(), Some(3));

        let field_names = spec.field_names();
        assert_eq!(
            field_names,
            vec!["chrom", "start", "end", "name", "score", "strand", "BED6+1", "BED6+2", "BED6+3"]
        );
    }

    #[test]
    fn test_bed_schema_bedn_plus() {
        let spec: BedSchema = "bed6+".parse().unwrap();
        assert_eq!(spec.standard_field_count(), 6);
        assert_eq!(spec.custom_field_count(), None);

        let field_names = spec.field_names();
        assert_eq!(
            field_names,
            vec!["chrom", "start", "end", "name", "score", "strand", "rest"]
        );
    }

    #[test]
    fn test_bed_schema_bedgraph() {
        let spec: BedSchema = "bedgraph".parse().unwrap();
        assert_eq!(spec.standard_field_count(), 3);
        assert_eq!(spec.custom_field_count(), Some(1));
        assert_eq!(spec.custom_fields()[0].name, "value");
        assert_eq!(spec.custom_fields()[0].ty, FieldType::Float);
    }

    #[test]
    fn test_bed_schema_invalid() {
        let result = BedSchema::new(0, None);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Invalid BED schema: n < 3 (n=0)"
        );

        let result = BedSchema::new(13, None);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Invalid BED schema: n > 12 (n=13)"
        );

        let result: Result<BedSchema, _> = "invalid".parse();
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Invalid BED format specifier: invalid"
        );
    }

    #[test]
    fn test_bed_schema_display() {
        let spec: BedSchema = "bed6".parse().unwrap();
        assert_eq!(spec.to_string(), "bed6");

        let spec: BedSchema = "bed6+3".parse().unwrap();
        assert_eq!(spec.to_string(), "bed6+3");

        let spec: BedSchema = "bed6+".parse().unwrap();
        assert_eq!(spec.to_string(), "bed6+");

        let spec: BedSchema = "bed12".parse().unwrap();
        assert_eq!(spec.to_string(), "bed12");
    }
}

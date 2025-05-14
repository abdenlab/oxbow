use std::io;
use std::str::FromStr;

use super::field::{bed_standard_fields, FieldDef, FieldType};

/// Represents a typed BED schema using AutoSql-based field definitions.
///
/// A BED schema be created in several ways:
///
/// * Defined from `n` standard fields and a vec of custom [`FieldDef`].
/// * Defined from `n` and `Option<m>` parameters.
/// * Parsing a `BED`*n\[+\[m\]\]* or `bedGraph` specifier.
///
/// # Examples
///
/// ```
/// use oxbow::bbi::model::base::schema::BedSchema;
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
    fields: Vec<FieldDef>,
}

impl BedSchema {
    /// Define a BED schema with `n` standard fields and an optional list of custom field definitions.
    ///
    /// The `n` standard fields have pre-defined types. See the BED format specification for details.
    pub fn new(n: usize, custom: Option<Vec<FieldDef>>) -> io::Result<Self> {
        let bed_standard_fields = bed_standard_fields();
        if n < 3 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid BED schema: n < 3 (n={})", n),
            ));
        } else if n > 12 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid BED schema: n > 12 (n={})", n),
            ));
        }

        let mut fields: Vec<FieldDef> = bed_standard_fields
            .iter()
            .take(n)
            .map(FieldDef::try_from)
            .collect::<Result<Vec<_>, _>>()?;

        let m = match custom {
            // extension m is 0 or more
            Some(custom) => {
                let len = custom.len();
                if !custom.is_empty() {
                    fields.extend(custom);
                }
                Some(len)
            }
            // extension m is undefined: lump into "rest"
            None => {
                fields.push(FieldDef::new("rest".to_string(), FieldType::String));
                None
            }
        };

        Ok(Self { n, m, fields })
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
    /// All standard fields are assigned their default types. All custom fields are assigned type
    /// `String`.
    pub fn new_from_nm(n: usize, m: Option<usize>) -> io::Result<Self> {
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
    ///
    /// # Note
    /// Strictly speaking, a bedGraph file should contain no overlapping intervals, but this schema
    /// object does not enforce any constraint on the contents of the data besides the fields
    /// and their types.
    pub fn new_bedgraph() -> io::Result<Self> {
        Self::new(
            3,
            Some(vec![FieldDef::new("value".to_string(), FieldType::Float)]),
        )
    }

    pub fn fields(&self) -> &Vec<FieldDef> {
        &self.fields
    }

    pub fn field_names(&self) -> Vec<String> {
        self.fields.iter().map(|field| field.name.clone()).collect()
    }

    pub fn standard_field_count(&self) -> usize {
        self.n
    }

    pub fn standard_fields(&self) -> Vec<FieldDef> {
        self.fields.iter().take(self.n).cloned().collect()
    }

    pub fn custom_field_count(&self) -> Option<usize> {
        self.m
    }

    pub fn custom_fields(&self) -> Vec<FieldDef> {
        self.fields.iter().skip(self.n).cloned().collect()
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
    type Err = io::Error;
    /// Define a BED schema using the shorthand BED*n+m* notation.
    ///
    /// The specifier can be one of the following (case-insensitive):
    ///
    /// - `BED`: Equivalent to `BED6`.
    /// - `BED{n}`: `n` standard fields and 0 custom fields.
    /// - `BED{n}+{m}`: `n` standard fields followed by `m` custom fields.
    /// - `BED{n}+`: `n` standard fields followed by an undefined number of custom fields.
    /// - `bedGraph`: special case of a BED3+1, where the fourth field is a 32-bit floating point
    ///    field named `value`.
    ///
    /// # Notes
    /// - For `BED{n}+m`, custom fields are named `BED{n}+1`, `BED{n}+2`, ..., `BED{n}+m`.
    /// - For `BED{n}+`, custom fields are collapsed into a single field named `rest`.
    /// - Since, *n+m* notation does not specify the types of custom fields, they are assigned type
    ///   `String`.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_ascii_lowercase();

        if s == "bed" {
            // Interpret as BED6
            return Self::new_from_nm(6, Some(0));
        }

        fn parse_error(s: &str) -> io::Error {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid BED format specifier: {}", s),
            )
        }

        if let Some(rest) = s.strip_prefix("bed") {
            if rest == "graph" {
                Self::new_bedgraph()
            } else if let Some(n) = rest.strip_suffix('+') {
                // BEDn+
                let n = n.parse::<usize>().map_err(|_| parse_error(&s))?;
                Self::new_from_nm(n, None)
            } else if let Some(pos) = rest.find('+') {
                // BEDn+m
                let n = rest[..pos].parse::<usize>().map_err(|_| parse_error(&s))?;
                let m = rest[pos + 1..]
                    .parse::<usize>()
                    .map_err(|_| parse_error(&s))?;
                Self::new_from_nm(n, Some(m))
            } else {
                // BEDn
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

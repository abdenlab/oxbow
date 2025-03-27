use std::io;
use std::str::FromStr;

use super::field::DEFAULT_FIELD_NAMES;

/// Represents the shorthand *n+m* notation for a BED schema.
///
/// It can be one of the following (case-insensitive):
///
/// - `BED`: Equivalent to `BED6`.
/// - `BED{n}`: `n` standard fields and 0 custom fields.
/// - `BED{n}+{m}`: `n` standard fields followed by `m` custom fields.
/// - `BED{n}+`: `n` standard fields followed by an undefined number of custom fields.
///
/// While the 12 standard fields have defined types, this notation does not specify the types of
/// custom fields.
///
/// For `BED{n}+m`, custom fields are named `BED{n}+1`, `BED{n}+2`, ..., `BED{n}+m`.
///
/// For `BED{n}+`, custom fields are collapsed into a single field named `rest`.
///
/// # Examples
///
/// ```
/// use std::str::FromStr;
/// use crate::bed::model::BedSchema;
///
/// let spec: BedSchema = "bed".parse()?;
/// assert_eq!(spec.standard_field_count(), 6);
/// assert_eq!(spec.custom_field_count(), 0);
///
/// let spec: BedSchema = "bed12".parse()?;
/// assert_eq!(spec.standard_field_count(), 12);
/// assert_eq!(spec.custom_field_count(), 0);
///
/// let spec: BedSchema = "bed6+3".parse()?;
/// assert_eq!(spec.standard_field_count(), 6);
/// assert_eq!(spec.custom_field_count(), Some(3));
///
/// let spec: BedSchema = "bed6+".parse()?;
/// assert_eq!(spec.standard_field_count(), 6);
/// assert_eq!(spec.custom_field_count(), None);
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
}

impl BedSchema {
    pub fn new(n: usize, m: Option<usize>) -> io::Result<Self> {
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
        Ok(Self { n, m })
    }

    pub fn standard_field_count(&self) -> usize {
        self.n
    }

    pub fn custom_field_count(&self) -> Option<usize> {
        self.m
    }

    pub fn field_names(&self) -> Vec<String> {
        let mut names = DEFAULT_FIELD_NAMES
            .iter()
            .take(self.n)
            .map(|name| name.to_string())
            .collect::<Vec<String>>();
        if let Some(m) = self.m {
            for i in 1..=m {
                names.push(format!("BED{}+{}", self.n, i));
            }
        } else {
            names.push("rest".to_string());
        }
        names
    }
}

impl FromStr for BedSchema {
    type Err = io::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.to_ascii_lowercase();

        if s == "bed" {
            return Self::new(6, Some(0));
        }

        fn parse_error(s: &str) -> io::Error {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Invalid BED format specifier: {}", s),
            )
        }

        if s.starts_with("bed") {
            let rest = &s[3..];
            if rest.ends_with('+') {
                // BEDn+
                let n = rest[..rest.len() - 1]
                    .parse::<usize>()
                    .map_err(|_| parse_error(&s))?;
                Self::new(n, None)
            } else if let Some(pos) = rest.find('+') {
                // BEDn+m
                let n = rest[..pos].parse::<usize>().map_err(|_| parse_error(&s))?;
                let m = rest[pos + 1..]
                    .parse::<usize>()
                    .map_err(|_| parse_error(&s))?;
                Self::new(n, Some(m))
            } else {
                // BEDn
                let n = rest.parse::<usize>().map_err(|_| parse_error(&s))?;
                Self::new(n, Some(0))
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
}

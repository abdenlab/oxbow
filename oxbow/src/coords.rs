//! Genomic coordinate systems and region types.

use crate::{OxbowError, Result};

/// Genomic coordinate system.
///
/// The notation `XY` encodes the base of the start coordinate (`X`) and the
/// base of the end coordinate (`Y`):
///
/// - `"11"` — 1-based start, 1-based end (closed; SAM/VCF/GFF convention)
/// - `"01"` — 0-based start, 1-based end (half-open; BED/BBI convention)
///
/// End coordinates are numerically identical in both systems; only start
/// positions differ. Use [`CoordSystem::start_offset_from`] to get the
/// additive offset needed to convert a start value from one system to another.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CoordSystem {
    /// 1-based start, closed end.
    OneClosed,
    /// 0-based start, half-open end.
    ZeroHalfOpen,
}

impl CoordSystem {
    /// Returns the additive offset to apply to a start coordinate when
    /// converting from `source_cs` to `self`.
    ///
    /// - `OneClosed` → `ZeroHalfOpen`: `-1`
    /// - `ZeroHalfOpen` → `OneClosed`: `+1`
    /// - same → same: `0`
    pub fn start_offset_from(self, source_cs: CoordSystem) -> i32 {
        match (source_cs, self) {
            (CoordSystem::OneClosed, CoordSystem::ZeroHalfOpen) => -1,
            (CoordSystem::ZeroHalfOpen, CoordSystem::OneClosed) => 1,
            _ => 0,
        }
    }
}

impl std::fmt::Display for CoordSystem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CoordSystem::OneClosed => write!(f, "11"),
            CoordSystem::ZeroHalfOpen => write!(f, "01"),
        }
    }
}

impl std::str::FromStr for CoordSystem {
    type Err = OxbowError;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "11" => Ok(CoordSystem::OneClosed),
            "01" => Ok(CoordSystem::ZeroHalfOpen),
            other => Err(OxbowError::invalid_input(format!(
                "invalid coordinate system '{other}'; expected \"01\" or \"11\""
            ))),
        }
    }
}

/// A genomic region.
///
/// Represents a query region on a named reference sequence. Internally,
/// coordinates are always stored as **0-based half-open** `[start, end)`.
///
/// # Parsing
///
/// Regions can be parsed from strings in two styles:
///
/// **UCSC notation** — coordinate system is ambiguous and must be supplied:
/// ```text
/// "chr1"                   → whole chromosome (start=0, end=None)
/// "chr1:10000-20000"       → depends on coord_system
/// "chr1:10,000-20,000"     → same, with comma separators
/// "chr1:10_000-20_000"     → same, with underscore separators
/// ```
///
/// **Explicit bracket notation** — self-describing coordinate system:
/// ```text
/// "chr1:[10000,20000)"     → 0-based half-open
/// "chr1:[10001,20000]"     → 1-based closed (normalized to 0-based internally)
/// ```
///
/// Use [`Region::parse`] with a default [`CoordSystem`] for UCSC notation,
/// or [`Region::from_str`] which assumes 1-based closed for bare UCSC strings.
///
/// # Examples
///
/// ```
/// use oxbow::{CoordSystem, Region};
///
/// // Construct directly (0-based half-open)
/// let r = Region::new("chr1", Some(10000), Some(20000));
/// assert_eq!(r.start, 10000);
/// assert_eq!(r.end, Some(20000));
///
/// let r = Region::new("chr1", None, None);
/// assert_eq!(r.start, 0);
/// assert_eq!(r.end, None);
///
/// // UCSC notation with explicit coord system
/// let r = Region::parse("chr1:10,001-20,000", CoordSystem::OneClosed).unwrap();
/// assert_eq!(r.start, 10000);  // normalized to 0-based
/// assert_eq!(r.end, Some(20000));
///
/// // Explicit bracket notation (coord system in the string itself)
/// let r: Region = "chr1:[10000,20000)".parse().unwrap();
/// assert_eq!(r.start, 10000);
/// assert_eq!(r.end, Some(20000));
///
/// let r: Region = "chr1:[10001,20000]".parse().unwrap();
/// assert_eq!(r.start, 10000);  // 1-based 10001 → 0-based 10000
/// assert_eq!(r.end, Some(20000));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    /// Reference sequence name.
    pub name: String,
    /// 0-based start position (inclusive).
    pub start: u64,
    /// 0-based end position (exclusive). `None` means to the end of the
    /// reference sequence.
    pub end: Option<u64>,
}

impl Region {
    /// Create a new region with 0-based half-open coordinates.
    ///
    /// `start` defaults to 0 if `None`.
    pub fn new(name: impl Into<String>, start: Option<u64>, end: Option<u64>) -> Self {
        Self {
            name: name.into(),
            start: start.unwrap_or(0),
            end,
        }
    }

    /// Parse a region string using the given coordinate system for UCSC
    /// notation. Explicit bracket notation overrides `coord_system`.
    pub fn parse(s: &str, coord_system: CoordSystem) -> Result<Self> {
        // Try explicit bracket notation first.
        if let Some(result) = Self::try_parse_bracket(s) {
            return result;
        }
        // Fall back to UCSC notation with the provided coord system.
        Self::parse_ucsc(s, coord_system)
    }

    /// Parse UCSC-style `name[:start[-end]]`.
    fn parse_ucsc(s: &str, coord_system: CoordSystem) -> Result<Self> {
        if s.is_empty() {
            return Err(OxbowError::invalid_input("empty region string"));
        }

        let (name, interval) = match s.rsplit_once(':') {
            Some((name, "")) => (name, None),
            Some((name, suffix)) => (name, Some(suffix)),
            None => (s, None),
        };

        if name.is_empty() {
            return Err(OxbowError::invalid_input("empty reference name"));
        }

        let (start, end) = match interval {
            None => (None, None),
            Some(iv) => {
                let parts: Vec<&str> = iv.splitn(2, '-').collect();
                let start = parse_number(parts[0])?;
                let end = if parts.len() == 2 {
                    Some(parse_number(parts[1])?)
                } else {
                    None
                };
                (Some(start), end)
            }
        };

        // Normalize to 0-based half-open.
        let (start, end) = match coord_system {
            CoordSystem::OneClosed => {
                // 1-based closed → 0-based half-open: start -= 1, end unchanged
                (start.map(|s| s.saturating_sub(1)), end)
            }
            CoordSystem::ZeroHalfOpen => (start, end),
        };

        Ok(Self::new(name, start, end))
    }

    /// Try to parse explicit bracket notation `name:[start,end)` or
    /// `name:[start,end]`. Returns `None` if the string doesn't match
    /// bracket notation.
    fn try_parse_bracket(s: &str) -> Option<Result<Self>> {
        let (name, rest) = s.rsplit_once(':')?;
        if !rest.starts_with('[') {
            return None;
        }

        let result = (|| {
            let rest = &rest[1..]; // strip leading '['

            let (half_open, body) = if let Some(body) = rest.strip_suffix(')') {
                (true, body)
            } else if let Some(body) = rest.strip_suffix(']') {
                (false, body)
            } else {
                return Err(OxbowError::invalid_input(format!(
                    "bracket notation must end with ')' or ']': '{s}'"
                )));
            };

            // Strip underscores first (commas are ambiguous in bracket notation).
            let body: String = body.chars().filter(|c| *c != '_').collect();
            let (start_str, end_str) = body.split_once(',').ok_or_else(|| {
                OxbowError::invalid_input(format!("bracket notation requires 'start,end': '{s}'"))
            })?;

            let start = start_str.parse::<u64>().map_err(|_| {
                OxbowError::invalid_input(format!("invalid start in bracket notation: '{s}'"))
            })?;
            let end = end_str.parse::<u64>().map_err(|_| {
                OxbowError::invalid_input(format!("invalid end in bracket notation: '{s}'"))
            })?;

            // Normalize to 0-based half-open.
            let (start, end) = if half_open {
                // [start, end) — 0-based half-open, already normalized
                (start, end)
            } else {
                // [start, end] — 1-based closed
                // start: 1-based → 0-based = start - 1
                // end: closed → half-open = end (numerically same)
                (start.saturating_sub(1), end)
            };

            Ok(Self::new(name, Some(start), Some(end)))
        })();

        Some(result)
    }

    /// Convert to a noodles `Region` for index-based seeking.
    ///
    /// Noodles regions are 1-based with inclusive bounds.
    pub fn to_noodles(&self) -> std::result::Result<noodles::core::Region, OxbowError> {
        use noodles::core::Position;

        match (self.start, self.end) {
            (0, None) => Ok(noodles::core::Region::new(self.name.as_str(), ..)),
            (s, None) => {
                let start = Position::try_from(s as usize + 1)
                    .map_err(|_| OxbowError::invalid_input("start position out of range"))?;
                Ok(noodles::core::Region::new(self.name.as_str(), start..))
            }
            (s, Some(e)) => {
                let start = Position::try_from(s as usize + 1)
                    .map_err(|_| OxbowError::invalid_input("start position out of range"))?;
                let end = Position::try_from(e as usize)
                    .map_err(|_| OxbowError::invalid_input("end position out of range"))?;
                Ok(noodles::core::Region::new(self.name.as_str(), start..=end))
            }
        }
    }
}

impl std::str::FromStr for Region {
    type Err = OxbowError;

    /// Parse a region string. Bracket notation is self-describing; bare
    /// UCSC notation assumes 1-based closed (the most common convention).
    fn from_str(s: &str) -> Result<Self> {
        Self::parse(s, CoordSystem::OneClosed)
    }
}

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.name)?;
        match (self.start, self.end) {
            (0, None) => {}
            (s, None) => write!(f, ":[{s},)")?,
            (s, Some(e)) => write!(f, ":[{s},{e})")?,
        }
        Ok(())
    }
}

/// Parse a number string, stripping `,` and `_` thousands separators.
fn parse_number(s: &str) -> Result<u64> {
    let cleaned: String = s.chars().filter(|c| *c != ',' && *c != '_').collect();
    cleaned
        .parse::<u64>()
        .map_err(|_| OxbowError::invalid_input(format!("invalid number: '{s}'")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let r = Region::new("chr1", Some(100), Some(200));
        assert_eq!(r.name, "chr1");
        assert_eq!(r.start, 100);
        assert_eq!(r.end, Some(200));
    }

    #[test]
    fn test_new_defaults() {
        let r = Region::new("chr1", None, None);
        assert_eq!(r.start, 0);
        assert_eq!(r.end, None);
    }

    #[test]
    fn test_ucsc_one_closed() {
        let r = Region::parse("chr1:10001-20000", CoordSystem::OneClosed).unwrap();
        assert_eq!(r.name, "chr1");
        assert_eq!(r.start, 10000);
        assert_eq!(r.end, Some(20000));
    }

    #[test]
    fn test_ucsc_zero_half_open() {
        let r = Region::parse("chr1:10000-20000", CoordSystem::ZeroHalfOpen).unwrap();
        assert_eq!(r.name, "chr1");
        assert_eq!(r.start, 10000);
        assert_eq!(r.end, Some(20000));
    }

    #[test]
    fn test_ucsc_whole_chrom() {
        let r = Region::parse("chr1", CoordSystem::OneClosed).unwrap();
        assert_eq!(r.name, "chr1");
        assert_eq!(r.start, 0);
        assert_eq!(r.end, None);
    }

    #[test]
    fn test_ucsc_start_only() {
        let r = Region::parse("chr1:5000", CoordSystem::OneClosed).unwrap();
        assert_eq!(r.start, 4999);
        assert_eq!(r.end, None);
    }

    #[test]
    fn test_ucsc_thousands_separators() {
        let r = Region::parse("chr1:10,001-20,000", CoordSystem::OneClosed).unwrap();
        assert_eq!(r.start, 10000);
        assert_eq!(r.end, Some(20000));

        let r = Region::parse("chr1:10_001-20_000", CoordSystem::OneClosed).unwrap();
        assert_eq!(r.start, 10000);
        assert_eq!(r.end, Some(20000));
    }

    #[test]
    fn test_bracket_half_open() {
        let r: Region = "chr1:[10000,20000)".parse().unwrap();
        assert_eq!(r.start, 10000);
        assert_eq!(r.end, Some(20000));
    }

    #[test]
    fn test_bracket_closed() {
        let r: Region = "chr1:[10001,20000]".parse().unwrap();
        assert_eq!(r.start, 10000);
        assert_eq!(r.end, Some(20000));
    }

    #[test]
    fn test_bracket_overrides_coord_system() {
        // Even with ZeroHalfOpen context, bracket notation is self-describing
        let r = Region::parse("chr1:[10001,20000]", CoordSystem::ZeroHalfOpen).unwrap();
        assert_eq!(r.start, 10000); // interpreted as 1-based closed
    }

    #[test]
    fn test_bracket_with_separators() {
        let r: Region = "chr1:[10_000,20_000)".parse().unwrap();
        assert_eq!(r.start, 10000);
        assert_eq!(r.end, Some(20000));
    }

    #[test]
    fn test_display_roundtrip() {
        let r = Region::new("chr1", Some(10000), Some(20000));
        assert_eq!(r.to_string(), "chr1:[10000,20000)");

        let parsed: Region = r.to_string().parse().unwrap();
        assert_eq!(r, parsed);
    }

    #[test]
    fn test_display_whole_chrom() {
        let r = Region::new("chr1", None, None);
        assert_eq!(r.to_string(), "chr1");
    }

    #[test]
    fn test_to_noodles_full_range() {
        let r = Region::new("chr1", Some(10000), Some(20000));
        let nr = r.to_noodles().unwrap();
        assert_eq!(nr.name(), &b"chr1"[..]);
        // 0-based 10000 → 1-based 10001
        let start = noodles::core::Position::try_from(10001).unwrap();
        let end = noodles::core::Position::try_from(20000).unwrap();
        assert_eq!(nr.start(), std::ops::Bound::Included(start));
        assert_eq!(nr.end(), std::ops::Bound::Included(end));
    }

    #[test]
    fn test_to_noodles_whole_chrom() {
        let r = Region::new("chr1", None, None);
        let nr = r.to_noodles().unwrap();
        assert_eq!(nr.start(), std::ops::Bound::Unbounded);
        assert_eq!(nr.end(), std::ops::Bound::Unbounded);
    }

    #[test]
    fn test_empty_string_errors() {
        assert!(Region::parse("", CoordSystem::OneClosed).is_err());
    }

    #[test]
    fn test_invalid_bracket_notation() {
        assert!("chr1:[10000,20000".parse::<Region>().is_err());
        assert!("chr1:[10000)".parse::<Region>().is_err());
    }
}

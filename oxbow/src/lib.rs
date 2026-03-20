//! # oxbow
//!
//! **`oxbow`** reads genomic data formats 🧬 as Apache Arrow 🏹.
//!
//! With the oxbow Rust library, you can serialize native formats into [Arrow IPC](https://arrow.apache.org/docs/python/ipc.html)
//! , stream larger-than-memory files as Arrow [RecordBatches](https://docs.rs/arrow/latest/arrow/record_batch/struct.RecordBatch.html)
//! with zero-copy over FFI, and more!
//!
//! ⚠️ The Rust API is under active development and is not yet stable. The API may change in future releases.
//!
//! [Source on GitHub](https://github.com/abdenlab/oxbow).
//!
//!
//! ## Features
//!
//! - 🚀 Supports commonly used file formats from the [htslib/GA4GH](https://www.htslib.org/) and the
//!   [UCSC](https://genome.ucsc.edu/FAQ/FAQformat.html) ecosystems.
//! - 🔍 Support for compression, indexing, column projection, and genomic range querying.
//! - 🔧 Support for nested fields and complex, typed schemas (e.g., SAM tags,
//!   VCF `INFO` and `FORMAT` fields, AutoSql, etc.).
//!
//!
//! ## Scanners
//!
//! The main interface to read files are the scanners. Each scanner is a parser for a specific
//! format and provides scanning methods that return an iterator implementing the
//! [`arrow::record_batch::RecordBatchReader`] trait.
//!
//! ### Sequence formats
//!
//! - [`fasta`](crate::sequence::FastaScanner): Scan FASTA files as Arrow RecordBatches.
//! - [`fastq`](crate::sequence::FastqScanner): Scan FASTQ files as Arrow RecordBatches.
//!
//! ### Alignment formats
//! - [`sam`](crate::alignment::SamScanner): Scan SAM files as Arrow RecordBatches.
//! - [`bam`](crate::alignment::BamScanner): Scan BAM files as Arrow RecordBatches.
//! - [`cram`](crate::alignment::CramScanner): Scan CRAM files as Arrow RecordBatches.
//!
//! ### Variant formats
//! - [`vcf`](crate::variant::VcfScanner): Scan VCF files as Arrow RecordBatches.
//! - [`bcf`](crate::variant::BcfScanner): Scan BCF files as Arrow RecordBatches.
//!
//! ### Interval feature formats
//! - [`bed`](crate::bed::BedScanner): Scan BED files as Arrow RecordBatches.
//! - [`gtf`](crate::gxf::GtfScanner): Scan GXF files as Arrow RecordBatches.
//! - [`gff`](crate::gxf::GffScanner): Scan GFF files as Arrow RecordBatches.
//!
//! ### UCSC Big Binary Indexed (BBI) formats
//! - [`bigbed`](crate::bbi::BigBedScanner): Scan BigBed files as Arrow RecordBatches.
//! - [`bigwig`](crate::bbi::BigWigScanner): Scan BigWig files as Arrow RecordBatches.
//! - [`BBI zoom`](crate::bbi::BBIZoomScanner): Scan zoom level summary statistics from
//!   BigWig/BigBed as Arrow RecordBatches.
//!
//!
//! ## License
//!
//! Licensed under MIT or Apache-2.0.
//!
/// The version of the oxbow core library.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub mod alignment;
pub mod batch;
pub mod bbi;
pub mod bed;
pub mod error;
pub mod gxf;
pub mod sequence;
pub mod util;
pub mod variant;

pub use error::{OxbowError, Result};

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

#[derive(Debug, Clone)]
pub enum Select<T> {
    /// Select specific items explicitly
    Some(Vec<T>),
    /// Omit (explicitly empty)
    Omit,
    /// Select all items (wildcard)
    All,
}

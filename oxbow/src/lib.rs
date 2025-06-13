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
pub mod alignment;
pub mod bbi;
pub mod bed;
pub mod gxf;
pub mod sequence;
pub mod util;
pub mod variant;

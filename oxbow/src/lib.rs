//! **oxbow** reads specialized bioinformatic file formats as Apache Arrow IPC.
//!
//! # Examples
//!
//! ## Read all records in a BAM File.
//!
//! ```no_run
//! use oxbow::bam::BamReader;
//!
//! let mut reader = BamReader::new("sample.bam").unwrap();
//! let ipc = reader.records_to_ipc(None).unwrap();
//! ```
//!
//! ## Query records
//!
//! Querying allows filtering records by region. It requires an associated BAM index (BAI).
//!
//! ```no_run
//! use oxbow::bam::BamReader;
//!
//! let mut reader = BamReader::new("sample.bam").unwrap();
//! let ipc = reader.records_to_ipc(Some("chr1:1-100000")).unwrap();
//! ```
//!

mod batch_builder;
pub mod fasta;
pub mod fastq;
pub mod vpos;
pub mod bam;
// pub mod cram;
pub mod vcf;
pub mod bcf;

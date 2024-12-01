use extendr_api::prelude::*;
use oxbow::bam::BamReader;
use oxbow::fasta::FastaReader;
use oxbow::fastq::FastqReader;
// use oxbow::cram::CramReader;
use oxbow::bcf::BcfReader;
use oxbow::vcf::VcfReader;
use oxbow::vpos;

/// Return Arrow IPC format from a FASTA file.
/// @export
#[extendr]
fn read_fasta(path: &str, region: Option<&str>) -> Vec<u8> {
    let mut reader = FastaReader::new(path).unwrap();
    reader.records_to_ipc(region).unwrap()
}

/// Return Arrow IPC format from a FASTQ file.
/// @export
#[extendr]
fn read_fastq(path: &str) -> Vec<u8> {
    let mut reader = FastqReader::new_from_path(path).unwrap();
    reader.records_to_ipc().unwrap()
}

/// Return Arrow IPC format from a BAM file.
/// @export
#[extendr]
fn read_bam(path: &str, region: Option<&str>) -> Vec<u8> {
    let mut reader = BamReader::new_from_path(path).unwrap();
    reader.records_to_ipc(region).unwrap()
}

/// Return Arrow IPC format from a BAM file.
/// @export
#[extendr]
fn read_bam_vpos(path: &str, cpos_lo: u64, upos_lo: u16, cpos_hi: u64, upos_hi: u16) -> Vec<u8> {
    let mut reader = BamReader::new_from_path(path).unwrap();
    reader
        .records_to_ipc_from_vpos((cpos_lo, upos_lo), (cpos_hi, upos_hi))
        .unwrap()
}

/// Return Arrow IPC format from a VCF file.
/// @export
#[extendr]
fn read_vcf(path: &str, region: Option<&str>) -> Vec<u8> {
    let mut reader = VcfReader::new_from_path(path).unwrap();
    reader.records_to_ipc(region).unwrap()
}

/// Return Arrow IPC format from a VCF file.
/// @export
#[extendr]
fn read_vcf_vpos(path: &str, cpos_lo: u64, upos_lo: u16, cpos_hi: u64, upos_hi: u16) -> Vec<u8> {
    let mut reader = VcfReader::new_from_path(path).unwrap();
    reader
        .records_to_ipc_from_vpos((cpos_lo, upos_lo), (cpos_hi, upos_hi))
        .unwrap()
}

/// Return Arrow IPC format from a BCF file.
/// @export
#[extendr]
fn read_bcf(path: &str, region: Option<&str>) -> Vec<u8> {
    let mut reader = BcfReader::new_from_path(path).unwrap();
    reader.records_to_ipc(region).unwrap()
}

/// Return Arrow IPC format from a BCF file.
/// @export
#[extendr]
fn read_bcf_vpos(path: &str, cpos_lo: u64, upos_lo: u16, cpos_hi: u64, upos_hi: u16) -> Vec<u8> {
    let mut reader = BcfReader::new_from_path(path).unwrap();
    reader
        .records_to_ipc_from_vpos((cpos_lo, upos_lo), (cpos_hi, upos_hi))
        .unwrap()
}

/// Return a virtual position partition with an approximate uncompressed spacing.
/// @export
#[extendr]
fn partition_from_index_file(path: &str, chunksize: u64) -> List {
    let pos = vpos::partition_from_index_file(path, chunksize);
    list!(
        compressed_offset = pos.iter().map(|x| x.0).collect::<Vec<_>>(),
        bin_index = pos.iter().map(|x| x.1).collect::<Vec<_>>()
    )
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod oxbow;
    fn read_fasta;
    fn read_fastq;
    fn read_bam;
    fn read_bam_vpos;
    fn read_vcf;
    fn read_vcf_vpos;
    fn read_bcf;
    fn read_bcf_vpos;
}

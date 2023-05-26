use extendr_api::prelude::*;
use oxbow::bam::BamReader;
use oxbow::vcf::VcfReader;

/// Return Arrow IPC format from a BAM file.
/// @export
#[extendr]
fn read_bam(path: &str, region: Option<&str>) -> Vec<u8> {
    let mut reader = BamReader::new(path).unwrap();
    reader.records_to_ipc(region).unwrap()
}

/// Return Arrow IPC format from a VCF file.
/// @export
#[extendr]
fn read_vcf(path: &str, region: Option<&str>) -> Vec<u8> {
    let mut reader = VcfReader::new(path).unwrap();
    reader.records_to_ipc(region).unwrap()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod oxbow;
    fn read_bam;
    fn read_vcf;
}

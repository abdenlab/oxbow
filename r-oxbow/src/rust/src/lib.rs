use std::io::BufReader;

use extendr_api::prelude::*;

use flate2::bufread::MultiGzDecoder;
use noodles::bgzf::IndexedReader as IndexedBgzfReader;
use noodles::core::Region;

use oxbow::sequence::{FastaScanner, FastqScanner};
use oxbow::util::batches_to_ipc;

pub const BUFFER_SIZE_BYTES: usize = const { 1024 * 1024 };

/// Return Arrow IPC format from a FASTQ file.
/// @export
#[extendr]
fn read_fastq(path: &str, fields: Option<Vec<String>>) -> Vec<u8> {
    let compressed = path.ends_with(".gz");
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();
    let scanner = FastqScanner::new();

    let ipc = if compressed {
        let gz_reader = std::io::BufReader::new(MultiGzDecoder::new(reader));
        let fmt_reader = noodles::fastq::io::Reader::new(gz_reader);
        let batches = scanner.scan(fmt_reader, fields, None, None).unwrap();
        batches_to_ipc(batches)
    } else {
        let fmt_reader = noodles::fastq::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, fields, None, None).unwrap();
        batches_to_ipc(batches)
    };

    ipc.unwrap()
}

/// Return Arrow IPC format from a FASTA file.
/// @export
#[extendr]
fn read_fasta(
    path: &str,
    regions: Option<Vec<String>>,
    index: Option<String>,
    gzi: Option<String>,
    fields: Option<Vec<String>>,
) -> Vec<u8> {
    let compressed = path.ends_with(".gz");
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();
    let scanner = FastaScanner::new();

    let ipc = if let Some(regions) = regions {
        let index_path = index.unwrap_or(format!("{}.fai", path));
        let index =
            noodles::fasta::fai::read(index_path).expect("Could not read FASTA index file.");
        let regions: Vec<Region> = regions
            .into_iter()
            .map(|s| s.parse::<Region>().unwrap())
            .collect();
        if compressed {
            let gzi_path = gzi.unwrap_or(format!("{}.gzi", path));
            let gzindex =
                noodles::bgzf::gzi::read(gzi_path).expect("Could not read GZI index file.");
            let bgzf_reader = IndexedBgzfReader::new(reader, gzindex);
            let fmt_reader = noodles::fasta::io::Reader::new(bgzf_reader);
            let batches = scanner
                .scan_query(fmt_reader, regions, index, fields, None)
                .unwrap();
            batches_to_ipc(batches)
        } else {
            let fmt_reader = noodles::fasta::io::Reader::new(reader);
            let batches = scanner
                .scan_query(fmt_reader, regions, index, fields, None)
                .unwrap();
            batches_to_ipc(batches)
        }
    } else {
        let fmt_reader = noodles::fasta::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, fields, None, None).unwrap();
        batches_to_ipc(batches)
    };

    ipc.unwrap()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod oxbow;
    fn read_fasta;
    fn read_fastq;
}

use std::io::BufReader;

use extendr_api::prelude::*;

use flate2::bufread::MultiGzDecoder;
use noodles::bgzf::IndexedReader as IndexedBgzfReader;
use noodles::core::Region;

use oxbow::alignment::{BamScanner, SamScanner};
use oxbow::bbi::{BigBedScanner, BigWigScanner};
use oxbow::bed::BedScanner;
use oxbow::gxf::{GffScanner, GtfScanner};
use oxbow::sequence::{FastaScanner, FastqScanner};
use oxbow::util::batches_to_ipc;
use oxbow::variant::{BcfScanner, GenotypeBy, VcfScanner};

pub const BUFFER_SIZE_BYTES: usize = const { 1024 * 1024 };

/// Return Arrow IPC format from a FASTQ file.
#[extendr]
fn read_fastq_impl(path: &str, fields: Option<Vec<String>>) -> Vec<u8> {
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
#[extendr]
fn read_fasta_impl(
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

/// Return Arrow IPC format from a SAM file.
#[extendr]
pub fn read_sam_impl(
    path: &str,
    region: Option<String>,
    index: Option<String>,
    fields: Option<Vec<String>>,
    scan_rows: Option<usize>,
) -> Vec<u8> {
    let compressed = path.ends_with(".gz");
    let scan_rows = Some(scan_rows.unwrap_or(1024));
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();

    let ipc = if let Some(region) = region {
        let index_path = index.unwrap_or(format!("{}.tbi", path));
        let index = noodles::tabix::read(index_path).expect("Could not read TBI index file.");
        let region = region.parse::<Region>().unwrap();
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::sam::io::Reader::new(bgzf_reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = SamScanner::new(header);
        let tag_defs = scanner.tag_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan_query(
                fmt_reader,
                region,
                index,
                fields,
                Some(tag_defs),
                None,
                None,
            )
            .unwrap();
        batches_to_ipc(batches)
    } else if compressed {
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::sam::io::Reader::new(bgzf_reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = SamScanner::new(header);
        let tag_defs = scanner.tag_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan(fmt_reader, fields, Some(tag_defs), None, None)
            .unwrap();
        batches_to_ipc(batches)
    } else {
        let mut fmt_reader = noodles::sam::io::Reader::new(reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = SamScanner::new(header);
        let tag_defs = scanner.tag_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan(fmt_reader, fields, Some(tag_defs), None, None)
            .unwrap();
        batches_to_ipc(batches)
    };

    ipc.unwrap()
}

/// Return Arrow IPC format from a BAM file.
#[extendr]
pub fn read_bam_impl(
    path: &str,
    region: Option<String>,
    index: Option<String>,
    fields: Option<Vec<String>>,
    scan_rows: Option<usize>,
) -> Vec<u8> {
    let compressed = true;
    let scan_rows = Some(scan_rows.unwrap_or(1024));
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();

    let ipc = if let Some(region) = region {
        let index_path = index.unwrap_or(format!("{}.bai", path));
        let index = noodles::bam::bai::read(index_path).expect("Could not read BAI index file.");
        let region = region.parse::<Region>().unwrap();
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::bam::io::Reader::from(bgzf_reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = BamScanner::new(header);
        let tag_defs = scanner.tag_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan_query(
                fmt_reader,
                region,
                index,
                fields,
                Some(tag_defs),
                None,
                None,
            )
            .unwrap();
        batches_to_ipc(batches)
    } else if compressed {
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::bam::io::Reader::from(bgzf_reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = BamScanner::new(header);
        let tag_defs = scanner.tag_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan(fmt_reader, fields, Some(tag_defs), None, None)
            .unwrap();
        batches_to_ipc(batches)
    } else {
        let mut fmt_reader = noodles::bam::io::Reader::from(reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = BamScanner::new(header);
        let tag_defs = scanner.tag_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan(fmt_reader, fields, Some(tag_defs), None, None)
            .unwrap();
        batches_to_ipc(batches)
    };

    ipc.unwrap()
}

/// Return Arrow IPC format from a VCF file.
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn read_vcf_impl(
    path: &str,
    region: Option<String>,
    index: Option<String>,
    fields: Option<Vec<String>>,
    info_fields: Option<Vec<String>>,
    genotype_fields: Option<Vec<String>>,
    samples: Option<Vec<String>>,
    genotype_by: Option<String>,
) -> Vec<u8> {
    let compressed = path.ends_with(".gz");
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();
    let genotype_by = match genotype_by.unwrap().as_str() {
        "sample" => GenotypeBy::Sample,
        "field" => GenotypeBy::Field,
        _ => panic!("Invalid value for `genotype_by`. Must be 'sample' or 'field'."),
    };

    let ipc = if let Some(region) = region {
        let index_path = index.unwrap_or(format!("{}.tbi", path));
        let index = noodles::tabix::read(index_path).expect("Could not read TBI index file.");
        let region = region.parse::<Region>().unwrap();
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = VcfScanner::new(header);
        let batches = scanner
            .scan_query(
                fmt_reader,
                region,
                index,
                fields,
                info_fields,
                genotype_fields,
                samples,
                Some(genotype_by),
                None,
                None,
            )
            .unwrap();
        batches_to_ipc(batches)
    } else if compressed {
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = VcfScanner::new(header);
        let batches = scanner
            .scan(
                fmt_reader,
                fields,
                info_fields,
                genotype_fields,
                samples,
                Some(genotype_by),
                None,
                None,
            )
            .unwrap();
        batches_to_ipc(batches)
    } else {
        let mut fmt_reader = noodles::vcf::io::Reader::new(reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = VcfScanner::new(header);
        let batches = scanner
            .scan(
                fmt_reader,
                fields,
                info_fields,
                genotype_fields,
                samples,
                Some(genotype_by),
                None,
                None,
            )
            .unwrap();
        batches_to_ipc(batches)
    };

    ipc.unwrap()
}

/// Return Arrow IPC format from a BCF file.
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn read_bcf_impl(
    path: &str,
    region: Option<String>,
    index: Option<String>,
    fields: Option<Vec<String>>,
    info_fields: Option<Vec<String>>,
    genotype_fields: Option<Vec<String>>,
    samples: Option<Vec<String>>,
    genotype_by: Option<String>,
) -> Vec<u8> {
    let compressed = true;
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();
    let genotype_by = match genotype_by.unwrap().as_str() {
        "sample" => GenotypeBy::Sample,
        "field" => GenotypeBy::Field,
        _ => panic!("Invalid value for `genotype_by`. Must be 'sample' or 'field'."),
    };

    let ipc = if let Some(region) = region {
        let index_path = index.unwrap_or(format!("{}.csi", path));
        let index = noodles::csi::read(index_path).expect("Could not read CSI index file.");
        let region = region.parse::<Region>().unwrap();
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = BcfScanner::new(header);
        let batches = scanner
            .scan_query(
                fmt_reader,
                region,
                index,
                fields,
                info_fields,
                genotype_fields,
                samples,
                Some(genotype_by),
                None,
                None,
            )
            .unwrap();
        batches_to_ipc(batches)
    } else if compressed {
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = BcfScanner::new(header);
        let batches = scanner
            .scan(
                fmt_reader,
                fields,
                info_fields,
                genotype_fields,
                samples,
                Some(genotype_by),
                None,
                None,
            )
            .unwrap();
        batches_to_ipc(batches)
    } else {
        let mut fmt_reader = noodles::bcf::io::Reader::from(reader);
        let header = fmt_reader.read_header().unwrap();
        let scanner = BcfScanner::new(header);
        let batches = scanner
            .scan(
                fmt_reader,
                fields,
                info_fields,
                genotype_fields,
                samples,
                Some(genotype_by),
                None,
                None,
            )
            .unwrap();
        batches_to_ipc(batches)
    };

    ipc.unwrap()
}

/// Return Arrow IPC format from a GTF file.
#[extendr]
pub fn read_gtf_impl(
    path: &str,
    region: Option<String>,
    index: Option<String>,
    fields: Option<Vec<String>>,
    scan_rows: Option<usize>,
) -> Vec<u8> {
    let compressed = path.ends_with(".gz");
    let scan_rows = Some(scan_rows.unwrap_or(1024));
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();

    let ipc = if let Some(region) = region {
        let index_path = index.unwrap_or(format!("{}.tbi", path));
        let index = noodles::tabix::read(index_path).expect("Could not read TBI index file.");
        let region = region.parse::<Region>().unwrap();
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
        let scanner = GtfScanner::new(None);
        let attr_defs = scanner.attribute_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan_query(
                fmt_reader,
                region,
                index,
                fields,
                Some(attr_defs),
                None,
                None,
            )
            .unwrap();
        batches_to_ipc(batches)
    } else if compressed {
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
        let scanner = GtfScanner::new(None);
        let attr_defs = scanner.attribute_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan(fmt_reader, fields, Some(attr_defs), None, None)
            .unwrap();
        batches_to_ipc(batches)
    } else {
        let mut fmt_reader = noodles::gtf::io::Reader::new(reader);
        let scanner = GtfScanner::new(None);
        let attr_defs = scanner.attribute_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan(fmt_reader, fields, Some(attr_defs), None, None)
            .unwrap();
        batches_to_ipc(batches)
    };

    ipc.unwrap()
}

/// Return Arrow IPC format from a GFF file.
#[extendr]
pub fn read_gff_impl(
    path: &str,
    region: Option<String>,
    index: Option<String>,
    fields: Option<Vec<String>>,
    scan_rows: Option<usize>,
) -> Vec<u8> {
    let compressed = path.ends_with(".gz");
    let scan_rows = Some(scan_rows.unwrap_or(1024));
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();

    let ipc = if let Some(region) = region {
        let index_path = index.unwrap_or(format!("{}.tbi", path));
        let index = noodles::tabix::read(index_path).expect("Could not read TBI index file.");
        let region = region.parse::<Region>().unwrap();
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
        let scanner = GffScanner::new(None);
        let attr_defs = scanner.attribute_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan_query(
                fmt_reader,
                region,
                index,
                fields,
                Some(attr_defs),
                None,
                None,
            )
            .unwrap();
        batches_to_ipc(batches)
    } else if compressed {
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let mut fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
        let scanner = GffScanner::new(None);
        let attr_defs = scanner.attribute_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan(fmt_reader, fields, Some(attr_defs), None, None)
            .unwrap();
        batches_to_ipc(batches)
    } else {
        let mut fmt_reader = noodles::gff::io::Reader::new(reader);
        let scanner = GffScanner::new(None);
        let attr_defs = scanner.attribute_defs(&mut fmt_reader, scan_rows).unwrap();
        let batches = scanner
            .scan(fmt_reader, fields, Some(attr_defs), None, None)
            .unwrap();
        batches_to_ipc(batches)
    };

    ipc.unwrap()
}

/// Return Arrow IPC format from a BED file.
#[extendr]
pub fn read_bed_impl(
    path: &str,
    bed_schema: &str,
    region: Option<String>,
    index: Option<String>,
    fields: Option<Vec<String>>,
) -> Vec<u8> {
    let compressed = path.ends_with(".gz");
    let bed_schema = bed_schema.parse().unwrap();
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();

    let ipc = if let Some(region) = region {
        let index_path = index.unwrap_or(format!("{}.tbi", path));
        let index = noodles::tabix::read(index_path).expect("Could not read TBI index file.");
        let region = region.parse::<Region>().unwrap();
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
        let scanner = BedScanner::new(bed_schema);
        let batches = scanner
            .scan_query(fmt_reader, region, index, fields, None, None)
            .unwrap();
        batches_to_ipc(batches)
    } else if compressed {
        let bgzf_reader = noodles::bgzf::Reader::new(reader);
        let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
        let scanner = BedScanner::new(bed_schema);
        let batches = scanner.scan(fmt_reader, fields, None, None).unwrap();
        batches_to_ipc(batches)
    } else {
        let fmt_reader = noodles::bed::io::Reader::new(reader);
        let scanner = BedScanner::new(bed_schema);
        let batches = scanner.scan(fmt_reader, fields, None, None).unwrap();
        batches_to_ipc(batches)
    };

    ipc.unwrap()
}

/// Return Arrow IPC format from a BigWig file.
#[extendr]
pub fn read_bigwig_impl(
    path: &str,
    region: Option<String>,
    fields: Option<Vec<String>>,
) -> Vec<u8> {
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();

    let ipc = if let Some(region) = region {
        let region = region.parse::<Region>().unwrap();
        let fmt_reader = bigtools::BigWigRead::open(reader).unwrap();
        let info = fmt_reader.info().clone();
        let scanner = BigWigScanner::new(info);
        let batches = scanner
            .scan_query(fmt_reader, region, fields, None, None)
            .unwrap();
        batches_to_ipc(batches)
    } else {
        let fmt_reader = bigtools::BigWigRead::open(reader).unwrap();
        let info = fmt_reader.info().clone();
        let scanner = BigWigScanner::new(info);
        let batches = scanner.scan(fmt_reader, fields, None, None).unwrap();
        batches_to_ipc(batches)
    };

    ipc.unwrap()
}

/// Return Arrow IPC format from a BigBed file.
#[extendr]
pub fn read_bigbed_impl(
    path: &str,
    bed_schema: &str,
    region: Option<String>,
    fields: Option<Vec<String>>,
) -> Vec<u8> {
    let bed_schema = bed_schema.parse().unwrap();
    let reader = std::fs::File::open(path)
        .map(|f| BufReader::with_capacity(BUFFER_SIZE_BYTES, f))
        .unwrap();

    let ipc = if let Some(region) = region {
        let region = region.parse::<Region>().unwrap();
        let fmt_reader = bigtools::BigBedRead::open(reader).unwrap();
        let info = fmt_reader.info().clone();
        let scanner = BigBedScanner::new(bed_schema, info);
        let batches = scanner
            .scan_query(fmt_reader, region, fields, None, None)
            .unwrap();
        batches_to_ipc(batches)
    } else {
        let fmt_reader = bigtools::BigBedRead::open(reader).unwrap();
        let info = fmt_reader.info().clone();
        let scanner = BigBedScanner::new(bed_schema, info);
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
    fn read_fasta_impl;
    fn read_fastq_impl;
    fn read_sam_impl;
    fn read_bam_impl;
    fn read_vcf_impl;
    fn read_bcf_impl;
    fn read_gtf_impl;
    fn read_gff_impl;
    fn read_bed_impl;
    fn read_bigwig_impl;
    fn read_bigbed_impl;
}

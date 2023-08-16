use std::collections::HashSet;
use std::fs::File;
use std::io;

use noodles::bam::bai;
use noodles::csi::index::ReferenceSequence;
use noodles::{bgzf, csi, tabix};
// use noodles::cram::crai;

fn get_ref_last_position(rseq: &ReferenceSequence) -> bgzf::VirtualPosition {
    let rend = rseq
        .bins()
        .values()
        .map(|bin| bin.chunks().iter().map(|chunk| chunk.end()).max().unwrap())
        .max()
        .unwrap();
    rend
}

fn get_offsets_from_linear_index(rseq: &ReferenceSequence) -> Vec<(u64, u16)> {
    let mut offsets: Vec<(u64, u16)> = Vec::new();
    let rend = get_ref_last_position(&rseq);
    for offset in rseq.linear_index() {
        let a = offset.compressed();
        let b = offset.uncompressed();
        if a == 0 && b == 0 {
            continue;
        }
        offsets.push((a, b));
    }
    offsets.push((rend.compressed(), rend.uncompressed()));

    // Remove duplicates and sort by compressed, then uncompressed
    let offsets_uniq: HashSet<(u64, u16)> = offsets.drain(..).collect();
    let mut offsets_uniq: Vec<(u64, u16)> = offsets_uniq.into_iter().collect();
    offsets_uniq.sort_by_key(|&(a, b)| (a, b));
    offsets_uniq
}

fn get_offsets_from_binning_index(rseq: &ReferenceSequence) -> Vec<(u64, u16)> {
    let mut offsets: Vec<(u64, u16)> = Vec::new();
    let rend = get_ref_last_position(&rseq);
    for bin in rseq.bins().values() {
        for chunk in bin.chunks() {
            offsets.push((chunk.start().compressed(), chunk.start().uncompressed()));
            offsets.push((chunk.end().compressed(), chunk.end().uncompressed()));
        }
    }
    offsets.push((rend.compressed(), rend.uncompressed()));

    // Remove duplicates and sort by compressed, then uncompressed
    let offsets_uniq: HashSet<(u64, u16)> = offsets.drain(..).collect();
    let mut offsets_uniq: Vec<(u64, u16)> = offsets_uniq.into_iter().collect();
    offsets_uniq.sort_by_key(|&(a, b)| (a, b));
    offsets_uniq
}

fn consolidate_chunks(offsets: &Vec<(u64, u16)>, chunksize: u64) -> Vec<(u64, u16)> {
    let mut consolidated = Vec::new();
    let mut last_offset = 0;

    consolidated.push(offsets[0]);
    for &(offset, value) in &offsets[1..offsets.len() - 1] {
        if offset >= last_offset + chunksize {
            consolidated.push((offset, value));
            last_offset = offset;
        }
    }
    consolidated.push(offsets[offsets.len() - 1]);

    consolidated
}

fn partition_from_index(index: &csi::Index, chunksize: u64) -> Vec<(u64, u16)> {
    let mut partition: Vec<(u64, u16)> = Vec::new();
    for rseq in index.reference_sequences() {
        if rseq.bins().is_empty() {
            continue;
        }
        let offsets = if rseq.linear_index().is_empty() {
            get_offsets_from_binning_index(&rseq)
        } else {
            get_offsets_from_linear_index(&rseq)
        };
        let consolidated = consolidate_chunks(&offsets, chunksize);
        partition.extend(consolidated);
    }

    // Remove duplicates and sort by compressed, then uncompressed
    let partition: HashSet<(u64, u16)> = partition.drain(..).collect();
    let mut partition: Vec<(u64, u16)> = partition.into_iter().collect();
    partition.sort_by_key(|&(a, b)| (a, b));
    partition
}

fn parse_index_file(path: &str) -> io::Result<csi::Index> {
    if path.ends_with(".csi") {
        let mut reader = File::open(path).map(csi::Reader::new)?;
        reader.read_index()
    } else if path.ends_with(".tbi") {
        let mut reader = File::open(path).map(tabix::Reader::new)?;
        reader.read_index()
    } else if path.ends_with(".bai") {
        let mut reader = File::open(path).map(bai::Reader::new)?;
        reader.read_header()?;
        reader.read_index()
    // } else if path.ends_with(".crai") {
    //     let mut reader = File::open(path).map(crai::Reader::new)?;
    //     reader.read_header()?;
    //     return reader.read_index();
    } else {
        panic!("Index file must end with .csi, .tbi, .bai");
    }
}

pub fn partition_from_index_file(path: &str, chunksize: u64) -> Vec<(u64, u16)> {
    let index = parse_index_file(path).unwrap();
    let partition = partition_from_index(&index, chunksize);
    partition
}

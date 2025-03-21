use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use noodles::bam::bai;
use noodles::core::region::Interval;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::csi::binning_index::index::reference_sequence::index::{
    BinnedIndex as BinnedRefIndex, Index as RefIndex, LinearIndex as LinearRefIndex,
};
use noodles::csi::binning_index::index::reference_sequence::ReferenceSequence;
use noodles::csi::binning_index::index::Index;
use noodles::csi::BinningIndex;
use noodles::{bgzf, csi, tabix};

/// An enum representing the two types of binning indexes.
///
/// A binning index contains a sub-index for each reference sequence.
pub enum IndexType {
    Linear(Index<LinearRefIndex>),
    Binned(Index<BinnedRefIndex>),
}

impl IndexType {
    pub fn into_boxed(self) -> Box<dyn BinningIndex> {
        match self {
            IndexType::Linear(index) => Box::new(index),
            IndexType::Binned(index) => Box::new(index),
        }
    }

    pub fn query(
        &self,
        reference_sequence_id: usize,
        interval: Interval,
    ) -> io::Result<Vec<Chunk>> {
        match self {
            IndexType::Linear(index) => index.query(reference_sequence_id, interval),
            IndexType::Binned(index) => index.query(reference_sequence_id, interval),
        }
    }
}

/// Get the last virtual position of a reference sequence's chunks.
fn get_ref_last_position(rseq: &ReferenceSequence<impl RefIndex>) -> bgzf::VirtualPosition {
    let rend = rseq
        .bins()
        .values()
        .map(|bin| bin.chunks().iter().map(|chunk| chunk.end()).max().unwrap())
        .max()
        .unwrap();
    rend
}

/// Get all the virtual position offsets from a linear reference sequence index.
fn get_offsets_from_linear_index(rseq: &ReferenceSequence<LinearRefIndex>) -> Vec<(u64, u16)> {
    let mut offsets: Vec<(u64, u16)> = Vec::new();
    let rend = get_ref_last_position(rseq);
    for offset in rseq.index() {
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

/// Get the virtual positions of the chunk bounds of a binned reference sequence index.
fn get_offsets_from_binned_index(rseq: &ReferenceSequence<BinnedRefIndex>) -> Vec<(u64, u16)> {
    let mut offsets: Vec<(u64, u16)> = Vec::new();
    let rend = get_ref_last_position(rseq);
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

/// Merge overlapping or consecutive chunks until each chunk exceeds `chunksize`.
fn consolidate_chunks(offsets: &[(u64, u16)], chunksize: u64) -> Vec<(u64, u16)> {
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

/// Partition a BGZF file into roughly equal-sized chunks based on its index.
pub fn partition_from_index(index: &IndexType, chunksize: u64) -> Vec<(u64, u16)> {
    let mut partition: Vec<(u64, u16)> = Vec::new();

    match index {
        IndexType::Linear(ix) => {
            for rseq in ix.reference_sequences() {
                if rseq.bins().is_empty() {
                    continue;
                }
                let offsets = get_offsets_from_linear_index(rseq);
                let consolidated = consolidate_chunks(&offsets, chunksize);
                partition.extend(consolidated);
            }
        }
        IndexType::Binned(ix) => {
            for rseq in ix.reference_sequences() {
                if rseq.bins().is_empty() {
                    continue;
                }
                let offsets = get_offsets_from_binned_index(rseq);
                let consolidated = consolidate_chunks(&offsets, chunksize);
                partition.extend(consolidated);
            }
        }
    }

    // Remove duplicates and sort by compressed, then uncompressed
    let partition: HashSet<(u64, u16)> = partition.drain(..).collect();
    let mut partition: Vec<(u64, u16)> = partition.into_iter().collect();
    partition.sort_by_key(|&(a, b)| (a, b));
    partition
}

fn parse_index_file(path: &str) -> io::Result<IndexType> {
    if path.ends_with(".csi") {
        let mut reader = File::open(path).map(csi::Reader::new)?;
        Ok(IndexType::Binned(reader.read_index()?))
    } else if path.ends_with(".tbi") {
        let mut reader = File::open(path).map(tabix::Reader::new)?;
        Ok(IndexType::Linear(reader.read_index()?))
    } else if path.ends_with(".bai") {
        let mut reader = File::open(path).map(bai::Reader::new)?;
        Ok(IndexType::Linear(reader.read_index()?))
    // } else if path.ends_with(".crai") {
    //     let mut reader = File::open(path).map(crai::Reader::new)?;
    //     reader.read_header()?;
    //     return reader.read_index();
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Index file must end with .csi, .tbi, or .bai",
        ))
    }
}

pub fn partition_from_index_file(path: &str, chunksize: u64) -> Vec<(u64, u16)> {
    let index = parse_index_file(path).unwrap();
    partition_from_index(&index, chunksize)
}

pub fn index_from_reader<R>(mut read: R) -> io::Result<IndexType>
where
    R: Read + Seek,
{
    // Unlike .tbi and .csi, .bai is not bgzf-compressed
    // so we read off the magic directly.
    let mut magic = [0; 4];
    read.read_exact(&mut magic)?;
    read.seek(SeekFrom::Start(0))?;
    if magic == b"BAI\x01" as &[u8] {
        let mut bai_reader = noodles::bam::bai::Reader::new(read);
        let index = bai_reader.read_index()?;
        Ok(IndexType::Linear(index))
    } else {
        let mut csi_reader = noodles::csi::Reader::new(read);
        match csi_reader.read_index() {
            Ok(index) => Ok(IndexType::Binned(index)),
            Err(_) => {
                let mut read = csi_reader.into_inner().into_inner();
                read.seek(SeekFrom::Start(0))?;
                let mut tabix_reader = noodles::tabix::Reader::new(read);
                let index = tabix_reader.read_index()?;
                Ok(IndexType::Linear(index))
            }
        }
    }
}

pub fn index_from_path(path: &str) -> io::Result<IndexType> {
    let bai_path = format!("{}.bai", path);
    let csi_path = format!("{}.csi", path);
    let tbi_path = format!("{}.tbi", path);
    let index = if Path::new(&bai_path).exists() {
        IndexType::Linear(noodles::bam::bai::read(bai_path)?)
    } else if Path::new(&csi_path).exists() {
        IndexType::Binned(noodles::csi::read(csi_path)?)
    } else if Path::new(&tbi_path).exists() {
        IndexType::Linear(noodles::tabix::read(tbi_path)?)
    } else {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "Could not find a .bai or .csi index file for the given file.",
        ));
    };
    Ok(index)
}

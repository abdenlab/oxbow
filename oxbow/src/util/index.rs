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
    pub fn from_path(path: &str) -> io::Result<Self> {
        if path.ends_with(".csi") {
            let mut reader = File::open(path).map(csi::io::Reader::new)?;
            Ok(IndexType::Binned(reader.read_index()?))
        } else if path.ends_with(".tbi") {
            let mut reader = File::open(path).map(tabix::io::Reader::new)?;
            Ok(IndexType::Linear(reader.read_index()?))
        } else if path.ends_with(".bai") {
            let mut reader = File::open(path).map(bai::io::Reader::new)?;
            Ok(IndexType::Linear(reader.read_index()?))
        } else {
            Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Index file must end with .csi, .tbi, or .bai",
            ))
        }
    }

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

/// Get the sorted decoded virtual positions from a linear reference sequence index.
fn get_offsets_from_linear_index(rseq: &ReferenceSequence<LinearRefIndex>) -> Vec<(u64, u16)> {
    let mut offsets: Vec<(u64, u16)> = Vec::new();
    let rend = get_ref_last_position(rseq);
    for offset in rseq.index() {
        let u = offset.compressed();
        let c = offset.uncompressed();
        if u == 0 && c == 0 {
            continue;
        }
        offsets.push((u, c));
    }
    offsets.push((rend.compressed(), rend.uncompressed()));

    // Remove duplicates and sort by compressed, then uncompressed
    let offsets_uniq: HashSet<(u64, u16)> = offsets.drain(..).collect();
    let mut offsets_uniq: Vec<(u64, u16)> = offsets_uniq.into_iter().collect();
    offsets_uniq.sort_by_key(|&(u, c)| (u, c));
    offsets_uniq
}

/// Get the sorted decoded virtual positions from the chunk bounds of a binned reference sequence index.
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
    offsets_uniq.sort_by_key(|&(u, c)| (u, c));
    offsets_uniq
}

/// Prune out boundaries from a virtual position partition until pieces meet or exceed `chunksize`
/// compressed bytes.
fn prune_partition(vpositions: &[(u64, u16)], chunksize: u64) -> Vec<(u64, u16)> {
    let mut pruned = Vec::new();
    let mut last_compressed = 0;

    pruned.push(vpositions[0]);
    for &(compressed, uncompressed) in &vpositions[1..vpositions.len() - 1] {
        if compressed >= last_compressed + chunksize {
            pruned.push((compressed, uncompressed));
            last_compressed = compressed;
        }
    }
    pruned.push(vpositions[vpositions.len() - 1]);

    pruned
}

/// Partition a BGZF file into roughly equal chunks based on the known record bounds in its index.
///
/// If `chunksize` is 0, all index boundary positions are returned (sorted and deduplicated).
pub fn partition_from_index(index: &IndexType, chunksize: u64) -> Vec<(u64, u16)> {
    let mut partition: Vec<(u64, u16)> = Vec::new();

    match index {
        IndexType::Linear(ix) => {
            for rseq in ix.reference_sequences() {
                if rseq.bins().is_empty() {
                    continue;
                }
                let offsets = get_offsets_from_linear_index(rseq);
                let pruned_offsets = if chunksize > 0 {
                    prune_partition(&offsets, chunksize)
                } else {
                    offsets
                };
                partition.extend(pruned_offsets);
            }
        }
        IndexType::Binned(ix) => {
            for rseq in ix.reference_sequences() {
                if rseq.bins().is_empty() {
                    continue;
                }
                let offsets = get_offsets_from_binned_index(rseq);
                let pruned_offsets = if chunksize > 0 {
                    prune_partition(&offsets, chunksize)
                } else {
                    offsets
                };
                partition.extend(pruned_offsets);
            }
        }
    }

    // Remove duplicates and sort by compressed, then uncompressed across all reference sequences
    let partition: HashSet<(u64, u16)> = partition.drain(..).collect();
    let mut partition: Vec<(u64, u16)> = partition.into_iter().collect();
    partition.sort_by_key(|&(u, c)| (u, c));

    partition
}

/// Load an index from a file reader, returning the appropriate IndexType.
pub fn index_from_reader<R>(mut read: R) -> io::Result<IndexType>
where
    R: Read + Seek,
{
    // BAI is not BGZF-compressed so we read off its candidate magic directly.
    let mut magic = [0; 4];
    read.read_exact(&mut magic)?;
    read.seek(SeekFrom::Start(0))?;

    if magic == b"BAI\x01" as &[u8] {
        let mut bai_reader = noodles::bam::bai::io::Reader::new(read);
        return Ok(IndexType::Linear(bai_reader.read_index()?));
    }

    // Try CSI, then TBI if CSI fails.
    let mut csi_reader = noodles::csi::io::Reader::new(read);
    match csi_reader.read_index() {
        Ok(index) => Ok(IndexType::Binned(index)),
        Err(_) => {
            let mut read = csi_reader.into_inner().into_inner();
            read.seek(SeekFrom::Start(0))?;
            let mut tabix_reader = noodles::tabix::io::Reader::new(read);
            match tabix_reader.read_index() {
                Ok(index) => Ok(IndexType::Linear(index)),
                Err(_) => Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Failed to read index from reader.",
                )),
            }
        }
    }
}

/// Find an index file for a given source path and load it.
pub fn index_from_source_path(path: &str) -> io::Result<IndexType> {
    let bai_path = format!("{}.bai", path);
    let csi_path = format!("{}.csi", path);
    let tbi_path = format!("{}.tbi", path);
    let index = if Path::new(&bai_path).exists() {
        IndexType::from_path(&bai_path)?
    } else if Path::new(&csi_path).exists() {
        IndexType::from_path(&csi_path)?
    } else if Path::new(&tbi_path).exists() {
        IndexType::from_path(&tbi_path)?
    } else {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "Could not find a .bai, .csi, or .tbi index file for the given file.",
        ));
    };
    Ok(index)
}

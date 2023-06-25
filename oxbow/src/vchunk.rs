#![allow(dead_code)]
use std::ops::RangeFull;

use noodles::{bgzf, core, csi::{self, index::reference_sequence::bin::Chunk}};
use noodles::core::Position;

// http://www.htslib.org/doc/bgzip.html
const ESTIMATED_BLOCK_SIZE: usize = 64 * 1024; // 64 KB

pub struct IndexOptimizer {
    index: csi::Index,
    max_block_size: usize,
}

impl IndexOptimizer {

    /// Convert a query interval to a sequence of index chunks.
    ///
    /// The query interval is converted to index bins, which are then flattened into a
    /// sequence of chunks. The chunks are then combined (optimized) to be larger within
    /// a memory limit.
    ///
    /// Adapted from https://github.com/zaeleus/noodles/blob/ec91bf5a5adff9bf4f2dc02c08b815d02fddad18/noodles-csi/src/index.rs#L114-L140
    pub fn query<I>(
        &self,
        reference_sequence_id: usize,
        interval: I,
    ) -> std::io::Result<Vec<Chunk>>
    where
        I: Into<core::region::Interval>,
    {
        let interval = interval.into();
        let reference_sequence = match self.index.reference_sequences().get(reference_sequence_id) {
            Some(reference_sequence) => reference_sequence,
            None => return Ok(Vec::new()),
        };

        let query_bins = reference_sequence
            .query(self.index.min_shift(), self.index.depth(), interval)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;

        let chunks: Vec<_> = query_bins
            .iter()
            .flat_map(|bin| bin.chunks())
            .copied()
            .collect();

        // https://github.com/zaeleus/noodles/blob/ec91bf5a5adff9bf4f2dc02c08b815d02fddad18/noodles-csi/src/index.rs#L161
        let start = interval.start().unwrap_or(Position::MIN);
        let min_offset =
            reference_sequence.min_offset(self.index.min_shift(), self.index.depth(), start);

        Ok(optimize_chunks(&chunks, min_offset, self.max_block_size))
    }

    pub fn query_all(&self) -> std::io::Result<Vec<Chunk>> {
        Ok(self
            .index
            .reference_sequences()
            .iter()
            .enumerate()
            .map(|(i, _)| self.query(i, RangeFull))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
            .collect())
    }
}

impl From<csi::Index> for IndexOptimizer {
    fn from(index: csi::Index) -> Self {
        Self {
            index,
            max_block_size: 1024 * 1024 * 100, // 100 MB
        }
    }
}

// Adapted from https://github.com/zaeleus/noodles/blob/ec91bf5a5adff9bf4f2dc02c08b815d02fddad18/noodles-csi/src/binning_index.rs#L72
pub fn optimize_chunks(
    chunks: &[Chunk],
    min_offset: bgzf::VirtualPosition,
    max_block_size: usize,
) -> Vec<Chunk> {
    let max_bgzf_blocks_per_chunk = max_block_size / ESTIMATED_BLOCK_SIZE;
    let mut chunks: Vec<_> = chunks
        .iter()
        .filter(|c| c.end() > min_offset)
        .copied()
        .collect();

    if chunks.is_empty() {
        return chunks;
    }

    chunks.sort_unstable_by_key(|c| c.start());

    // At worst, no chunks are merged, and the resulting list will be the same size as the input.
    let mut merged_chunks = Vec::with_capacity(chunks.len());

    // `chunks` is guaranteed to be non-empty.
    let mut current_chunk = chunks[0];

    let mut i = 0;
    for next_chunk in chunks.iter().skip(1) {
        if next_chunk.start() > current_chunk.end() || i >= max_bgzf_blocks_per_chunk {
            merged_chunks.push(current_chunk);
            current_chunk = *next_chunk;
            i = 0;
        } else if current_chunk.end() < next_chunk.end() {
            current_chunk = Chunk::new(current_chunk.start(), next_chunk.end());
            i += 1;
        }
    }

    merged_chunks.push(current_chunk);
    merged_chunks
}


#[cfg(test)]
mod tests {
    use super::*;
    use noodles::{bam,sam};
    // use noodles::core::Region;
    use std::io::{Read, Seek};

    type BamReader = bam::Reader<bgzf::Reader<std::io::Cursor<Vec<u8>>>>;

    fn get_reader() -> std::io::Result<(sam::Header, BamReader)> {
        let mut path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("fixtures/sample.bam");

        let header = std::fs::File::open(path.clone())
            .map(bam::Reader::new)?
            .read_header()?;

        let index = bam::bai::read(path.with_extension("bam.bai"))?;
        let chunks = IndexOptimizer::from(index).query_all()?;
        let chunk = chunks.first().unwrap();

        let mut buf: Vec<u8> = Vec::new();
        let mut file = std::fs::File::open(path)?;
        let (cstart, _upos) = chunk.start().into();
        let (cend, _uend) = chunk.end().into();
        file.seek(std::io::SeekFrom::Start(cstart))?;
        buf.resize((cend - cstart) as usize, 0);
        file.read_exact(&mut buf)?;

        Ok((header, bam::Reader::new(std::io::Cursor::new(buf))))
    }

    #[test]
    fn it_works() {
        let (header, mut reader) = get_reader().unwrap();
        let mut record = noodles::sam::alignment::Record::default();
        let size = reader.read_record(&header, &mut record).unwrap();
        println!("{:?}, size: {}", record, size);
        let size = reader.read_record(&header, &mut record).unwrap();
        println!("{:?}, size: {}", record, size);
        let size = reader.read_record(&header, &mut record).unwrap();
        println!("{:?}, size: {}", record, size);

        assert!(noodles::sam::alignment::Record::default() != record);
    }

}

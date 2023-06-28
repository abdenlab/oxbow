#![allow(dead_code)]
use std::ops::RangeFull;

use noodles::{bgzf, core, csi::{self, index::reference_sequence::bin::Chunk}};
use noodles::core::Position;

// http://www.htslib.org/doc/bgzip.html
const ESTIMATED_BLOCK_SIZE: usize = 64 * 1024; // 64 KB

pub struct IndexOptimizer<'a> {
    index: &'a csi::Index,
    max_block_size: usize,
}

impl IndexOptimizer<'_> {

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

        Ok(chunks)

        // Ok(optimize_chunks(&chunks, min_offset, self.max_block_size))
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

impl <'a> IndexOptimizer<'a> {
    fn new(index: &'a csi::Index) -> Self {
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
    use std::io::{self, Read, Seek};
    use std::fs;
    use noodles::{bam, sam};

    const REGION: &str = "chr1:100000-8000000";

    fn indexed_reader() -> io::Result<bam::IndexedReader<bgzf::Reader<fs::File>>> {
        let path: std::path::PathBuf = [env!("CARGO_MANIFEST_DIR"), "fixtures/sample.bam"].iter().collect();
        noodles::bam::indexed_reader::Builder::default().build_from_path(path)
    }

    fn read_range<R: Read + Seek>(target: &mut R, chunk: &Chunk) -> io::Result<Vec<u8>> {
        let mut buf: Vec<u8> = Vec::new();
        println!("reading range {:#?}", chunk);
        target.seek(std::io::SeekFrom::Start(chunk.start().compressed()))?;
        buf.resize((chunk.end().compressed() - chunk.start().compressed()) as usize, 0);
        target.read_exact(&mut buf)?;
        Ok(buf)
    }

    #[test]
    fn it_works() {
        let mut indexed_reader = indexed_reader().unwrap();
        let header = indexed_reader.read_header().unwrap();

        let region: noodles::core::Region = REGION.parse().unwrap();

        let chunks = IndexOptimizer::new(indexed_reader.index())
            .query(
                header.reference_sequences().get_index_of(region.name()).unwrap(),
                region.interval(),
            )
            .unwrap();

        chunks
            .iter()
            .for_each(|c| {
                println!(
                    "start: {:?}, end: {:?}",
                    (c.start().compressed(), c.start().uncompressed()),
                    (c.end().compressed(), c.end().uncompressed()),
                )
            });

        let chunk = chunks.first().unwrap();
        let buf = read_range(indexed_reader.get_mut().get_mut(), chunk).unwrap();

        let (_cpos, upos) = chunk.start().into();

        {
            let mut reader = bam::Reader::new(io::Cursor::new(buf));
            reader.seek(bgzf::VirtualPosition::try_from((0, upos)).unwrap()).unwrap();

            let mut record = sam::alignment::Record::default();
            let size = reader.read_record(&header, &mut record).unwrap();
            println!("{:?}, size: {}", record, size);
            let size = reader.read_record(&header, &mut record).unwrap();
            println!("{:?}, size: {}", record, size);
            let size = reader.read_record(&header, &mut record).unwrap();
            println!("{:?}, size: {}", record, size);
        }

        assert!(false);
    }

    #[test]
    fn test2() {
        let mut indexed_reader = indexed_reader().unwrap();
        let header = indexed_reader.read_header().unwrap();
        let query = indexed_reader.query(&header, &REGION.parse().unwrap()).unwrap();
        assert_eq!(query.count(), 2);
    }

}

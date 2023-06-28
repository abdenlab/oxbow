#![allow(dead_code)]
use noodles::{bgzf, core, csi::{self, index::reference_sequence::bin::Chunk}};
use noodles::core::Position;

// http://www.htslib.org/doc/bgzip.html
const ESTIMATED_BLOCK_SIZE: usize = 64 * 1024; // 64 KB

pub struct IndexOptimizer<'a> {
    index: &'a csi::Index,
    max_block_size: usize,
}

impl <'a> IndexOptimizer<'a> {
    fn new(index: &'a csi::Index) -> Self {
        Self {
            index,
            max_block_size: 1024 * 1024 * 100, // 100 MB
        }
    }
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

        let reference_sequence = self
            .index
            .reference_sequences()
            .get(reference_sequence_id)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!("invalid reference sequence ID: {reference_sequence_id}"),
                )
            })?;

        let query_bins = reference_sequence
            .query(self.index.min_shift(), self.index.depth(), interval)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;

        let chunks: Vec<_> = query_bins
            .iter()
            .flat_map(|bin| bin.chunks())
            .copied()
            .collect();

        println!("Chunks:");
        chunks.iter().for_each(print_chunk);

        let (start, _) = resolve_interval(self.index.min_shift(), self.index.depth(), interval)?;
        let min_offset = reference_sequence.min_offset(self.index.min_shift(), self.index.depth(), start);
        let merged_chunks = optimize_chunks(&chunks, min_offset);

        println!("Merged chunks:");
        merged_chunks.iter().for_each(print_chunk);

        println!("Finished");

        Ok(merged_chunks)
    }

}

pub fn optimize_chunks(chunks: &[Chunk], min_offset: bgzf::VirtualPosition) -> Vec<Chunk> {
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

    for next_chunk in chunks.iter().skip(1) {
        if next_chunk.start() > current_chunk.end() {
            merged_chunks.push(current_chunk);
            current_chunk = *next_chunk;
        } else if current_chunk.end() < next_chunk.end() {
            current_chunk = Chunk::new(current_chunk.start(), next_chunk.end());
        }
    }

    merged_chunks.push(current_chunk);

    merged_chunks
}

fn resolve_interval<I>(min_shift: u8, depth: u8, interval: I) -> std::io::Result<(Position, Position)>
where
    I: Into<core::region::Interval>,
{
    let interval = interval.into();

    let start = interval.start().unwrap_or(Position::MIN);

    let max_position = {
        let n = (1 << (usize::from(min_shift) + 3 * usize::from(depth))) - 1;
        Position::try_from(n).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))
    }?;

    if start > max_position {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "invalid start bound",
        ));
    }

    let end = interval.end().unwrap_or(max_position);

    if end > max_position {
        Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "invalid end bound",
        ))
    } else {
        Ok((start, end))
    }
}

fn print_chunk(chunk: &Chunk) {
    println!(
        "Chunk(start={:?}, end={:?})",
        (chunk.start().compressed(), chunk.start().uncompressed()),
        (chunk.end().compressed(), chunk.end().uncompressed()),
    );
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{self, Read, Seek};
    use std::{fs, path};
    use noodles::{bam, sam};

    const REGION: &str = "chr1:100000-1000000";

    fn indexed_reader() -> io::Result<bam::IndexedReader<bgzf::Reader<fs::File>>> {
        let mut path = path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        // Example data from ...
        // wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam -O fixtures/encode.bam
        // wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam.bai -O fixtures/encode.bam.bai
        path.push("fixtures/encode.bam");
        noodles::bam::indexed_reader::Builder::default().build_from_path(path)
    }

    fn read_range<R: Read + Seek>(target: &mut R, chunk: &Chunk) -> io::Result<Vec<u8>> {
        let mut buf: Vec<u8> = Vec::new();
        target.seek(std::io::SeekFrom::Start(chunk.start().compressed()))?;
        buf.resize((chunk.end().compressed() - chunk.start().compressed()) as usize, 0);
        target.read_exact(&mut buf)?;
        Ok(buf)
    }

    #[test]
    fn test_optimizer() {
        let mut indexed_reader = indexed_reader().unwrap();
        let header = indexed_reader.read_header().unwrap();
        let region: noodles::core::Region = REGION.parse().unwrap();
        let chunks = IndexOptimizer::new(indexed_reader.index())
            .query(
                header.reference_sequences().get_index_of(region.name()).unwrap(),
                region.interval(),
            )
            .unwrap();

        let total: usize = chunks.iter()
            .map(|chunk| {
                let bytes_reader = read_range(indexed_reader.get_mut().get_mut(), chunk)
                    .map(io::Cursor::new)
                    .unwrap();
                let mut reader = bam::Reader::new(bytes_reader);
                let pos = bgzf::VirtualPosition::try_from((0, chunk.start().uncompressed())).unwrap();
                reader.seek(pos).unwrap();
                let total = reader.records(&header).count();
                total
            })
            .sum();

        assert_eq!(total, 452);
    }

    #[test]
    fn test_regular() {
        let mut indexed_reader = indexed_reader().unwrap();
        let header = indexed_reader.read_header().unwrap();
        let query = indexed_reader.query(&header, &REGION.parse().unwrap()).unwrap();
        assert_eq!(query.count(), 452);
    }

}

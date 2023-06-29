#![allow(dead_code)]
use noodles::{bgzf, sam::{self, alignment::Record}, core::{self, region::Interval}, csi::{self, index::reference_sequence::bin::Chunk}};
use noodles::core::Position;

// http://www.htslib.org/doc/bgzip.html
const ESTIMATED_BLOCK_SIZE: usize = 64 * 1024; // 64 KB

pub struct IndexOptimizer<'a> {
    resolver: &'a dyn Resolver,
    index: &'a csi::Index,
    max_block_size: usize,
}

impl <'a> IndexOptimizer<'a> {
    fn new(resolver: &'a dyn Resolver, index: &'a csi::Index) -> Self {
        Self {
            resolver,
            index,
            max_block_size: 1024 * 1024 * 100, // 100 MB
        }
    }
}

pub trait Resolver {
    fn resolve(&self, region: &core::region::Region) -> usize;
}

impl Resolver for sam::Header {

    fn resolve(&self, region: &core::region::Region) -> usize {
        self.reference_sequences().get_index_of(region.name()).unwrap()
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
    pub fn query(&self, region: &core::region::Region) -> std::io::Result<Vec<Chunk>> {
        let reference_sequence_id = self.resolver.resolve(region);
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

        let interval = region.interval();
        let query_bins = reference_sequence
            .query(self.index.min_shift(), self.index.depth(), interval)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidInput, e))?;

        let chunks: Vec<_> = query_bins
            .iter()
            .flat_map(|bin| bin.chunks())
            .copied()
            .collect();

        let (start, _) = resolve_interval(self.index.min_shift(), self.index.depth(), interval)?;
        let min_offset = reference_sequence.min_offset(self.index.min_shift(), self.index.depth(), start);
        let merged_chunks = optimize_chunks(&chunks, min_offset);

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

fn intersects(
    record: &Record,
    reference_sequence_id: usize,
    region_interval: Interval
) -> bool {
    match (
        record.reference_sequence_id(),
        record.alignment_start(),
        record.alignment_end(),
    ) {
        (Some(id), Some(start), Some(end)) => {
            let alignment_interval = (start..=end).into();
            id == reference_sequence_id && region_interval.intersects(alignment_interval)
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{self, Read, Seek};
    use std::ops::Range;
    use std::{fs, path};
    use noodles::bam;

    const REGION: &str = "chr2:100000-2000000";

    fn indexed_reader() -> io::Result<bam::IndexedReader<bgzf::Reader<fs::File>>> {
        let mut path = path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        // Example data from ...
        // wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam -O fixtures/encode.bam
        // wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam.bai -O fixtures/encode.bam.bai
        path.push("fixtures/encode.bam");
        noodles::bam::indexed_reader::Builder::default().build_from_path(path)
    }

    fn read_range<R: Read + Seek>(target: &mut R, range: Range<u64>) -> io::Result<Vec<u8>> {
        let mut buf: Vec<u8> = Vec::new();
        target.seek(std::io::SeekFrom::Start(range.start))?;
        buf.resize((range.end - range.start) as usize, 0);
        target.read_exact(&mut buf)?;
        Ok(buf)
    }

    #[test]
    fn it_works() {
        let mut indexed_reader = indexed_reader().unwrap();
        let header = indexed_reader.read_header().unwrap();
        let region = REGION.parse().unwrap();
        let chunks = IndexOptimizer::new(&header, indexed_reader.index())
            .query(&region)
            .unwrap();

        let block_starts: std::collections::HashSet<_> = indexed_reader
            .index()
            .reference_sequences()
            .iter()
            .flat_map(|r| r.bins().values().flat_map(|b| b.chunks()))
            .map(|c| c.start().compressed())
            .collect();

        let mut block_starts = Vec::from_iter(block_starts);
        block_starts.sort_unstable();

        let total: usize = chunks
            .iter()
            .map(|chunk| {
                let values = &block_starts;
                let start = values.binary_search(&chunk.start().compressed()).unwrap_or_else(|x| x);
                let end = values.binary_search(&chunk.end().compressed()).unwrap_or_else(|x| x);
                (chunk, values[start]..(values[end + 1] + 1))
            })
            .map(|(chunk, range)| {
                let bytes_reader = read_range(indexed_reader.get_mut().get_mut(), range)
                    .map(io::Cursor::new)
                    .unwrap();
                let mut reader = bam::Reader::new(bytes_reader);
                let pos = bgzf::VirtualPosition::try_from((0, chunk.start().uncompressed())).unwrap();
                reader.seek(pos).unwrap();
                // let record = sam::alignment::Record::default();
                // match self.reader.read_record(self.header, &mut self.record) {
                //     Ok(0) => None,
                //     Ok(_) => Some(Ok(self.record.clone())),
                //     Err(e) => Some(Err(e)),
                // }
                let records: Vec<_> = reader.records(&header)
                    .filter(|r| match r {
                        // idk why but ignoring UnexpectedEof is necessary
                        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => false,
                        Err(e) => panic!("{:?}", e),
                        Ok(r) => intersects(r, header.resolve(&region), region.interval())
                    })
                    .collect();
                records.len()
            })
            .sum();

        assert_eq!(
            total,
            indexed_reader.query(&header, &region).unwrap().count(),
            "total number of records in chunks should be equal to total number of records in query"
        );
    }

}

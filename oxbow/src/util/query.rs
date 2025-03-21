use std::io::{self, BufRead, Read};
use std::vec;

use noodles::bgzf::VirtualPosition;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;

/// A reader that decodes a sequence of virtual position ranges ("chunks") in a BGZF-encoded file.
///
/// The chunks normally come from querying an index file with a genomic range.
///
/// # Examples
///
/// ```no_run
/// use oxbow::util::query::BgzfChunkReader;
/// use std::fs::File;
/// use std::io::{self, BufReader};
/// use noodles::bgzf::VirtualPosition;
/// use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
///
/// let inner = File::open("sample.bam").unwrap();
/// let bgzf_reader = noodles::bgzf::Reader::new(BufReader::new(inner));
/// let chunks = vec![
///     Chunk::new(VirtualPosition::from(0), VirtualPosition::from(6)),
///     Chunk::new(VirtualPosition::from(12), VirtualPosition::from(18)),
/// ];
/// let mut chunk_reader = BgzfChunkReader::new(bgzf_reader, chunks);
/// ```
///
/// ```no_run
/// use oxbow::util::query::BgzfChunkReader;
/// use std::fs::File;
/// use std::io::{self, BufReader};
/// use noodles::core::region::Interval;
/// use noodles::csi::BinningIndex;
///
/// let inner = File::open("sample.bam").unwrap();
/// let index = noodles::bam::bai::read("sample.bam.bai").unwrap();
/// let chrom_id = 4;
/// let interval = "101-200".parse::<Interval>().unwrap();
/// let bgzf_reader = noodles::bgzf::Reader::new(BufReader::new(inner));
/// let chunks = index.query(chrom_id, interval).unwrap();
/// let mut chunk_reader = BgzfChunkReader::new(bgzf_reader, chunks);
/// ```
pub struct BgzfChunkReader<R> {
    reader: R,
    chunks: vec::IntoIter<Chunk>,
    state: State,
}

enum State {
    Seek,
    Read(VirtualPosition),
    Done,
}

impl<R> BgzfChunkReader<R>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    pub fn new(reader: R, chunks: Vec<Chunk>) -> Self {
        Self {
            reader,
            chunks: chunks.into_iter(),
            state: State::Seek,
        }
    }
}

impl<R> Read for BgzfChunkReader<R>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }
}

impl<R> BufRead for BgzfChunkReader<R>
where
    R: noodles::bgzf::io::BufRead + noodles::bgzf::io::Seek,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        loop {
            match self.state {
                State::Seek => {
                    self.state = match self.chunks.next() {
                        Some(chunk) => {
                            self.reader.seek_to_virtual_position(chunk.start())?;
                            State::Read(chunk.end())
                        }
                        None => State::Done,
                    }
                }
                State::Read(chunk_end) => {
                    if self.reader.virtual_position() < chunk_end {
                        return self.reader.fill_buf();
                    } else {
                        self.state = State::Seek;
                    }
                }
                State::Done => return Ok(&[]),
            }
        }
    }

    fn consume(&mut self, amt: usize) {
        self.reader.consume(amt);
    }
}

use std::io::{self, BufRead, Read, Seek};
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

/// A reader that reads a sequence of byte ranges from a seekable file.
///
/// The ranges are specified as 2-tuples of (start, end) byte positions.
/// This is the non-BGZF equivalent of [`BgzfChunkReader`], designed for
/// reading discontinuous regions from uncompressed files.
///
/// # Examples
///
/// ```no_run
/// use oxbow::util::query::ByteRangeReader;
/// use std::fs::File;
/// use std::io::BufReader;
///
/// let inner = File::open("sample.txt").unwrap();
/// let buffered_reader = BufReader::new(inner);
/// let ranges = vec![
///     (0, 100),      // Read bytes 0-100
///     (200, 300),    // Read bytes 200-300
/// ];
/// let mut range_reader = ByteRangeReader::new(buffered_reader, ranges);
/// ```
pub struct ByteRangeReader<R> {
    reader: R,
    ranges: vec::IntoIter<(u64, u64)>,
    state: ByteRangeState,
    current_pos: u64,
}

enum ByteRangeState {
    Seek,
    Read(u64), // end position
    Done,
}

impl<R> ByteRangeReader<R>
where
    R: BufRead + Seek,
{
    pub fn new(reader: R, ranges: Vec<(u64, u64)>) -> Self {
        Self {
            reader,
            ranges: ranges.into_iter(),
            state: ByteRangeState::Seek,
            current_pos: 0,
        }
    }
}

impl<R> Read for ByteRangeReader<R>
where
    R: BufRead + Seek,
{
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let mut src = self.fill_buf()?;
        let amt = src.read(buf)?;
        self.consume(amt);
        Ok(amt)
    }
}

impl<R> BufRead for ByteRangeReader<R>
where
    R: BufRead + Seek,
{
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        loop {
            match self.state {
                ByteRangeState::Seek => match self.ranges.next() {
                    Some((start, end)) => {
                        self.reader.seek(io::SeekFrom::Start(start))?;
                        self.current_pos = start;
                        self.state = ByteRangeState::Read(end);
                    }
                    None => {
                        self.state = ByteRangeState::Done;
                        return Ok(&[]);
                    }
                },
                ByteRangeState::Read(range_end) => {
                    let remaining = range_end.saturating_sub(self.current_pos) as usize;
                    if remaining == 0 {
                        self.state = ByteRangeState::Seek;
                        continue;
                    }

                    let buf = self.reader.fill_buf()?;
                    let available = buf.len().min(remaining);
                    return Ok(&buf[..available]);
                }
                ByteRangeState::Done => return Ok(&[]),
            }
        }
    }

    fn consume(&mut self, amt: usize) {
        self.reader.consume(amt);
        self.current_pos += amt as u64;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_byte_range_reader_single_range() {
        let data = b"0123456789ABCDEFGHIJ";
        let cursor = Cursor::new(data);
        let ranges = vec![(5, 10)]; // Read "56789"
        let mut reader = ByteRangeReader::new(cursor, ranges);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).unwrap();
        assert_eq!(buf, b"56789");
    }

    #[test]
    fn test_byte_range_reader_multiple_ranges() {
        let data = b"0123456789ABCDEFGHIJ";
        let cursor = Cursor::new(data);
        let ranges = vec![
            (0, 3),   // "012"
            (5, 8),   // "567"
            (15, 20), // "FGHIJ"
        ];
        let mut reader = ByteRangeReader::new(cursor, ranges);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).unwrap();
        assert_eq!(buf, b"012567FGHIJ");
    }

    #[test]
    fn test_byte_range_reader_empty_ranges() {
        let data = b"0123456789";
        let cursor = Cursor::new(data);
        let ranges = vec![];
        let mut reader = ByteRangeReader::new(cursor, ranges);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).unwrap();
        assert_eq!(buf, b"");
    }

    #[test]
    fn test_byte_range_reader_zero_length_range() {
        let data = b"0123456789";
        let cursor = Cursor::new(data);
        let ranges = vec![(5, 5)]; // Zero-length range
        let mut reader = ByteRangeReader::new(cursor, ranges);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).unwrap();
        assert_eq!(buf, b"");
    }

    #[test]
    fn test_byte_range_reader_full_file() {
        let data = b"0123456789";
        let cursor = Cursor::new(data);
        let ranges = vec![(0, 10)]; // Entire file
        let mut reader = ByteRangeReader::new(cursor, ranges);

        let mut buf = Vec::new();
        reader.read_to_end(&mut buf).unwrap();
        assert_eq!(buf, b"0123456789");
    }

    #[test]
    fn test_byte_range_reader_with_bufread() {
        let data = b"0123456789ABCDEFGHIJ";
        let cursor = Cursor::new(data);
        let ranges = vec![(2, 7), (10, 15)]; // "23456" + "ABCDE"
        let mut reader = ByteRangeReader::new(cursor, ranges);

        // Use BufRead methods
        let mut line = Vec::new();
        reader.read_until(b'E', &mut line).unwrap();
        assert_eq!(line, b"23456ABCDE");
    }
}

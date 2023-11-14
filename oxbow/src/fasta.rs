use arrow::{
    array::{ArrayRef, GenericStringBuilder},
    error::ArrowError,
    record_batch::RecordBatch,
};
use noodles::core::Region;
use noodles::fasta::{self, fai};
use std::{
    fs::File,
    io::{self, BufReader, Read, Seek},
    iter,
    path::Path,
    str,
    sync::Arc,
};

use crate::batch_builder::{write_ipc, BatchBuilder};

pub fn index_from_reader<R>(read: R) -> io::Result<fai::Index>
where
    R: Read,
{
    let mut fai_reader = fai::Reader::new(BufReader::new(read));
    fai_reader.read_index()
}

pub fn index_from_path(path: &str) -> io::Result<fai::Index> {
    let fai_path = format!("{}.fai", path);
    let index = if Path::new(&fai_path).exists() {
        fai::read(fai_path)?
    } else {
        panic!("Could not find a .fai index file for the given fasta file.");
    };
    Ok(index)
}

/// A FASTA reader.
pub struct FastaReader<R> {
    reader: fasta::Reader<BufReader<R>>,
    index: fai::Index,
}

impl FastaReader<BufReader<File>> {
    /// Creates a Fasta reader from a given file path.
    pub fn new_from_path(path: &str) -> io::Result<Self> {
        let index = index_from_path(path)?;
        let file = File::open(path)?;
        let bufreader = BufReader::with_capacity(1024 * 1024, file);
        let reader = fasta::Reader::new(BufReader::new(bufreader));
        Ok(Self { reader, index })
    }
}

impl<R: Read + Seek> FastaReader<R> {
    /// Creates a Fasta Reader.
    pub fn new(read: R, index: fai::Index) -> io::Result<Self> {
        let reader = fasta::Reader::new(BufReader::new(read));
        Ok(Self { reader, index })
    }

    /// Returns the records in the given region as Apache Arrow IPC.
    ///
    /// If the region is `None`, all records are returned.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use oxbow::fasta::FastaReader;
    ///
    /// let mut reader = FastaReader::new_from_path("sample.fasta.gz").unwrap();
    /// let ipc = reader.records_to_ipc(Some("sq0")).unwrap();
    /// ```
    pub fn records_to_ipc(&mut self, region: Option<&str>) -> Result<Vec<u8>, ArrowError> {
        let batch_builder = FastaBatchBuilder::new(1024)?;
        if let Some(region) = region {
            let region: Region = region.parse().unwrap();
            let query = self.reader.query(&self.index, &region).unwrap();
            let iter = iter::once(query);
            return write_ipc(iter, batch_builder);
        }
        let records = self.reader.records().map(|r| r.unwrap());
        write_ipc(records, batch_builder)
    }
}

struct FastaBatchBuilder {
    name: GenericStringBuilder<i32>,
    sequence: GenericStringBuilder<i32>,
}

impl FastaBatchBuilder {
    pub fn new(_capacity: usize) -> Result<Self, ArrowError> {
        Ok(Self {
            name: GenericStringBuilder::<i32>::new(),
            sequence: GenericStringBuilder::<i32>::new(),
        })
    }
}

impl BatchBuilder for FastaBatchBuilder {
    type Record<'a> = &'a fasta::record::Record;

    fn push(&mut self, record: Self::Record<'_>) {
        let seq = record.sequence().as_ref();

        self.name.append_value(record.name());
        self.sequence.append_value(str::from_utf8(seq).unwrap());
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        RecordBatch::try_from_iter(vec![
            ("name", Arc::new(self.name.finish()) as ArrayRef),
            ("sequence", Arc::new(self.sequence.finish()) as ArrayRef),
        ])
    }
}

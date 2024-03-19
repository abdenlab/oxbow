use arrow::{
    array::{ArrayRef, GenericStringBuilder},
    error::ArrowError,
    record_batch::RecordBatch,
};
use noodles::fasta;
use std::{
    io::{self, BufReader},
    iter, str,
    sync::Arc,
};

use crate::batch_builder::{write_ipc, BatchBuilder};

pub fn new_from_path(
    path: &str,
) -> io::Result<fasta::IndexedReader<Box<dyn fasta::io::BufReadSeek>>> {
    // Also reads the index file and handles (b)gzipped files
    fasta::indexed_reader::Builder::default().build_from_path(path)
}

pub fn new_from_reader<R>(fasta: R, fai: R) -> io::Result<fasta::IndexedReader<BufReader<R>>>
where
    R: io::Read,
{
    fasta::indexed_reader::Builder::default()
        .set_index(fasta::fai::Reader::new(BufReader::new(fai)).read_index()?)
        .build_from_reader(BufReader::new(fasta))
}

/// Returns the records in the given region as Apache Arrow IPC.
///
/// If the region is `None`, all records are returned.
///
/// # Examples
///
/// ```no_run
/// use oxbow::fasta;
///
/// let mut reader = fasta::new_from_path("sample.fasta.gz").unwrap();
/// let ipc = fasta::records_to_ipc(reader, Some("sq0")).unwrap();
/// ```
pub fn records_to_ipc<R>(
    mut indexed_reader: fasta::IndexedReader<R>,
    region: Option<&str>,
) -> Result<Vec<u8>, ArrowError>
where
    R: fasta::io::BufReadSeek,
{
    let batch_builder = FastaBatchBuilder::new(1024)?;
    if let Some(region) = region {
        let region = region.parse().unwrap();
        let query = indexed_reader.query(&region)?;
        let record_iter = iter::once(query);
        return write_ipc(record_iter, batch_builder);
    } else {
        let mut reader = fasta::reader::Builder.build_from_reader(indexed_reader.into_inner())?;
        let records = reader.records().map(|r| r.unwrap());
        return write_ipc(records, batch_builder);
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

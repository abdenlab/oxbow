use arrow::array::{
    ArrayRef, GenericStringBuilder
};
use arrow::{
    error::ArrowError, record_batch::RecordBatch,
};
use noodles::core::Region;
use noodles::fasta;
use noodles::fasta::fai;
use std::sync::Arc;

use crate::batch_builder::{write_ipc, BatchBuilder};

type BufferedReader = std::io::BufReader<std::fs::File>;

/// A FASTA reader.
pub struct FastaReader {
    reader: fasta::IndexedReader<BufferedReader>,
    stream_reader: fasta::Reader<Box<dyn std::io::BufRead>>,
}

impl FastaReader {
    /// Creates a Fasta Reader.
    pub fn new(path: &str) -> std::io::Result<Self> {
        let index = fai::read(format!("{}.fai", path))?;
        let file = std::fs::File::open(path)?;
        let bufreader = std::io::BufReader::with_capacity(1024 * 1024, file);
        let reader = fasta::indexed_reader::Builder::default()
            .set_index(index)
            .build_from_reader(bufreader)?;
        let stream_reader = fasta::reader::Builder::default()
            .build_from_path(path)?;
        Ok(Self { reader, stream_reader })
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
    /// let mut reader = FastaReader::new("sample.fasta.gz").unwrap();
    /// let ipc = reader.records_to_ipc(Some("sq0")).unwrap();
    /// ```
    pub fn records_to_ipc(&mut self, region: Option<&str>) -> Result<Vec<u8>, ArrowError> {
        let batch_builder = FastaBatchBuilder::new(1024)?;

        if let Some(region) = region {
            let region: Region = region.parse().unwrap();
            let query = self
                .reader
                .query(&region)
                .unwrap();
            let iter = std::iter::once(query);
            return write_ipc(iter, batch_builder);
        }

        let records = self.stream_reader.records().map(|r| r.unwrap());
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
    type Record = fasta::record::Record;
    

    fn push(&mut self, record: &Self::Record) {
        let seq = record.sequence().as_ref();
        
        self.name.append_value(record.name());
        self.sequence.append_value(std::str::from_utf8(seq).unwrap());
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        RecordBatch::try_from_iter(vec![
            ("name", Arc::new(self.name.finish()) as ArrayRef),
            ("sequence", Arc::new(self.sequence.finish()) as ArrayRef),
        ])
    }
}

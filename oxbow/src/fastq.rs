use arrow::array::{
    ArrayRef, GenericStringBuilder
};
use arrow::{
    error::ArrowError, record_batch::RecordBatch,
};
// use noodles::core::Region;
use noodles::fastq;
// use noodles::fastq::fai;
use std::sync::Arc;

use crate::batch_builder::{write_ipc, BatchBuilder};

type BufferedReader = std::io::BufReader<std::fs::File>;

/// A FASTQ reader.
pub struct FastqReader {
    reader: fastq::Reader<BufferedReader>,
    // index: fai::Reader<std::io::BufReader<std::fs::File>>,
}

impl FastqReader {
    /// Creates a Fastq Reader.
    pub fn new(path: &str) -> std::io::Result<Self> {
        // let index_file = std::fs::File::open(format!("{}.fai", path))?;
        // let index_bufreader = std::io::BufReader::with_capacity(1024 * 1024, index_file);
        // let index = fai::Reader::new(index_bufreader);

        let file = std::fs::File::open(path)?;
        let bufreader = std::io::BufReader::with_capacity(1024 * 1024, file);
        let reader = fastq::reader::Reader::new(bufreader);

        Ok(Self { reader })
    }

    /// Returns the records in the given region as Apache Arrow IPC.
    ///
    /// If the region is `None`, all records are returned.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use oxbow::fastq::FastqReader;
    ///
    /// let mut reader = FastqReader::new("sample.fastq.gz").unwrap();
    /// let ipc = reader.records_to_ipc().unwrap();
    /// ```
    pub fn records_to_ipc(&mut self) -> Result<Vec<u8>, ArrowError> {
        let batch_builder = FastqBatchBuilder::new(1024)?;
        let records = self.reader.records().map(|r| r.unwrap());
        write_ipc(records, batch_builder)
    }
}

struct FastqBatchBuilder {
    name: GenericStringBuilder<i32>,
    description: GenericStringBuilder<i32>,
    sequence: GenericStringBuilder<i32>,
    quality_scores: GenericStringBuilder<i32>,
}

impl FastqBatchBuilder {
    pub fn new(_capacity: usize) -> Result<Self, ArrowError> {
        Ok(Self {
            name: GenericStringBuilder::<i32>::new(),
            description: GenericStringBuilder::<i32>::new(),
            sequence: GenericStringBuilder::<i32>::new(),
            quality_scores: GenericStringBuilder::<i32>::new(),
        })
    }
}

impl BatchBuilder for FastqBatchBuilder {
    type Record = fastq::Record;
    
    fn push(&mut self, record: &Self::Record) {
        self.name.append_value(std::str::from_utf8(record.name()).unwrap());
        self.description.append_value(std::str::from_utf8(record.description()).unwrap());
        self.sequence.append_value(std::str::from_utf8(record.sequence()).unwrap());
        self.quality_scores.append_value(std::str::from_utf8(record.quality_scores()).unwrap());
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        RecordBatch::try_from_iter(vec![
            ("name", Arc::new(self.name.finish()) as ArrayRef),
            ("description", Arc::new(self.description.finish()) as ArrayRef),
            ("sequence", Arc::new(self.sequence.finish()) as ArrayRef),
            ("quality_scores", Arc::new(self.quality_scores.finish()) as ArrayRef),
        ])
    }
}

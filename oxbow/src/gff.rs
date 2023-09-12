use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::sync::Arc;

use arrow::array::{ArrayRef, Float32Builder, GenericStringBuilder, Int32Builder};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use noodles::gff;

use crate::batch_builder::{write_ipc_err, BatchBuilder};

/// A GFF reader.
pub struct GffReader<R> {
    reader: gff::Reader<BufReader<R>>,
}

impl GffReader<BufReader<File>> {
    /// Creates a GFF reader from a given file path.
    pub fn new_from_path(path: &str) -> std::io::Result<Self> {
        let reader = File::open(path)
            .map(BufReader::new)
            .map(BufReader::new)
            .map(gff::Reader::new)?;
        Ok(Self { reader })
    }
}

impl<R: Read + Seek> GffReader<R> {
    /// Creates a GFF Reader.
    pub fn new(read: R) -> std::io::Result<Self> {
        let reader = gff::Reader::new(BufReader::new(read));
        Ok(Self { reader })
    }

    /// Returns the records in the given region as Apache Arrow IPC.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use oxbow::gff::GffReader;
    ///
    /// let mut reader = GffReader::new_from_path("sample.gff").unwrap();
    /// let ipc = reader.records_to_ipc().unwrap();
    /// ```
    pub fn records_to_ipc(&mut self) -> Result<Vec<u8>, ArrowError> {
        let batch_builder = GffBatchBuilder::new(1024)?;
        let records = self
            .reader
            .records()
            .map(|i| i.map_err(|e| ArrowError::ExternalError(e.into())));
        return write_ipc_err(records, batch_builder);
    }
}

struct GffBatchBuilder {
    reference_sequence_name: GenericStringBuilder<i32>,
    source: GenericStringBuilder<i32>,
    ty: GenericStringBuilder<i32>,
    start: Int32Builder,
    end: Int32Builder,
    score: Float32Builder,
    strand: GenericStringBuilder<i32>,
    phase: GenericStringBuilder<i32>,
    attributes: GenericStringBuilder<i32>,
}

impl GffBatchBuilder {
    pub fn new(capacity: usize) -> Result<Self, ArrowError> {
        Ok(Self {
            reference_sequence_name: GenericStringBuilder::<i32>::new(),
            source: GenericStringBuilder::<i32>::new(),
            ty: GenericStringBuilder::<i32>::new(),
            start: Int32Builder::with_capacity(capacity),
            end: Int32Builder::with_capacity(capacity),
            score: Float32Builder::new(),
            strand: GenericStringBuilder::<i32>::new(),
            phase: GenericStringBuilder::<i32>::new(),
            attributes: GenericStringBuilder::<i32>::new(),
        })
    }
}

impl BatchBuilder for GffBatchBuilder {
    type Record<'a> = &'a gff::record::Record;

    fn push(&mut self, record: Self::Record<'_>) {
        self.reference_sequence_name
            .append_value(record.reference_sequence_name().to_string());
        self.source.append_value(record.source().to_string());
        self.ty.append_value(record.ty().to_string());
        self.start.append_value(usize::from(record.start()) as i32);
        self.end.append_value(usize::from(record.end()) as i32);
        match record.score() {
            Some(score) => self.score.append_value(score),
            None => self.score.append_null(),
        }
        self.strand.append_value(record.strand().to_string());
        match record.phase() {
            Some(phase) => self.phase.append_value(phase),
            None => self.phase.append_null(),
        }
        self.attributes
            .append_value(record.attributes().to_string());
    }

    fn finish(mut self) -> Result<RecordBatch, ArrowError> {
        RecordBatch::try_from_iter(vec![
            (
                "seqid",
                Arc::new(self.reference_sequence_name.finish()) as ArrayRef,
            ),
            ("source", Arc::new(self.source.finish()) as ArrayRef),
            ("type", Arc::new(self.ty.finish()) as ArrayRef),
            ("start", Arc::new(self.start.finish()) as ArrayRef),
            ("end", Arc::new(self.end.finish()) as ArrayRef),
            ("score", Arc::new(self.score.finish()) as ArrayRef),
            ("strand", Arc::new(self.strand.finish()) as ArrayRef),
            ("phase", Arc::new(self.phase.finish()) as ArrayRef),
            ("attributes", Arc::new(self.attributes.finish()) as ArrayRef),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::ipc::reader::FileReader;
    use arrow::record_batch::RecordBatch;

    fn read_record_batch() -> RecordBatch {
        let mut dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        dir.push("../fixtures/example.gff");
        let mut reader = GffReader::new_from_path(dir.to_str().unwrap()).unwrap();
        let ipc = reader.records_to_ipc().unwrap();
        let cursor = std::io::Cursor::new(ipc);
        let mut arrow_reader = FileReader::try_new(cursor, None).unwrap();
        // make sure we have one batch
        assert_eq!(arrow_reader.num_batches(), 1);
        arrow_reader.next().unwrap().unwrap()
    }

    #[test]
    fn test_read_all() {
        let record_batch = read_record_batch();
        assert_eq!(record_batch.num_rows(), 6);
    }
}

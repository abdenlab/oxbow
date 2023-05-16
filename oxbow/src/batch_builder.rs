use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;

pub trait BatchBuilder {
    type Record;
    fn push(&mut self, record: &Self::Record);
    fn finish(self) -> Result<RecordBatch, ArrowError>;
}

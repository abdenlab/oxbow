use std::io;

use arrow::datatypes::SchemaRef;
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;

/// A builder for Arrow record batches.
pub trait RecordBatchBuilder {
    /// Returns the Arrow schema for the record batches produced by this builder.
    fn schema(&self) -> SchemaRef;

    /// Finalizes the current batch and returns it, resetting the builder for the next batch.
    fn finish(&mut self) -> Result<RecordBatch, ArrowError>;
}

/// Push a record into a batch builder.
pub trait Push<T> {
    fn push(&mut self, record: T) -> io::Result<()>;
}

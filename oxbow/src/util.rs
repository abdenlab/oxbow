use arrow::array::RecordBatchReader;
use arrow::error::ArrowError;
use arrow::ipc::writer::FileWriter;

pub mod index;
pub mod query;

/// Serializes a sequence of Arrow record batches into Arrow IPC bytes.
pub fn batches_to_ipc(batches: impl RecordBatchReader) -> Result<Vec<u8>, ArrowError> {
    let schema = batches.schema();
    let mut writer = FileWriter::try_new(Vec::new(), &schema)?;
    for batch in batches {
        writer.write(&batch?)?;
    }
    writer.finish()?;
    writer.into_inner()
}

use arrow::error::ArrowError;
use arrow::ipc::writer::FileWriter;
use arrow::record_batch::RecordBatch;

pub trait BatchBuilder {
    type Record<'a>;
    fn push(&mut self, record: Self::Record<'_>);
    fn finish(self) -> Result<RecordBatch, ArrowError>;
}

pub fn write_ipc_err<T>(
    records: impl Iterator<Item = Result<T, ArrowError>>,
    mut batch_builder: impl for<'a> BatchBuilder<Record<'a> = &'a T>,
) -> Result<Vec<u8>, ArrowError> {
    for record in records {
        let record = record?;
        batch_builder.push(&record)
    }
    finish_batch(batch_builder)
}

pub fn write_ipc<T>(
    records: impl Iterator<Item = T>,
    mut batch_builder: impl for<'a> BatchBuilder<Record<'a> = &'a T>,
) -> Result<Vec<u8>, ArrowError> {
    for record in records {
        let record = record;
        batch_builder.push(&record)
    }
    finish_batch(batch_builder)
}

pub fn finish_batch(batch_builder: impl BatchBuilder) -> Result<Vec<u8>, ArrowError> {
    let batch = batch_builder.finish()?;
    let mut writer = FileWriter::try_new(Vec::new(), &batch.schema())?;
    writer.write(&batch)?;
    writer.finish()?;
    writer.into_inner()
}

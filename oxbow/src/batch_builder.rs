use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use arrow::ipc::writer::FileWriter;

pub trait BatchBuilder {
    type Record<'a>;
    fn push(&mut self, record: Self::Record<'_>);
    fn finish(self) -> Result<RecordBatch, ArrowError>;
}

pub fn write_ipc<T>(
    records: impl Iterator<Item = T>,
    mut batch_builder: impl for<'a> BatchBuilder<Record<'a> = &'a T>,
) -> Result<Vec<u8>, ArrowError> {
    records.for_each(|record| batch_builder.push(&record));
    let batch = batch_builder.finish()?;
    let mut writer = FileWriter::try_new(Vec::new(), &batch.schema())?;
    writer.write(&batch)?;
    writer.finish()?;
    writer.into_inner()
}

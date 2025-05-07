use arrow::array::{DictionaryArray, RecordBatchReader, StringArray, StringDictionaryBuilder};
use arrow::datatypes::ArrowDictionaryKeyType;
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

// Build a DictionaryArray and reinitialize the builder with the same dictionary values.
pub(crate) fn reset_dictarray_builder<T: ArrowDictionaryKeyType>(
    builder: &mut StringDictionaryBuilder<T>,
) -> DictionaryArray<T> {
    let array = builder.finish();
    let dict_values = array
        .values()
        .as_any()
        .downcast_ref::<StringArray>()
        .expect("Failed to downcast to StringArray")
        .clone();
    *builder =
        StringDictionaryBuilder::<T>::new_with_dictionary(array.len(), &dict_values).unwrap();
    array
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{ArrayRef, DictionaryArray, Int32Array, StringArray};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    #[test]
    fn test_batches_to_ipc_with_dictionary() {
        // Create a schema and record batch
        let schema = Arc::new(Schema::new(vec![Field::new(
            "dict_col",
            DataType::Dictionary(Box::new(DataType::Int32), Box::new(DataType::Utf8)),
            false,
        )]));

        let array1: ArrayRef = {
            let keys = Int32Array::from(vec![Some(0), Some(1), Some(0), Some(2)]);
            let values = Arc::new(StringArray::from(vec!["a", "b", "c"])) as ArrayRef;
            let array = DictionaryArray::try_new(keys, values).unwrap();
            Arc::new(array)
        };
        let batch1 = RecordBatch::try_new(schema.clone(), vec![array1.clone()]).unwrap();

        let array2: ArrayRef = {
            let keys = Int32Array::from(vec![Some(0), Some(0), Some(0), Some(2)]);
            let values = Arc::new(StringArray::from(vec!["a", "b", "c"])) as ArrayRef;
            let array = DictionaryArray::try_new(keys, values).unwrap();
            Arc::new(array)
        };
        let batch2 = RecordBatch::try_new(schema.clone(), vec![array2.clone()]).unwrap();

        // Serialize the record batch to IPC
        let batches = vec![batch1, batch2];
        let reader = arrow::record_batch::RecordBatchIterator::new(
            batches.into_iter().map(Ok),
            schema.clone(),
        );
        let ipc_bytes = batches_to_ipc(reader).expect("Serialization failed");

        // Deserialize the IPC bytes back into record batches
        let cursor = std::io::Cursor::new(ipc_bytes);
        let reader =
            arrow::ipc::reader::FileReader::try_new(cursor, None).expect("Deserialization failed");

        // Verify the deserialized data matches the original
        let deserialized_batches: Vec<RecordBatch> = reader.collect::<Result<_, _>>().unwrap();
        assert_eq!(deserialized_batches.len(), 2);
        assert_eq!(deserialized_batches[0].schema(), schema);
        assert_eq!(deserialized_batches[0].column(0).as_ref(), array1.as_ref());
    }
}

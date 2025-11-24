use std::{panic, panic::AssertUnwindSafe};

use arrow::error::ArrowError;
use arrow::record_batch::RecordBatchReader;

/// An Arrow RecordBatchReader wrapper that converts panics to errors during iteration.
///
/// Once a panic is caught, the reader becomes unusable and will return `None` on all
/// subsequent calls to `next()`, as the inner reader's state may be corrupted.
pub(crate) struct UnwindCatchingRecordBatchReader<R>
where
    R: RecordBatchReader,
{
    inner: R,
    panicked: bool,
}

impl<R> UnwindCatchingRecordBatchReader<R>
where
    R: RecordBatchReader,
{
    /// Creates a new UnwindCatchingRecordBatchReader wrapping the given reader.
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            panicked: false,
        }
    }
}

impl<R> RecordBatchReader for UnwindCatchingRecordBatchReader<R>
where
    R: RecordBatchReader,
{
    fn schema(&self) -> arrow::datatypes::SchemaRef {
        self.inner.schema()
    }
}

impl<R> Iterator for UnwindCatchingRecordBatchReader<R>
where
    R: RecordBatchReader,
{
    type Item = Result<arrow::record_batch::RecordBatch, ArrowError>;

    fn next(&mut self) -> Option<Self::Item> {
        // If we've already caught a panic, the reader is poisoned and unusable
        if self.panicked {
            return None;
        }

        match panic::catch_unwind(AssertUnwindSafe(|| self.inner.next())) {
            Ok(opt) => opt,
            Err(payload) => {
                self.panicked = true;
                let message = if let Some(s) = payload.downcast_ref::<&str>() {
                    format!("Panic in RecordBatchReader: {}", s)
                } else if let Some(s) = payload.downcast_ref::<String>() {
                    format!("Panic in RecordBatchReader: {}", s)
                } else {
                    "Panic in RecordBatchReader (non-string payload)".to_string()
                };
                Some(Err(ArrowError::ExternalError(Box::new(
                    std::io::Error::other(message),
                ))))
            }
        }
    }
}

/// Wraps a RecordBatchReader to catch unwinding panics during iteration and convert them to errors.
///
/// Use this function to wrap a [`RecordBatchReader`] before passing it to Python via
/// [`pyo3_arrow::PyRecordBatchReader`] to make it safe to use across other native extensions.
///
/// # Purpose
///
/// Generally, if Rust is called from foreign code and a panic occurs that unwinds across an FFI
/// boundary, the process it is embedded in will abort. For our [`RecordBatchReader`]
/// implementations, which return Arrow record batches via iteration, this is undesirable behavior:
/// for example, a panic during iteration in a Jupyter environment would kill the Python kernel
/// process. Therefore, we would like to propagate panics in a defined way to the FFI caller.
///
/// We export Oxbow's RecordBatchReaders to Python using the [`pyo3_arrow::PyRecordBatchReader`]
/// wrapper. If we only cared about the Rust-Python boundary, we wouldn't need to worry because PyO3
/// already catches panics and converts them into Python exceptions [`pyo3::panic::PanicException`].
///
/// However, [`pyo3_arrow::PyRecordBatchReader`] also implements the
/// [Arrow C Stream Interface](https://arrow.apache.org/docs/format/CStreamInterface.html), which
/// provides a "capsule" object that can be consumed by code within other Python extension modules,
/// such as `pyarrow`, written in other languages (e.g., C/C++). If a panic were to occur while such
/// code reads record batches through the Arrow C Stream API, it would unwind across a different FFI
/// boundary without getting caught, which again would abort the process.
///
/// To prevent this from happening, the returned wrapper uses [`std::panic::catch_unwind`] to
/// intercept any unwinding panics that occur while reading record batches, converting them into
/// [`arrow::error::ArrowError`], as required by the [`RecordBatchReader`] trait. A consumer of the
/// RecordBatchReader through the Arrow C Stream Interface will be able to receive information
/// about the error and handle or propagate it appropriately, while a Python caller using the
/// exported PyRecordBatchReader directly would still receive a Python exception, as desired;
/// however, rather than pyo3's `PanicException`, the exception raised would be derived from
/// `ArrowError`. The iterator is poisoned after a panic so that subsequent calls to `next()` will
/// return `None`.
///
/// # Arguments
///
/// * `reader` - The RecordBatchReader to wrap
///
/// # Returns
///
/// A boxed [`UnwindCatchingRecordBatchReader`]
///
/// # Example
///
/// ```ignore
/// let reader = scanner.scan(file, None, None, None, None)?;
/// let safe_reader = err_on_unwind(reader);
/// let py_reader = PyRecordBatchReader::new(safe_reader);
/// ```
pub(crate) fn err_on_unwind<R>(reader: R) -> Box<UnwindCatchingRecordBatchReader<R>>
where
    R: RecordBatchReader,
{
    Box::new(UnwindCatchingRecordBatchReader::new(reader))
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{Int32Array, RecordBatch};
    use arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    /// A mock RecordBatchReader that yields a fixed number of batches without panicking
    struct MockReader {
        schema: arrow::datatypes::SchemaRef,
        remaining: usize,
    }

    impl MockReader {
        fn new(num_batches: usize) -> Self {
            let schema = Arc::new(Schema::new(vec![Field::new("a", DataType::Int32, false)]));
            Self {
                schema,
                remaining: num_batches,
            }
        }
    }

    impl Iterator for MockReader {
        type Item = Result<RecordBatch, ArrowError>;

        fn next(&mut self) -> Option<Self::Item> {
            if self.remaining == 0 {
                return None;
            }
            self.remaining -= 1;
            let array = Int32Array::from(vec![1, 2, 3]);
            let batch = RecordBatch::try_new(self.schema.clone(), vec![Arc::new(array)]).unwrap();
            Some(Ok(batch))
        }
    }

    impl RecordBatchReader for MockReader {
        fn schema(&self) -> arrow::datatypes::SchemaRef {
            self.schema.clone()
        }
    }

    /// A mock RecordBatchReader that panics on the first call to next()
    struct PanickingReader {
        schema: arrow::datatypes::SchemaRef,
        panic_message: String,
    }

    impl PanickingReader {
        fn new(panic_message: String) -> Self {
            let schema = Arc::new(Schema::new(vec![Field::new("a", DataType::Int32, false)]));
            Self {
                schema,
                panic_message,
            }
        }
    }

    impl Iterator for PanickingReader {
        type Item = Result<RecordBatch, ArrowError>;

        fn next(&mut self) -> Option<Self::Item> {
            panic!("{}", self.panic_message);
        }
    }

    impl RecordBatchReader for PanickingReader {
        fn schema(&self) -> arrow::datatypes::SchemaRef {
            self.schema.clone()
        }
    }

    #[test]
    fn test_err_on_unwind_normal_reader() {
        // Test that a normal reader works through the wrapper
        let reader = MockReader::new(3);
        let mut safe_reader = err_on_unwind(reader);

        // Should be able to read all 3 batches
        assert!(safe_reader.next().is_some());
        assert!(safe_reader.next().is_some());
        assert!(safe_reader.next().is_some());
        assert!(safe_reader.next().is_none());
    }

    #[test]
    fn test_err_on_unwind_catches_panic() {
        // Test that a panicking reader has its panic caught and converted to an error
        let reader = PanickingReader::new("test panic message".to_string());
        let mut safe_reader = err_on_unwind(reader);

        // The first call should return Some(Err(...)) instead of panicking
        let result = safe_reader.next();
        assert!(result.is_some());

        let err = result.unwrap();
        assert!(err.is_err());

        // Check that the error message contains our panic message
        let err_msg = err.unwrap_err().to_string();
        assert!(err_msg.contains("test panic message"));
    }

    #[test]
    fn test_reader_unusable_after_panic() {
        // Test that after catching a panic, the reader returns None on subsequent calls
        let reader = PanickingReader::new("initial panic".to_string());
        let mut safe_reader = err_on_unwind(reader);

        // First call returns the error
        let result = safe_reader.next();
        assert!(result.is_some());
        assert!(result.unwrap().is_err());

        // All subsequent calls should return None (reader is poisoned)
        assert!(safe_reader.next().is_none());
        assert!(safe_reader.next().is_none());
        assert!(safe_reader.next().is_none());
    }
}

use std::{future::Future, pin::Pin};

use arrow_array::RecordBatch;
use futures::Stream;

pub type AsyncBatchReader = Pin<Box<dyn Stream<Item = RecordBatch>>>;

pub trait AsyncScanner {
    fn scan(
        &self,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> Pin<Box<dyn Future<Output = AsyncBatchReader> + Send>>;
}

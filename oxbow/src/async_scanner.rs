use std::pin::Pin;

use crate::Result;
use arrow_array::RecordBatch;
use futures::Stream;
use tokio::io::AsyncBufRead;

pub trait AsyncScanner {
    fn scan<R: AsyncBufRead + Unpin + Send + 'static>(
        &self,
        reader: R,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> Pin<Box<dyn Stream<Item = Result<RecordBatch>>>>;
}

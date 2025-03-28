mod query;
mod range;
mod stream;

pub use query::BatchIterator as QueryBatchIterator;
pub use range::BatchIterator as RangeBatchIterator;
pub use stream::BatchIterator;

pub mod batch_iterator;
pub mod format;
pub mod model;

pub use format::bam::Scanner as BamScanner;
pub use format::cram::Scanner as CramScanner;
pub use format::sam::Scanner as SamScanner;

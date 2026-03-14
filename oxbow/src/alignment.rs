pub mod model;
pub mod scanner;

pub use model::Model as AlignmentModel;
pub use scanner::bam::Scanner as BamScanner;
pub use scanner::cram::Scanner as CramScanner;
pub use scanner::sam::Scanner as SamScanner;

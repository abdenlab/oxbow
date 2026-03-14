pub mod model;
pub mod scanner;

pub use model::Model as SequenceModel;
pub use scanner::fasta::Scanner as FastaScanner;
pub use scanner::fastq::Scanner as FastqScanner;

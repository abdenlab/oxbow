pub mod batch_iterator;
pub mod format;
pub mod model;

pub use format::bcf::Scanner as BcfScanner;
pub use format::vcf::Scanner as VcfScanner;
pub use model::GenotypeBy;

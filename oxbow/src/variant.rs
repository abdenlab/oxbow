pub mod model;
pub mod scanner;

pub use model::GenotypeBy;
pub use model::Model as VariantModel;
pub use scanner::bcf::Scanner as BcfScanner;
pub use scanner::vcf::Scanner as VcfScanner;

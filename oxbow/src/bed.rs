pub mod model;
pub mod scanner;

pub use model::{BedSchema, FieldDef, FieldType};
pub use scanner::bed::Scanner as BedScanner;

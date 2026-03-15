pub mod model;
pub mod scanner;

pub use crate::bed::model::schema::BedSchema;
pub use model::base::Model as BBIBaseModel;
pub use model::zoom::Model as BBIZoomModel;
pub use scanner::bbizoom::Scanner as BBIZoomScanner;
pub use scanner::bigbed::Scanner as BigBedScanner;
pub use scanner::bigwig::Scanner as BigWigScanner;
pub use scanner::BBIReader;

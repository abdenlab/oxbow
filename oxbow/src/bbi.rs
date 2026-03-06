pub mod model;
pub mod scanner;

pub use model::base::BedSchema;
pub use scanner::bbizoom::Scanner as BBIZoomScanner;
pub use scanner::bigbed::Scanner as BigBedScanner;
pub use scanner::bigwig::Scanner as BigWigScanner;
pub use scanner::BBIReader;

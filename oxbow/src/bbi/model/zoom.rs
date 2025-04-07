pub mod batch_builder;
pub mod field;

pub use batch_builder::BatchBuilder;
pub use batch_builder::Push;

pub struct BBIZoomRecord<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
    pub bases_covered: u64,
    pub min: f64,
    pub max: f64,
    pub sum: f64,
    pub sum_squares: f64,
}

impl<'a> BBIZoomRecord<'a> {
    pub fn new(chrom: &'a str, entry: &'a bigtools::ZoomRecord) -> Self {
        Self {
            chrom,
            start: entry.start,
            end: entry.end,
            bases_covered: entry.summary.bases_covered,
            min: entry.summary.min_val,
            max: entry.summary.max_val,
            sum: entry.summary.sum,
            sum_squares: entry.summary.sum_squares,
        }
    }
}

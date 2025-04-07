pub mod batch_builder;
pub mod field;
pub mod schema;

pub use batch_builder::BatchBuilder;
pub use batch_builder::Push;
pub use schema::BedSchema;

pub struct BigBedRecord<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
    pub rest: &'a String,
}

impl<'a> BigBedRecord<'a> {
    pub fn new(chrom: &'a str, entry: &'a bigtools::BedEntry) -> Self {
        Self {
            chrom,
            start: entry.start,
            end: entry.end,
            rest: &entry.rest,
        }
    }
}

pub struct BigWigRecord<'a> {
    pub chrom: &'a str,
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

impl<'a> BigWigRecord<'a> {
    pub fn new(chrom: &'a str, entry: &'a bigtools::Value) -> Self {
        Self {
            chrom,
            start: entry.start,
            end: entry.end,
            value: entry.value,
        }
    }
}

pub mod bbizoom;
pub mod bigbed;
pub mod bigwig;

use bigtools::{BigBedRead, BigWigRead};

pub enum BBIReader<R> {
    BigWig(BigWigRead<R>),
    BigBed(BigBedRead<R>),
}

pub enum BBIFileType {
    BigWig,
    BigBed,
}

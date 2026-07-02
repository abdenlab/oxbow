use std::future::Future;

use crate::variant::model::{BatchBuilder, GenotypeBy, Model};

pub struct AsyncVcfScanner {
    header: noodles::vcf::Header,
    model: Model,
}
impl AsyncVcfScanner {
    pub fn builder() -> AsyncVcfScannerBuilder {
        AsyncVcfScannerBuilder::default()
    }
}

#[derive(Default)]
pub struct AsyncVcfScannerBuilder {}

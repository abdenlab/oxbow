use pyo3::prelude::*;

mod filelike;
mod util;

mod sequence;
use crate::sequence::{read_fasta, read_fastq, PyFastaScanner, PyFastqScanner};

/////////////////
/// Python module.
/////////////////

#[pymodule]
#[pyo3(name = "oxbow")]
fn oxbow_sandbox(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFastaScanner>()?;
    m.add_class::<PyFastqScanner>()?;
    m.add_function(wrap_pyfunction!(read_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(read_fastq, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

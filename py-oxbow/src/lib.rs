use pyo3::prelude::*;

mod filelike;
mod util;

mod alignment;
mod sequence;
use crate::alignment::{read_bam, read_sam, PyBamScanner, PySamScanner};
use crate::sequence::{read_fasta, read_fastq, PyFastaScanner, PyFastqScanner};

/////////////////
/// Python module.
/////////////////

#[pymodule]
#[pyo3(name = "oxbow")]
fn oxbow_sandbox(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFastaScanner>()?;
    m.add_class::<PyFastqScanner>()?;
    m.add_class::<PySamScanner>()?;
    m.add_class::<PyBamScanner>()?;
    m.add_function(wrap_pyfunction!(read_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(read_fastq, m)?)?;
    m.add_function(wrap_pyfunction!(read_sam, m)?)?;
    m.add_function(wrap_pyfunction!(read_bam, m)?)?;
    Ok(())
}

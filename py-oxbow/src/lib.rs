use pyo3::prelude::*;

mod filelike;
mod util;

mod alignment;
mod bbi;
mod bed;
mod error;
mod gxf;
mod sequence;
mod variant;
use crate::alignment::{read_bam, read_cram, read_sam, PyBamScanner, PyCramScanner, PySamScanner};
use crate::bbi::{
    read_bigbed, read_bigwig, PyBBIFileType, PyBBIZoomScanner, PyBigBedScanner, PyBigWigScanner,
};
use crate::bed::{read_bed, PyBedScanner};
use crate::gxf::{read_gff, read_gtf, PyGffScanner, PyGtfScanner};
use crate::sequence::{read_fasta, read_fastq, PyFastaScanner, PyFastqScanner};
use crate::util::partition_from_index;
use crate::variant::{read_bcf, read_vcf, PyBcfScanner, PyVcfScanner};

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
    m.add_class::<PyCramScanner>()?;
    m.add_class::<PyVcfScanner>()?;
    m.add_class::<PyBcfScanner>()?;
    m.add_class::<PyGffScanner>()?;
    m.add_class::<PyGtfScanner>()?;
    m.add_class::<PyBedScanner>()?;
    m.add_class::<PyBigBedScanner>()?;
    m.add_class::<PyBigWigScanner>()?;
    m.add_class::<PyBBIFileType>()?;
    m.add_class::<PyBBIZoomScanner>()?;
    m.add_function(wrap_pyfunction!(read_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(read_fastq, m)?)?;
    m.add_function(wrap_pyfunction!(read_sam, m)?)?;
    m.add_function(wrap_pyfunction!(read_bam, m)?)?;
    m.add_function(wrap_pyfunction!(read_cram, m)?)?;
    m.add_function(wrap_pyfunction!(read_gff, m)?)?;
    m.add_function(wrap_pyfunction!(read_gtf, m)?)?;
    m.add_function(wrap_pyfunction!(read_vcf, m)?)?;
    m.add_function(wrap_pyfunction!(read_bcf, m)?)?;
    m.add_function(wrap_pyfunction!(read_bed, m)?)?;
    m.add_function(wrap_pyfunction!(read_bigbed, m)?)?;
    m.add_function(wrap_pyfunction!(read_bigwig, m)?)?;
    m.add_function(wrap_pyfunction!(partition_from_index, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

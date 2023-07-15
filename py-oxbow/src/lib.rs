use pyo3::prelude::*;
use pyo3::types::PyBytes;

use oxbow::fasta::FastaReader;
use oxbow::fastq::FastqReader;
use oxbow::bam::BamReader;
use oxbow::cram::CramReader;
use oxbow::vcf::VcfReader;
use oxbow::bcf::BcfReader;

#[pyfunction]
fn read_fasta(path: &str, region: Option<&str>) -> PyObject {
    let mut reader = FastaReader::new(path).unwrap();
    let ipc = reader.records_to_ipc(region).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pyfunction]
fn read_fastq(path: &str) -> PyObject {
    let mut reader = FastqReader::new(path).unwrap();
    let ipc = reader.records_to_ipc().unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pyfunction]
fn read_bam(path: &str, region: Option<&str>) -> PyObject {
    let mut reader = BamReader::new(path).unwrap();
    let ipc = reader.records_to_ipc(region).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pyfunction]
fn read_cram(path: &str, fasta_path: &str, region: Option<&str>) -> PyObject {
    let mut reader = CramReader::new(path, fasta_path).unwrap();
    let ipc = reader.records_to_ipc(region).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pyfunction]
fn read_vcf(path: &str, region: Option<&str>) -> PyObject {
    let mut reader = VcfReader::new(path).unwrap();
    let ipc = reader.records_to_ipc(region).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pyfunction]
fn read_bcf(path: &str, region: Option<&str>) -> PyObject {
    let mut reader = BcfReader::new(path).unwrap();
    let ipc = reader.records_to_ipc(region).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pymodule]
#[pyo3(name = "oxbow")]
fn py_oxbow(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(read_fastq, m)?)?;
    m.add_function(wrap_pyfunction!(read_bam, m)?)?;
    m.add_function(wrap_pyfunction!(read_cram, m)?)?;
    m.add_function(wrap_pyfunction!(read_vcf, m)?)?;
    m.add_function(wrap_pyfunction!(read_bcf, m)?)?;
    Ok(())
}

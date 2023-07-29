use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::types::PyList;

use oxbow::fasta::FastaReader;
use oxbow::fastq::FastqReader;
use oxbow::bam::BamReader;
// use oxbow::cram::CramReader;
use oxbow::vcf::VcfReader;
use oxbow::bcf::BcfReader;

use oxbow::vpos;


#[pyfunction]
fn partition_from_index_file(path: &str, chunksize: u64) -> PyObject {
    let voffsets = vpos::partition_from_index_file(path, chunksize);
    Python::with_gil(|py| PyList::new(py, &voffsets).into())
}

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
fn read_bam_vpos(path: &str, pos_lo: (u64, u16), pos_hi: (u64, u16)) -> PyObject {
    let mut reader = BamReader::new(path).unwrap();
    let ipc = reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

// #[pyfunction]
// fn read_cram(path: &str, fasta_path: &str, region: Option<&str>) -> PyObject {
//     let mut reader = CramReader::new(path, fasta_path).unwrap();
//     let ipc = reader.records_to_ipc(region).unwrap();
//     Python::with_gil(|py| PyBytes::new(py, &ipc).into())
// }

#[pyfunction]
fn read_vcf(path: &str, region: Option<&str>) -> PyObject {
    let mut reader = VcfReader::new(path).unwrap();
    let ipc = reader.records_to_ipc(region).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pyfunction]
fn read_vcf_vpos(path: &str, pos_lo: (u64, u16), pos_hi: (u64, u16)) -> PyObject {
    let mut reader = VcfReader::new(path).unwrap();
    let ipc = reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pyfunction]
fn read_bcf(path: &str, region: Option<&str>) -> PyObject {
    let mut reader = BcfReader::new(path).unwrap();
    let ipc = reader.records_to_ipc(region).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pyfunction]
fn read_bcf_vpos(path: &str, pos_lo: (u64, u16), pos_hi: (u64, u16)) -> PyObject {
    let mut reader = BcfReader::new(path).unwrap();
    let ipc = reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pymodule]
#[pyo3(name = "oxbow")]
fn py_oxbow(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(read_fastq, m)?)?;
    m.add_function(wrap_pyfunction!(partition_from_index_file, m)?)?;
    m.add_function(wrap_pyfunction!(read_bam, m)?)?;
    m.add_function(wrap_pyfunction!(read_bam_vpos, m)?)?;
    // m.add_function(wrap_pyfunction!(read_cram, m)?)?;
    // m.add_function(wrap_pyfunction!(read_cram_vpos, m)?)?;
    m.add_function(wrap_pyfunction!(read_vcf, m)?)?;
    m.add_function(wrap_pyfunction!(read_vcf_vpos, m)?)?;
    m.add_function(wrap_pyfunction!(read_bcf, m)?)?;
    m.add_function(wrap_pyfunction!(read_bcf_vpos, m)?)?;
    Ok(())
}

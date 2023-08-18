use pyo3::prelude::*;
use pyo3::types::PyBytes;
use pyo3::types::PyList;
use pyo3::types::PyString;

use oxbow::bam;
use oxbow::bam::BamReader;
use oxbow::bigbed::BigBedReader;
use oxbow::bigwig::BigWigReader;
use oxbow::fasta::FastaReader;
use oxbow::fastq::FastqReader;
// use oxbow::cram::CramReader;
use oxbow::bcf::BcfReader;
use oxbow::vcf::VcfReader;

use oxbow::vpos;

mod file_like;

use file_like::PyFileLikeObject;

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
<<<<<<< HEAD
fn read_bam(
    py: Python,
    path_or_file_like: PyObject,
    index: Option<PyObject>,
    region: Option<&str>,
) -> PyObject {
=======
fn read_bam(py: Python, path_or_file_like: PyObject, region: Option<&str>, index: Option<PyObject>) -> PyObject {
>>>>>>> 5252a90 (Turn from_path back into method)
    if let Ok(string_ref) = path_or_file_like.downcast::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = BamReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        let ipc = reader.records_to_ipc(region).unwrap();
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    } else {
        // Otherwise, treat it as file-like
        let file_like = match PyFileLikeObject::new(path_or_file_like, true, false, true) {
            Ok(file_like) => file_like,
            Err(_) => panic!("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object."),
        };
        let index_file_like = match PyFileLikeObject::new(index.unwrap(), true, false, true) {
            Ok(file_like) => file_like,
            Err(_) => panic!("Unknown argument for `index`. Not a file path string or url, and not a file-like object."),
        };
        let index = bam::index_from_reader(index_file_like).unwrap();
        let mut reader = BamReader::new(file_like, index).unwrap();
        let ipc = reader.records_to_ipc(region).unwrap();
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    }
}

#[pyfunction]
fn read_bam_vpos(
    py: Python,
    path_or_file_like: PyObject,
    pos_lo: (u64, u16),
    pos_hi: (u64, u16),
    index: Option<PyObject>,
) -> PyObject {
    if let Ok(string_ref) = path_or_file_like.downcast::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = BamReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        let ipc = reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap();
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    } else {
        // Otherwise, treat it as file-like
        let file_like = match PyFileLikeObject::new(path_or_file_like, true, false, true) {
            Ok(file_like) => file_like,
            Err(_) => panic!("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object."),
        };
        let index_file_like = match PyFileLikeObject::new(index.unwrap(), true, false, true) {
            Ok(file_like) => file_like,
            Err(_) => panic!("Unknown argument for `index`. Not a file path string or url, and not a file-like object."),
        };
        let index = bam::index_from_reader(index_file_like).unwrap();
        let mut reader = BamReader::new(file_like, index).unwrap();
        let ipc = reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap();
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    }
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

#[pyfunction]
fn read_bigwig(py: Python, path_or_file_like: PyObject, region: Option<&str>) -> PyObject {
    if let Ok(string_ref) = path_or_file_like.downcast::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = BigWigReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        let ipc = reader.records_to_ipc(region).unwrap();
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    } else {
        // Otherwise, treat it as file-like
        let file_like = match PyFileLikeObject::new(path_or_file_like, true, false, true) {
            Ok(file_like) => file_like,
            Err(_) => panic!("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object."),
        };
        let mut reader = BigWigReader::new(file_like).unwrap();
        let ipc = reader.records_to_ipc(region).unwrap();
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    }
}

#[pyfunction]
fn read_bigbed(py: Python, path_or_file_like: PyObject, region: Option<&str>) -> PyObject {
    if let Ok(string_ref) = path_or_file_like.downcast::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = BigBedReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        let ipc = reader.records_to_ipc(region).unwrap();
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    } else {
        // Otherwise, treat it as file-like
        let file_like = match PyFileLikeObject::new(path_or_file_like, true, false, true) {
            Ok(file_like) => file_like,
            Err(_) => panic!("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object."),
        };
        let mut reader = BigBedReader::new(file_like).unwrap();
        let ipc = reader.records_to_ipc(region).unwrap();
        Python::with_gil(|py| PyBytes::new(py, &ipc).into())
    }
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
    m.add_function(wrap_pyfunction!(read_bigwig, m)?)?;
    m.add_function(wrap_pyfunction!(read_bigbed, m)?)?;
    Ok(())
}

use std::collections::HashSet;
use std::io::BufReader;

use pyo3::prelude::*;
use pyo3::types::PyString;

use oxbow::bam;
use oxbow::bam::BamReader;
use oxbow::bcf;
use oxbow::bigbed::BigBedReader;
use oxbow::bigwig::BigWigReader;
use oxbow::fasta::FastaReader;
use oxbow::fastq::FastqReader;
use oxbow::vcf;
// use oxbow::cram::CramReader;
use oxbow::bcf::BcfReader;
use oxbow::gff::GffReader;
use oxbow::gtf::GtfReader;
use oxbow::vcf::VcfReader;

use oxbow::vpos;

mod file_like;

use file_like::PyFileLikeObject;

/// Wraps a `PyFileLikeObject` in a `BufReader` with a 1MB buffer.
/// Ensures consistent buffering for Python file-like objects.
fn buffered_file_like(path_or_file_like: PyObject) -> PyResult<BufReader<PyFileLikeObject>> {
    PyFileLikeObject::new(path_or_file_like, true, false, true)
        // Alternative to `std::io::BufReader::new` with a larger buffer size (1MB instead of 8KB).
        .map(|file_like| BufReader::with_capacity(const { 1024 * 1024 }, file_like))
}

#[pyfunction]
fn partition_from_index_file(path: &str, chunksize: u64) -> Vec<(u64, u16)> {
    vpos::partition_from_index_file(path, chunksize)
}

#[pyfunction]
#[pyo3(signature = (path, region=None))]
fn read_fasta(path: &str, region: Option<&str>) -> Vec<u8> {
    let mut reader = FastaReader::new(path).unwrap();
    reader.records_to_ipc(region).unwrap()
}

#[pyfunction]
fn read_fastq(py: Python, path_or_file_like: PyObject) -> Vec<u8> {
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it like a path
        let mut reader = FastqReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        reader.records_to_ipc().unwrap()
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let mut reader = FastqReader::new(file_like).unwrap();
        reader.records_to_ipc().unwrap()
    }
}

#[pyfunction]
#[pyo3(signature = (path_or_file_like, region=None, index=None))]
fn read_bam(
    py: Python,
    path_or_file_like: PyObject,
    region: Option<&str>,
    index: Option<PyObject>,
) -> Vec<u8> {
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = BamReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        reader.records_to_ipc(region).unwrap()
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let index_file_like = buffered_file_like(index.unwrap())
            .expect("Unknown argument for `index`. Not a file path string or url, and not a file-like object.");
        let index = bam::index_from_reader(index_file_like).unwrap();
        let mut reader = BamReader::new(file_like, index).unwrap();
        reader.records_to_ipc(region).unwrap()
    }
}

#[pyfunction]
#[pyo3(signature = (path_or_file_like, pos_lo, pos_hi, index=None))]
fn read_bam_vpos(
    py: Python,
    path_or_file_like: PyObject,
    pos_lo: (u64, u16),
    pos_hi: (u64, u16),
    index: Option<PyObject>,
) -> Vec<u8> {
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = BamReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap()
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let index_file_like = buffered_file_like(index.unwrap())
            .expect("Unknown argument for `index`. Not a file path string or url, and not a file-like object.");
        let index = bam::index_from_reader(index_file_like).unwrap();
        let mut reader = BamReader::new(file_like, index).unwrap();
        reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap()
    }
}

#[pyfunction]
#[pyo3(signature = (path_or_file_like, region=None, index=None))]
fn read_vcf(
    py: Python,
    path_or_file_like: PyObject,
    region: Option<&str>,
    index: Option<PyObject>,
) -> Vec<u8> {
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = VcfReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        reader.records_to_ipc(region).unwrap()
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let index_file_like = buffered_file_like(index.unwrap())
            .expect("Unknown argument for `index`. Not a file path string or url, and not a file-like object.");
        let index = vcf::index_from_reader(index_file_like).unwrap();
        let mut reader = VcfReader::new(file_like, index).unwrap();
        reader.records_to_ipc(region).unwrap()
    }
}

#[pyfunction]
#[pyo3(signature = (path_or_file_like, pos_lo, pos_hi, index=None))]
fn read_vcf_vpos(
    py: Python,
    path_or_file_like: PyObject,
    pos_lo: (u64, u16),
    pos_hi: (u64, u16),
    index: Option<PyObject>,
) -> Vec<u8> {
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = VcfReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap()
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let index_file_like = buffered_file_like(index.unwrap())
            .expect("Unknown argument for `index`. Not a file path string or url, and not a file-like object.");
        let index = vcf::index_from_reader(index_file_like).unwrap();
        let mut reader = VcfReader::new(file_like, index).unwrap();
        reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap()
    }
}

#[pyfunction]
#[pyo3(signature = (path_or_file_like, region=None, index=None))]
fn read_bcf(
    py: Python,
    path_or_file_like: PyObject,
    region: Option<&str>,
    index: Option<PyObject>,
) -> Vec<u8> {
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = BcfReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        reader.records_to_ipc(region).unwrap()
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let index_file_like = buffered_file_like(index.unwrap())
            .expect("Unknown argument for `index`. Not a file path string or url, and not a file-like object.");
        let index = bcf::index_from_reader(index_file_like).unwrap();
        let mut reader = BcfReader::new(file_like, index).unwrap();
        reader.records_to_ipc(region).unwrap()
    }
}

#[pyfunction]
#[pyo3(signature = (path_or_file_like, pos_lo, pos_hi, index=None))]
fn read_bcf_vpos(
    py: Python,
    path_or_file_like: PyObject,
    pos_lo: (u64, u16),
    pos_hi: (u64, u16),
    index: Option<PyObject>,
) -> Vec<u8> {
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = VcfReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap()
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let index_file_like = buffered_file_like(index.unwrap())
            .expect("Unknown argument for `index`. Not a file path string or url, and not a file-like object.");
        let index = bcf::index_from_reader(index_file_like).unwrap();
        let mut reader = VcfReader::new(file_like, index).unwrap();
        reader.records_to_ipc_from_vpos(pos_lo, pos_hi).unwrap()
    }
}

#[pyfunction]
#[pyo3(signature = (path_or_file_like, region=None, zoom_level=None, zoom_summary_columns=None))]
fn read_bigwig(
    py: Python,
    path_or_file_like: PyObject,
    region: Option<&str>,
    zoom_level: Option<u32>,
    zoom_summary_columns: Option<HashSet<String>>,
) -> Vec<u8> {
    let zoom_summary_columns_ref = zoom_summary_columns
        .as_ref()
        .map(|h| h.iter().map(String::as_str).collect::<HashSet<&str>>());
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = BigWigReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        match zoom_level {
            Some(zoom_level) => reader
                .zoom_records_to_ipc(region, zoom_level, zoom_summary_columns_ref)
                .unwrap(),
            None => reader.records_to_ipc(region).unwrap(),
        }
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let mut reader = BigWigReader::new(file_like).unwrap();
        match zoom_level {
            Some(zoom_level) => reader
                .zoom_records_to_ipc(region, zoom_level, zoom_summary_columns_ref)
                .unwrap(),
            None => reader.records_to_ipc(region).unwrap(),
        }
    }
}

#[pyfunction]
#[pyo3(signature = (path_or_file_like, region=None, fields=None))]
fn read_bigbed(
    py: Python,
    path_or_file_like: PyObject,
    region: Option<&str>,
    fields: Option<HashSet<String>>,
) -> Vec<u8> {
    let fields_ref = fields
        .as_ref()
        .map(|h| h.iter().map(String::as_str).collect::<HashSet<&str>>());
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = BigBedReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        reader.records_to_ipc(region, fields_ref).unwrap()
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let mut reader = BigBedReader::new(file_like).unwrap();
        reader.records_to_ipc(region, fields_ref).unwrap()
    }
}

#[pyfunction]
fn read_gff(py: Python, path_or_file_like: PyObject) -> Vec<u8> {
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = GffReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        reader.records_to_ipc().unwrap()
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let mut reader = GffReader::new(file_like).unwrap();
        reader.records_to_ipc().unwrap()
    }
}

#[pyfunction]
fn read_gtf(py: Python, path_or_file_like: PyObject) -> Vec<u8> {
    if let Ok(string_ref) = path_or_file_like.downcast_bound::<PyString>(py) {
        // If it's a string, treat it as a path
        let mut reader = GtfReader::new_from_path(string_ref.to_str().unwrap()).unwrap();
        reader.records_to_ipc().unwrap()
    } else {
        // Otherwise, treat it as file-like
        let file_like = buffered_file_like(path_or_file_like)
            .expect("Unknown argument for `path_url_or_file_like`. Not a file path string or url, and not a file-like object.");
        let mut reader = GtfReader::new(file_like).unwrap();
        reader.records_to_ipc().unwrap()
    }
}

#[pymodule]
#[pyo3(name = "oxbow")]
fn py_oxbow(m: &Bound<'_, PyModule>) -> PyResult<()> {
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
    m.add_function(wrap_pyfunction!(read_gff, m)?)?;
    m.add_function(wrap_pyfunction!(read_gtf, m)?)?;
    Ok(())
}

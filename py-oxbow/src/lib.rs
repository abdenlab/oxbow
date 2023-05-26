use pyo3::prelude::*;
use pyo3::types::PyBytes;

use oxbow::bam::BamReader;
use oxbow::vcf::VcfReader;

#[pyfunction]
fn read_bam(path: &str, region: Option<&str>) -> PyObject {
    let mut reader = BamReader::new(path).unwrap();
    let ipc = reader.records_to_ipc(region).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pyfunction]
fn read_vcf(path: &str, region: Option<&str>) -> PyObject {
    let mut reader = VcfReader::new(path).unwrap();
    let ipc = reader.records_to_ipc(region).unwrap();
    Python::with_gil(|py| PyBytes::new(py, &ipc).into())
}

#[pymodule]
#[pyo3(name = "oxbow")]
fn py_oxbow(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_bam, m)?)?;
    m.add_function(wrap_pyfunction!(read_vcf, m)?)?;
    Ok(())
}

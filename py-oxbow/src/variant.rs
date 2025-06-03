use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use pyo3_arrow::PyRecordBatchReader;
use pyo3_arrow::PySchema;

use noodles::core::Region;

use crate::util::{pyobject_to_bufreader, resolve_index, Reader};
use oxbow::util::batches_to_ipc;
use oxbow::util::index::IndexType;
use oxbow::variant::{BcfScanner, GenotypeBy, VcfScanner};

/// A VCF file scanner.
///
/// Parameters
/// ----------
/// src : PyObject
///    The path to the VCF file or a file-like object.
/// compressed : bool, optional [default: False]
///    Whether the source is BGZF-compressed.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyVcfScanner {
    src: PyObject,
    reader: Reader,
    scanner: VcfScanner,
    compressed: bool,
}

#[pymethods]
impl PyVcfScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false))]
    fn new(py: Python, src: PyObject, compressed: bool) -> PyResult<Self> {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let mut fmt_reader = noodles::vcf::io::Reader::new(reader);
        let header = fmt_reader.read_header()?;
        let reader = fmt_reader.into_inner();
        let scanner = VcfScanner::new(header);
        Ok(Self {
            src,
            reader,
            scanner,
            compressed,
        })
    }

    fn __getstate__(&self, py: Python<'_>) -> PyResult<PyObject> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let args = (self.src.clone_ref(py), self.compressed.into_py_any(py)?);
        let kwargs = PyDict::new(py);
        Ok((args.into_py_any(py)?, kwargs.into_py_any(py)?))
    }

    /// Return the names of the references sequences.
    fn chrom_names(&self) -> Vec<String> {
        self.scanner.chrom_names()
    }

    /// Return the names of the reference sequences and their lengths in bp.
    fn chrom_sizes(&self) -> Vec<(String, u32)> {
        self.scanner.chrom_sizes()
    }

    /// Return the names of the fixed fields.
    fn field_names(&self) -> Vec<String> {
        self.scanner.field_names()
    }

    /// Return the definitions of the INFO fields.
    fn info_field_defs(&self) -> Vec<(String, String, String)> {
        self.scanner.info_field_defs()
    }

    /// Return the names of the INFO fields.
    fn info_field_names(&self) -> Vec<String> {
        self.scanner.info_field_names()
    }

    /// Return the definitions of the FORMAT fields.
    fn genotype_field_defs(&self) -> Vec<(String, String, String)> {
        self.scanner.genotype_field_defs()
    }

    /// Return the names of the FORMAT fields.
    fn genotype_field_names(&self) -> Vec<String> {
        self.scanner.genotype_field_names()
    }

    /// Return the names of the samples.
    fn sample_names(&self) -> Vec<String> {
        self.scanner.sample_names()
    }

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///    Names of the fixed fields to project.
    /// info_fields : list[str], optional
    ///    Names of the INFO fields to project.
    /// genotype_fields : list[str], optional
    ///    Names of the sample-specific genotype fields to project.
    /// samples : list[str], optional
    ///    Names of the samples to include in the genotype fields.
    /// genotype_by : Literal["sample", "field"], optional [default: "sample"]
    ///    How to project the genotype fields. If "sample", the columns
    ///    correspond to the samples. If "field", the columns correspond to
    ///    the genotype fields.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    #[pyo3(signature = (fields=None, info_fields=None, genotype_fields=None, samples=None, genotype_by=None))]
    fn schema(
        &self,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<String>,
    ) -> PyResult<PySchema> {
        let genotype_by = resolve_genotype_by(genotype_by)?;
        let schema =
            self.scanner
                .schema(fields, info_fields, genotype_fields, samples, genotype_by)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Scan batches of records from the file.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// info_fields : list[str], optional
    ///     Names of the INFO fields to project.
    /// genotype_fields : list[str], optional
    ///     Names of the sample-specific genotype fields to project.
    /// samples : list[str], optional
    ///     Names of the samples to include in the genotype fields.
    /// genotype_by : Literal["sample", "field"], optional [default: "sample"]
    ///     How to project the genotype fields. If "sample", the columns
    ///     correspond to the samples. If "field", the columns correspond to
    ///     the genotype fields.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan. If None, records are scanned
    ///     until EOF.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (fields=None, info_fields=None, genotype_fields=None, samples=None, genotype_by=None, batch_size=1024, limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn scan(
        &mut self,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<String>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let genotype_by = resolve_genotype_by(genotype_by)?;
        let reader = self.reader.clone();
        let fmt_reader = noodles::vcf::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan(
                fmt_reader,
                fields,
                info_fields,
                genotype_fields,
                samples,
                genotype_by,
                batch_size,
                limit,
            )
            .map_err(PyErr::new::<PyValueError, _>)?;
        let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
        Ok(py_batch_reader)
    }

    /// Scan batches of records from a genomic range query on a BGZF-encoded file.
    ///
    /// This operation requires an index file.
    ///
    /// Parameters
    /// ----------
    /// region : str
    ///     Genomic region in the format "chr:start-end".
    /// index : path or file-like, optional
    ///     The index file to use for querying the region. If None and the
    ///     source was provided as a path, we will attempt to load the index
    ///     from the same path with an additional extension.
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// info_fields : list[str], optional
    ///     Names of the INFO fields to project.
    /// genotype_fields : list[str], optional
    ///     Names of the sample-specific genotype fields to project.
    /// samples : list[str], optional
    ///     Names of the samples to include in the genotype fields.
    /// genotype_by : Literal["sample", "field"], optional [default: "sample"]
    ///     How to project the genotype fields. If "sample", the columns
    ///     correspond to the samples. If "field", the columns correspond to
    ///     the genotype fields.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (region, index=None, fields=None, info_fields=None, genotype_fields=None, samples=None, genotype_by=None, batch_size=1024, limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn scan_query(
        &mut self,
        py: Python,
        region: String,
        index: Option<PyObject>,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<String>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;
        let genotype_by = resolve_genotype_by(genotype_by)?;
        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, self.src.clone_ref(py), index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader,
                                region,
                                index,
                                fields,
                                info_fields,
                                genotype_fields,
                                samples,
                                genotype_by,
                                batch_size,
                                limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader,
                                region,
                                index,
                                fields,
                                info_fields,
                                genotype_fields,
                                samples,
                                genotype_by,
                                batch_size,
                                limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, self.src.clone_ref(py), index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader,
                                region,
                                index,
                                fields,
                                info_fields,
                                genotype_fields,
                                samples,
                                genotype_by,
                                batch_size,
                                limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader,
                                region,
                                index,
                                fields,
                                info_fields,
                                genotype_fields,
                                samples,
                                genotype_by,
                                batch_size,
                                limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            _ => Err(PyErr::new::<PyValueError, _>(
                "Scanning query ranges is only supported for indexed bgzf-compressed sources.",
            )),
        }
    }
}

/// A BCF file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the BCF file or a file-like object.
/// compressed : bool, optional [default: True]
///     Whether the source is BGZF-compressed. If None, it is assumed to be
///     compressed.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyBcfScanner {
    src: PyObject,
    reader: Reader,
    scanner: BcfScanner,
    compressed: bool,
}

#[pymethods]
impl PyBcfScanner {
    #[new]
    #[pyo3(signature = (src, compressed=true))]
    fn new(py: Python, src: PyObject, compressed: bool) -> PyResult<Self> {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let mut fmt_reader = noodles::bcf::io::Reader::from(reader);
        let header = fmt_reader.read_header()?;
        let reader = fmt_reader.into_inner();
        let scanner = BcfScanner::new(header);
        Ok(Self {
            src,
            reader,
            scanner,
            compressed,
        })
    }

    fn __getstate__(&self, py: Python<'_>) -> PyResult<PyObject> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let args = (self.src.clone_ref(py), self.compressed.into_py_any(py)?);
        let kwargs = PyDict::new(py);
        Ok((args.into_py_any(py)?, kwargs.into_py_any(py)?))
    }

    /// Return the names of the reference sequences.
    fn chrom_names(&self) -> Vec<String> {
        self.scanner.chrom_names()
    }

    /// Return the names of the reference sequences and their lengths in bp.
    fn chrom_sizes(&self) -> Vec<(String, u32)> {
        self.scanner.chrom_sizes()
    }

    /// Return the names of the fixed fields.
    fn field_names(&self) -> Vec<String> {
        self.scanner.field_names()
    }

    /// Return the definitions of the INFO fields.
    fn info_field_defs(&self) -> Vec<(String, String, String)> {
        self.scanner.info_field_defs()
    }

    /// Return the names of the INFO fields.
    fn info_field_names(&self) -> Vec<String> {
        self.scanner.info_field_names()
    }

    /// Return the definitions of the FORMAT fields.
    fn genotype_field_defs(&self) -> Vec<(String, String, String)> {
        self.scanner.genotype_field_defs()
    }

    /// Return the definitions of the FORMAT fields.
    fn genotype_field_names(&self) -> Vec<String> {
        self.scanner.genotype_field_names()
    }

    /// Return the names of the samples.
    fn sample_names(&self) -> Vec<String> {
        self.scanner.sample_names()
    }

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///    Names of the fixed fields to project.
    /// info_fields : list[str], optional
    ///    Names of the INFO fields to project.
    /// genotype_fields : list[str], optional
    ///    Names of the sample-specific genotype fields to project.
    /// samples : list[str], optional
    ///    Names of the samples to include in the genotype fields.
    /// genotype_by : Literal["sample", "field"], optional [default: "sample"]
    ///    How to project the genotype fields. If "sample", the columns
    ///    correspond to the samples. If "field", the columns correspond to
    ///    the genotype fields.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    #[pyo3(signature = (fields=None, info_fields=None, genotype_fields=None, samples=None, genotype_by=None))]
    fn schema(
        &self,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<String>,
    ) -> PyResult<PySchema> {
        let genotype_by = resolve_genotype_by(genotype_by)?;
        let schema =
            self.scanner
                .schema(fields, info_fields, genotype_fields, samples, genotype_by)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Scan batches of records from the file.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// info_fields : list[str], optional
    ///     Names of the INFO fields to project.
    /// genotype_fields : list[str], optional
    ///     Names of the sample-specific genotype fields to project.
    /// samples : list[str], optional
    ///     Names of the samples to include in the genotype fields.
    /// genotype_by : Literal["sample", "field"], optional [default: "sample"]
    ///     How to project the genotype fields. If "sample", the columns
    ///     correspond to the samples. If "field", the columns correspond to
    ///     the genotype fields.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan. If None, records are scanned
    ///     until EOF.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (fields=None, info_fields=None, genotype_fields=None, samples=None, genotype_by=None, batch_size=1024, limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn scan(
        &mut self,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<String>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let genotype_by = resolve_genotype_by(genotype_by)?;
        let reader = self.reader.clone();
        let fmt_reader = noodles::bcf::io::Reader::from(reader);
        let batch_reader = self
            .scanner
            .scan(
                fmt_reader,
                fields,
                info_fields,
                genotype_fields,
                samples,
                genotype_by,
                batch_size,
                limit,
            )
            .map_err(PyErr::new::<PyValueError, _>)?;
        let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
        Ok(py_batch_reader)
    }

    /// Scan batches of records from a genomic range query on a BGZF-encoded file.
    ///
    /// This operation requires an index file.
    ///
    /// Parameters
    /// ----------
    /// region : str
    ///     Genomic region in the format "chr:start-end".
    /// index : path or file-like, optional
    ///     The index file to use for querying the region. If None and the
    ///     source was provided as a path, we will attempt to load the index
    ///     from the same path with an additional extension.
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// info_fields : list[str], optional
    ///     Names of the INFO fields to project.
    /// genotype_fields : list[str], optional
    ///     Names of the sample-specific genotype fields to project.
    /// samples : list[str], optional
    ///     Names of the samples to include in the genotype fields.
    /// genotype_by : Literal["sample", "field"], optional [default: "sample"]
    ///     How to project the genotype fields. If "sample", the columns
    ///     correspond to the samples. If "field", the columns correspond to
    ///     the genotype fields.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (region, index=None, fields=None, info_fields=None, genotype_fields=None, samples=None, genotype_by=None, batch_size=1024, limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn scan_query(
        &mut self,
        py: Python,
        region: String,
        index: Option<PyObject>,
        fields: Option<Vec<String>>,
        info_fields: Option<Vec<String>>,
        genotype_fields: Option<Vec<String>>,
        samples: Option<Vec<String>>,
        genotype_by: Option<String>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;
        let genotype_by = resolve_genotype_by(genotype_by)?;

        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, self.src.clone_ref(py), index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader,
                                region,
                                index,
                                fields,
                                info_fields,
                                genotype_fields,
                                samples,
                                genotype_by,
                                batch_size,
                                limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader,
                                region,
                                index,
                                fields,
                                info_fields,
                                genotype_fields,
                                samples,
                                genotype_by,
                                batch_size,
                                limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, self.src.clone_ref(py), index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader,
                                region,
                                index,
                                fields,
                                info_fields,
                                genotype_fields,
                                samples,
                                genotype_by,
                                batch_size,
                                limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader,
                                region,
                                index,
                                fields,
                                info_fields,
                                genotype_fields,
                                samples,
                                genotype_by,
                                batch_size,
                                limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            _ => Err(PyErr::new::<PyValueError, _>(
                "Scanning query ranges is only supported for indexed bgzf-compressed sources.",
            )),
        }
    }
}

fn resolve_genotype_by(genotype_by: Option<String>) -> PyResult<Option<GenotypeBy>> {
    match genotype_by {
        Some(s) if s == *"sample" => Ok(Some(GenotypeBy::Sample)),
        Some(s) if s == *"field" => Ok(Some(GenotypeBy::Field)),
        None => Ok(None),
        _ => Err(PyErr::new::<PyValueError, _>(
            "genotype_by must be either 'sample' or 'field'.",
        )),
    }
}

/// Return Arrow IPC format from a VCF file.
///
/// Parameters
/// ----------
/// src: str or file-like
///     The path to the VCF file or a file-like object.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// info_fields : list[str], optional
///     Names of the INFO fields to project.
/// genotype_fields : list[str], optional
///     Names of the sample-specific genotype fields to project.
/// samples : list[str], optional
///     Names of the samples to include in the genotype fields.
/// genotype_by : Literal["sample", "field"], optional [default: "sample"]
///     How to project the genotype fields. If "sample", the columns
///     correspond to the samples. If "field", the columns correspond to
///     the genotype fields.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, region=None, index=None, fields=None, info_fields=None, genotype_fields=None, samples=None, genotype_by=None, compressed=false))]
#[allow(clippy::too_many_arguments)]
pub fn read_vcf(
    py: Python,
    src: PyObject,
    region: Option<String>,
    index: Option<PyObject>,
    fields: Option<Vec<String>>,
    info_fields: Option<Vec<String>>,
    genotype_fields: Option<Vec<String>>,
    samples: Option<Vec<String>>,
    genotype_by: Option<String>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let mut fmt_reader = noodles::vcf::io::Reader::new(reader);
    let header = fmt_reader.read_header()?;
    let scanner = VcfScanner::new(header);
    let reader = fmt_reader.into_inner();
    let genotype_by = resolve_genotype_by(genotype_by)?;

    let ipc = if let Some(region) = region {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, src.clone_ref(py), index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    info_fields,
                    genotype_fields,
                    samples,
                    genotype_by,
                    None,
                    None,
                )?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, src.clone_ref(py), index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    info_fields,
                    genotype_fields,
                    samples,
                    genotype_by,
                    None,
                    None,
                )?;
                batches_to_ipc(batches)
            }
            _ => {
                return Err(PyErr::new::<PyValueError, _>(
                    "Scanning query ranges is only supported for indexed BGZF-encoded sources.",
                ));
            }
        }
    } else {
        let fmt_reader = noodles::vcf::io::Reader::new(reader);
        let batches = scanner.scan(
            fmt_reader,
            fields,
            info_fields,
            genotype_fields,
            samples,
            genotype_by,
            None,
            None,
        )?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

/// Return Arrow IPC format from a BCF file.
///
/// Parameters
/// ----------
/// src: str or file-like
///     The path to the BCF file or a file-like object.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// info_fields : list[str], optional
///     Names of the INFO fields to project.
/// genotype_fields : list[str], optional
///     Names of the sample-specific genotype fields to project.
/// samples : list[str], optional
///     Names of the samples to include in the genotype fields.
/// genotype_by : Literal["sample", "field"], optional [default: "sample"]
///     How to project the genotype fields. If "sample", the columns
///     correspond to the samples. If "field", the columns correspond to
///     the genotype fields.
/// compressed : bool, optional [default: True]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, region=None, index=None, fields=None, info_fields=None, genotype_fields=None, samples=None, genotype_by=None, compressed=true))]
#[allow(clippy::too_many_arguments)]
pub fn read_bcf(
    py: Python,
    src: PyObject,
    region: Option<String>,
    index: Option<PyObject>,
    fields: Option<Vec<String>>,
    info_fields: Option<Vec<String>>,
    genotype_fields: Option<Vec<String>>,
    samples: Option<Vec<String>>,
    genotype_by: Option<String>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let mut fmt_reader = noodles::bcf::io::Reader::from(reader);
    let header = fmt_reader.read_header()?;
    let scanner = BcfScanner::new(header);
    let reader = fmt_reader.into_inner();
    let genotype_by = resolve_genotype_by(genotype_by)?;

    let ipc = if let Some(region) = region {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, src.clone_ref(py), index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    info_fields,
                    genotype_fields,
                    samples,
                    genotype_by,
                    None,
                    None,
                )?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, src.clone_ref(py), index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    info_fields,
                    genotype_fields,
                    samples,
                    genotype_by,
                    None,
                    None,
                )?;
                batches_to_ipc(batches)
            }
            _ => {
                return Err(PyErr::new::<PyValueError, _>(
                    "Scanning query ranges is only supported for indexed BGZF-encoded sources.",
                ));
            }
        }
    } else {
        let fmt_reader = noodles::bcf::io::Reader::from(reader);
        let batches = scanner.scan(
            fmt_reader,
            fields,
            info_fields,
            genotype_fields,
            samples,
            genotype_by,
            None,
            None,
        )?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

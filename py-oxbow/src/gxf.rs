use std::io::Seek;
use std::sync::Arc;

use noodles::bgzf::io::Seek as _;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use pyo3_arrow::PyRecordBatchReader;
use pyo3_arrow::PySchema;

use noodles::core::Region;

use crate::util::{pyobject_to_bufreader, resolve_index, Reader};
use oxbow::gxf::{GffScanner, GtfScanner};
use oxbow::util::batches_to_ipc;
use oxbow::util::index::IndexType;

/// A GTF file scanner.
///
/// Parameters
/// ----------
/// obj : str or file-like
///     The path to the GTF file or a file-like object.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed. If None, it is assumed to be
///     uncompressed.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyGtfScanner {
    src: PyObject,
    reader: Reader,
    scanner: GtfScanner,
    compressed: bool,
}

#[pymethods]
impl PyGtfScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false))]
    fn new(py: Python, src: PyObject, compressed: Option<bool>) -> PyResult<Self> {
        let compressed = compressed.unwrap_or(false);
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let scanner = GtfScanner::new(None);
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

    // fn chrom_names(&self) -> Vec<String> {
    //     let scanner = GtfScanner::new(None);
    //     scanner.chrom_names().unwrap()
    // }

    /// Return the names of the fixed fields.
    fn field_names(&self) -> Vec<String> {
        self.scanner.field_names()
    }

    /// Discover attribute definitions by sniffing `scan_rows` records.
    ///
    /// The reader stream is reset to its original position after scanning.
    ///
    /// Parameters
    /// ----------
    /// scan_rows : int, optional [default: 1024]
    ///    The number of records to scan. If None, all records are scanned.
    ///
    /// Returns
    /// -------
    /// list[tuple[str, str]]
    ///     A list of attribute definitions, where each definition is a tuple
    ///     of the attribute name and its type (always ``String`` for GTF).
    #[pyo3(signature = (scan_rows=1024))]
    fn attribute_defs(&mut self, scan_rows: Option<usize>) -> PyResult<Vec<(String, String)>> {
        let mut reader = self.reader.clone();
        match &mut reader {
            Reader::BgzfFile(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
                let defs = self.scanner.attribute_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
                let defs = self.scanner.attribute_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            _ => {
                let pos = reader.stream_position()?;
                let mut fmt_reader = noodles::gtf::io::Reader::new(reader);
                let defs = self.scanner.attribute_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader
                    .into_inner()
                    .seek(std::io::SeekFrom::Start(pos))?;
                Ok(defs)
            }
        }
    }

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///    Names of the fixed fields to project.
    /// attribute_defs : list[tuple[str, str]], optional
    ///    Definitions of attribute fields to project.
    ///
    /// Returns
    /// -------
    /// pyo3_arrow.PySchema
    #[pyo3(signature = (fields=None, attribute_defs=None))]
    fn schema(
        &self,
        fields: Option<Vec<String>>,
        attribute_defs: Option<Vec<(String, String)>>,
    ) -> PyResult<PySchema> {
        let schema = self.scanner.schema(fields, attribute_defs)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Scan batches of records from the file.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// attribute_defs : list[tuple[str, str]], optional
    ///    Definitions of attribute fields to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan. If None, records are scanned
    ///     until EOF.
    ///
    /// Returns
    /// -------
    /// pyo3_arrow.PyRecordBatchReader
    ///     A PyCapsule stream iterator for the record batches.
    #[pyo3(signature = (fields=None, attribute_defs=None,batch_size=1024, limit=None))]
    fn scan(
        &mut self,
        fields: Option<Vec<String>>,
        attribute_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let fmt_reader = noodles::gtf::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, fields, attribute_defs, batch_size, limit)
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
    /// attribute_defs : list[tuple[str, str]], optional
    ///    Definitions of attribute fields to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// pyo3_arrow.PyRecordBatchReader
    ///     A PyCapsule stream iterator for the record batches.
    #[pyo3(signature = (region, index=None, fields=None, attribute_defs=None, batch_size=1024, limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn scan_query(
        &mut self,
        py: Python,
        region: &str,
        index: Option<PyObject>,
        fields: Option<Vec<String>>,
        attribute_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
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
                                attribute_defs,
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
                                attribute_defs,
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
                let fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
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
                                attribute_defs,
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
                                attribute_defs,
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

/// A GFF file scanner.
///
/// Parameters
/// ----------
/// obj : str or file-like
///     The path to the GTF file or a file-like object.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed. If None, it is assumed to be
///     uncompressed.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyGffScanner {
    src: PyObject,
    reader: Reader,
    scanner: GffScanner,
    compressed: bool,
}

#[pymethods]
impl PyGffScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false))]
    fn new(py: Python, src: PyObject, compressed: Option<bool>) -> PyResult<Self> {
        let compressed = compressed.unwrap_or(false);
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed).unwrap();
        let scanner = GffScanner::new(None);
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

    /// Return the names of the fixed fields.
    fn field_names(&self) -> Vec<String> {
        self.scanner.field_names()
    }

    /// Discover attribute definitions by sniffing `scan_rows` records.
    ///
    /// The reader stream is reset to its original position after scanning.
    ///
    /// Parameters
    /// ----------
    /// scan_rows : int, optional [default: 1024]
    ///    The number of records to scan. If None, all records are scanned.
    ///
    /// Returns
    /// -------
    /// list[tuple[str, str]]
    ///     A list of attribute definitions, where each definition is a tuple
    ///     of the attribute name and its type (``String`` or ``Array``).
    #[pyo3(signature = (scan_rows=1024))]
    fn attribute_defs(&mut self, scan_rows: Option<usize>) -> PyResult<Vec<(String, String)>> {
        let mut reader = self.reader.clone();
        match &mut reader {
            Reader::BgzfFile(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
                let defs = self.scanner.attribute_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
                let defs = self.scanner.attribute_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            _ => {
                let pos = reader.stream_position()?;
                let mut fmt_reader = noodles::gff::io::Reader::new(reader);
                let defs = self.scanner.attribute_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader
                    .into_inner()
                    .seek(std::io::SeekFrom::Start(pos))?;
                Ok(defs)
            }
        }
    }

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///    Names of the fixed fields to project.
    /// attribute_defs : list[tuple[str, str]], optional
    ///    Definitions of attribute fields to project.
    ///
    /// Returns
    /// -------
    /// pyo3_arrow.PySchema
    #[pyo3(signature = (fields=None, attribute_defs=None))]
    fn schema(
        &self,
        fields: Option<Vec<String>>,
        attribute_defs: Option<Vec<(String, String)>>,
    ) -> PyResult<PySchema> {
        let schema = self.scanner.schema(fields, attribute_defs)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Scan batches of records from the file.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// attribute_defs : list[tuple[str, str]], optional
    ///    Definitions of attribute fields to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan. If None, records are scanned
    ///     until EOF.
    ///
    /// Returns
    /// -------
    /// pyo3_arrow.PyRecordBatchReader
    ///     A PyCapsule stream iterator for the record batches.
    #[pyo3(signature = (fields=None, attribute_defs=None, batch_size=1024, limit=None))]
    fn scan(
        &mut self,
        fields: Option<Vec<String>>,
        attribute_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let fmt_reader = noodles::gff::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, fields, attribute_defs, batch_size, limit)
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
    /// attribute_defs : list[tuple[str, str]], optional
    ///    Definitions of attribute fields to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// pyo3_arrow.PyRecordBatchReader
    ///     A PyCapsule stream iterator for the record batches.
    #[pyo3(signature = (region, index=None, fields=None, attribute_defs=None, batch_size=1024, limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn scan_query(
        &mut self,
        py: Python,
        region: &str,
        index: Option<PyObject>,
        fields: Option<Vec<String>>,
        attribute_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
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
                                attribute_defs,
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
                                attribute_defs,
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
                let fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
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
                                attribute_defs,
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
                                attribute_defs,
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

/// Return Arrow IPC format from a GTF file.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the source file or a file-like object.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// tag_defs : list[tuple[str, str]], optional
///    Definitions of attribute fields to project.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, region=None, index=None, fields=None, attr_defs=None, compressed=false))]
pub fn read_gtf(
    py: Python,
    src: PyObject,
    region: Option<String>,
    index: Option<PyObject>,
    fields: Option<Vec<String>>,
    attr_defs: Option<Vec<(String, String)>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let scanner = GtfScanner::new(None);

    let ipc = if let Some(region) = region {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, src.clone_ref(py), index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    attr_defs,
                    None,
                    None,
                )?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, src.clone_ref(py), index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    attr_defs,
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
        let fmt_reader = noodles::gtf::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, fields, None, None, None)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

/// Return Arrow IPC format from a GFF file.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the source file or a file-like object.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// attr_defs : list[tuple[str, str]], optional
///    Definitions of attribute fields to project.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, region=None, index=None, fields=None, attr_defs=None, compressed=false))]
pub fn read_gff(
    py: Python,
    src: PyObject,
    region: Option<String>,
    index: Option<PyObject>,
    fields: Option<Vec<String>>,
    attr_defs: Option<Vec<(String, String)>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let scanner = GffScanner::new(None);

    let ipc = if let Some(region) = region {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, src.clone_ref(py), index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    attr_defs,
                    None,
                    None,
                )?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, src.clone_ref(py), index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    attr_defs,
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
        let fmt_reader = noodles::gff::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, fields, None, None, None)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

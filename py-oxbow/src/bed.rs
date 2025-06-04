use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use pyo3_arrow::PyRecordBatchReader;
use pyo3_arrow::PySchema;

use noodles::core::Region;

use crate::util::{pyobject_to_bufreader, resolve_index, Reader};
use oxbow::bed::{BedScanner, BedSchema};
use oxbow::util::batches_to_ipc;
use oxbow::util::index::IndexType;

/// A BED file scanner.
///
/// Parameters
/// ----------
/// stc : str or file-like
///     The path to the BED file or a file-like object.
/// bed_schema : str
///     The BED schema specifier, e.g., "bed6+3".
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed. If None, it is assumed to be
///     uncompressed.
///
/// Notes
/// -----
/// The BED schema specifier can be one of the following (case-insensitive):
///
/// - ``bed``: Equivalent to ``BED6``.
/// - ``bed{n}``: `n` standard fields and 0 custom fields.
/// - ``bed{n}+{m}``: `n` standard fields followed by `m` custom fields.
/// - ``bed{n}+``: `n` standard fields followed by an undefined number of custom fields.
///
/// While the 12 standard fields have defined types, custom fields are
/// intepreted as text. ``bed{n}+`` custom fields are collapsed into a single
/// field named `rest`.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyBedScanner {
    src: PyObject,
    reader: Reader,
    scanner: BedScanner,
    bed_schema: String,
    compressed: bool,
}

#[pymethods]
impl PyBedScanner {
    #[new]
    #[pyo3(signature = (src, bed_schema, compressed=false))]
    fn new(py: Python, src: PyObject, bed_schema: String, compressed: bool) -> PyResult<Self> {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let _bed_schema: BedSchema = bed_schema.parse().unwrap();
        let scanner = BedScanner::new(_bed_schema);
        Ok(Self {
            src,
            reader,
            scanner,
            bed_schema,
            compressed,
        })
    }

    fn __getstate__(&self, py: Python) -> PyResult<PyObject> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python) -> PyResult<(PyObject, PyObject)> {
        let args = (
            self.src.clone_ref(py),
            self.bed_schema.clone().into_py_any(py)?,
            self.compressed.into_py_any(py)?,
        );
        let kwargs = PyDict::new(py);
        Ok((args.into_py_any(py)?, kwargs.into_py_any(py)?))
    }

    // fn chrom_names(&self) -> Vec<String> {
    //     let scanner = BedScanner::new(None);
    //     scanner.chrom_names().unwrap()
    // }

    /// Return the names of the BED fields.
    fn field_names(&self) -> Vec<String> {
        self.scanner.field_names()
    }

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///    Names of the BED fields to project.
    /// tag_defs : list[tuple[str, str]], optional
    ///    Definitions of tag fields to project.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    #[pyo3(signature = (fields=None))]
    fn schema(&self, fields: Option<Vec<String>>) -> PyResult<PySchema> {
        let schema = self.scanner.schema(fields)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Scan batches of records from the file.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the BED fields to project.
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
    #[pyo3(signature = (fields=None, batch_size=1024, limit=None))]
    fn scan(
        &mut self,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let fmt_reader = noodles::bed::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, fields, batch_size, limit)
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
    ///     Names of the BED fields to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (region, index=None, fields=None, batch_size=1024, limit=None))]
    fn scan_query(
        &mut self,
        py: Python,
        region: &str,
        index: Option<PyObject>,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, self.src.clone_ref(py), index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(fmt_reader, region, index, fields, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(fmt_reader, region, index, fields, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, self.src.clone_ref(py), index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(fmt_reader, region, index, fields, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(fmt_reader, region, index, fields, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            _ => Err(PyErr::new::<PyValueError, _>(
                "Scanning query ranges is only supported for bgzf-compressed sources.",
            )),
        }
    }
}

/// Return Arrow IPC format from a BED file.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the source file or a file-like object.
/// bed_schema : str
///     The BED schema specifier, e.g., "bed6+3".
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
///
/// Notes
/// -----
/// The BED schema specifier can be one of the following (case-insensitive):
///
/// - ``bed``: Equivalent to ``BED6``.
/// - ``bed{n}``: `n` standard fields and 0 custom fields.
/// - ``bed{n}+{m}``: `n` standard fields followed by `m` custom fields.
/// - ``bed{n}+``: `n` standard fields followed by an undefined number of custom fields.
///
/// While the 12 standard fields have defined types, custom fields are
/// intepreted as text. ``bed{n}+`` custom fields are collapsed into a single
/// field named `rest`.
#[pyfunction]
#[pyo3(signature = (src, bed_schema, region=None, index=None, fields=None, compressed=false))]
pub fn read_bed(
    py: Python,
    src: PyObject,
    bed_schema: String,
    region: Option<String>,
    index: Option<PyObject>,
    fields: Option<Vec<String>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let bed_schema: BedSchema = bed_schema.parse().unwrap();
    let scanner = BedScanner::new(bed_schema);

    let ipc = if let Some(region) = region {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, src.clone_ref(py), index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    None,
                    None,
                )?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, src.clone_ref(py), index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
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
        let fmt_reader = noodles::bed::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, fields, None, None)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

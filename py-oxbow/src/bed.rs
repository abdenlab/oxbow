use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use pyo3_arrow::PyRecordBatchReader;
use pyo3_arrow::PySchema;

use crate::error::{err_on_unwind, to_py};
use crate::util::resolve_coord_system;
use crate::util::{
    pyobject_to_bufreader, resolve_fields, resolve_index, PyVirtualPosition, Reader,
};
use oxbow::bed::{BedScanner, BedSchema, FieldDef, FieldType};
use oxbow::util::batches_to_ipc;
use oxbow::util::index::IndexType;
use oxbow::CoordSystem;

/// Extract custom field definitions from a Python list or dict.
fn extract_custom_defs(obj: &Bound<'_, PyAny>) -> PyResult<Vec<FieldDef>> {
    if let Ok(dict) = obj.cast::<pyo3::types::PyDict>() {
        return dict
            .iter()
            .map(|(k, v)| {
                let name: String = k.extract()?;
                let ty_str: String = v.extract()?;
                let ty: FieldType = ty_str.parse().map_err(to_py)?;
                Ok(FieldDef::new(name, ty))
            })
            .collect::<PyResult<Vec<_>>>();
    }

    if let Ok(list) = obj.cast::<pyo3::types::PyList>() {
        return list
            .iter()
            .map(|item| {
                let (name, ty_str): (String, String) = item.extract()?;
                let ty: FieldType = ty_str.parse().map_err(to_py)?;
                Ok(FieldDef::new(name, ty))
            })
            .collect::<PyResult<Vec<_>>>();
    }

    Err(PyErr::new::<PyValueError, _>(
        "Custom field definitions must be a list[tuple[str, str]] or dict[str, str]",
    ))
}

/// Parse a BED schema from a Python object.
///
/// Accepts:
/// - `str`: a BED schema specifier (e.g., "bed6+3", "bedgraph")
/// - `tuple[str, list | dict]`: a base specifier + custom field definitions
///   (e.g., `("bed6", [("signalValue", "float"), ...])`)
pub fn resolve_bed_schema(py: Python, obj: &Py<PyAny>) -> PyResult<BedSchema> {
    let obj = obj.bind(py);

    if let Ok(s) = obj.extract::<String>() {
        return s.parse::<BedSchema>().map_err(to_py);
    }

    if let Ok(tuple) = obj.cast::<pyo3::types::PyTuple>() {
        if tuple.len() != 2 {
            return Err(PyErr::new::<PyValueError, _>(
                "Schema tuple must have exactly 2 elements: (base_specifier, custom_defs)",
            ));
        }
        let base: String = tuple.get_item(0)?.extract()?;
        let base_schema: BedSchema = base.parse().map_err(to_py)?;
        let n = base_schema.standard_field_count();
        let custom_defs = extract_custom_defs(&tuple.get_item(1)?)?;
        return BedSchema::new(n, Some(custom_defs)).map_err(to_py);
    }

    Err(PyErr::new::<PyValueError, _>(
        "bed_schema must be a str or a tuple of (str, list | dict)",
    ))
}

/// A BED file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the BED file or a file-like object.
/// bed_schema : str, list[tuple[str, str]], or dict[str, str]
///     The BED schema. Can be a specifier string (e.g., "bed6+3"), a list
///     of (name, type) tuples, or a dict mapping names to types.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
/// fields : list[str], optional
///     Names of the BED fields to include in the schema.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyBedScanner {
    src: Py<PyAny>,
    reader: Reader,
    scanner: BedScanner,
    compressed: bool,
}

#[pymethods]
impl PyBedScanner {
    #[new]
    #[pyo3(signature = (src, bed_schema, compressed=false, fields=None, coords=None))]
    fn new(
        py: Python,
        src: Py<PyAny>,
        bed_schema: Py<PyAny>,
        compressed: bool,
        fields: Option<Py<PyAny>>,
        coords: Option<String>,
    ) -> PyResult<Self> {
        let coord_system = resolve_coord_system(coords)?;
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let parsed_schema = resolve_bed_schema(py, &bed_schema)?;
        let scanner = BedScanner::new(
            parsed_schema,
            resolve_fields(fields, py)?,
            coord_system.unwrap_or(CoordSystem::ZeroHalfOpen),
        )
        .map_err(to_py)?;
        Ok(Self {
            src,
            reader,
            scanner,
            compressed,
        })
    }

    fn __getstate__(&self, py: Python) -> PyResult<Py<PyAny>> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
        let model = self.scanner.model();
        let args = (
            self.src.clone_ref(py),
            model.bed_schema().to_string().into_py_any(py)?,
        );
        let kwargs = PyDict::new(py);
        kwargs.set_item("compressed", self.compressed)?;
        kwargs.set_item("fields", model.field_names())?;
        kwargs.set_item("coords", model.coord_system().to_string())?;
        Ok((args.into_py_any(py)?, kwargs.into_py_any(py)?))
    }

    /// Return the names of the BED fields.
    fn field_names(&self) -> Vec<String> {
        self.scanner.field_names()
    }

    /// Return the Arrow schema.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    fn schema(&self) -> PySchema {
        PySchema::new(Arc::new(self.scanner.schema().clone()))
    }

    /// Scan batches of records from the file.
    ///
    /// Parameters
    /// ----------
    /// columns : list[str], optional
    ///     Names of the columns to project.
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
    #[pyo3(signature = (columns=None, batch_size=1024, limit=None))]
    fn scan(
        &mut self,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let fmt_reader = noodles::bed::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, columns, batch_size, limit)
            .map_err(to_py)?;
        Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
    }

    /// Scan batches of records from specified byte ranges in the file.
    ///
    /// The byte positions must align with record boundaries.
    ///
    /// Parameters
    /// ----------
    /// byte_ranges : list[tuple[int, int]]
    ///     List of (start, end) byte position tuples to read from.
    /// columns : list[str], optional
    ///     Names of the columns to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan. If None, all records
    ///     in the specified ranges are scanned.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (byte_ranges, columns=None, batch_size=1024, limit=None))]
    fn scan_byte_ranges(
        &mut self,
        byte_ranges: Vec<(u64, u64)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let fmt_reader = noodles::bed::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan_byte_ranges(fmt_reader, byte_ranges, columns, batch_size, limit)
            .map_err(to_py)?;
        Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
    }

    /// Scan batches of records from virtual position ranges in a BGZF file.
    ///
    /// The virtual positions must align with record boundaries. That means
    /// that the compressed offset must point to the beginning of a BGZF block
    /// and the uncompressed offset must point to the beginning or end of a
    /// record decoded within the block.
    ///
    /// Parameters
    /// ----------
    /// vpos_ranges : list[tuple[vpos, vpos]]
    ///     List of virtual position ranges as pairs. Each virtual position can
    ///     be given as either a packed virtual position (int), or an unpacked
    ///     tuple of ints ``(c, u)`` specifying the compressed and uncompressed
    ///     offsets, respectively.
    /// columns : list[str], optional
    ///     Names of the columns to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan. If None, all records
    ///     in the specified ranges are scanned.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (vpos_ranges, columns=None, batch_size=1024, limit=None))]
    fn scan_virtual_ranges(
        &mut self,
        vpos_ranges: Vec<(PyVirtualPosition, PyVirtualPosition)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let vpos_ranges = vpos_ranges
            .into_iter()
            .map(|(start, end)| Ok((start.to_virtual_position()?, end.to_virtual_position()?)))
            .collect::<PyResult<Vec<_>>>()?;
        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
                let batch_reader = self
                    .scanner
                    .scan_virtual_ranges(fmt_reader, vpos_ranges, columns, batch_size, limit)
                    .map_err(to_py)?;
                Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
                let batch_reader = self
                    .scanner
                    .scan_virtual_ranges(fmt_reader, vpos_ranges, columns, batch_size, limit)
                    .map_err(to_py)?;
                Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
            }
            _ => Err(PyErr::new::<PyValueError, _>(
                "Scanning virtual position ranges is only supported for bgzf-compressed sources.",
            )),
        }
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
    /// columns : list[str], optional
    ///     Names of the columns to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan. If None, all records
    ///     intersecting the query range are scanned.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (region, index=None, columns=None, batch_size=1024, limit=None))]
    fn scan_query(
        &mut self,
        py: Python,
        region: &str,
        index: Option<Py<PyAny>>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region =
            oxbow::Region::parse(region, self.scanner.model().coord_system()).map_err(to_py)?;

        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &self.src, index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(fmt_reader, region, index, columns, batch_size, limit)
                            .map_err(to_py)?;
                        PyRecordBatchReader::new(err_on_unwind(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(fmt_reader, region, index, columns, batch_size, limit)
                            .map_err(to_py)?;
                        PyRecordBatchReader::new(err_on_unwind(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &self.src, index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(fmt_reader, region, index, columns, batch_size, limit)
                            .map_err(to_py)?;
                        PyRecordBatchReader::new(err_on_unwind(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(fmt_reader, region, index, columns, batch_size, limit)
                            .map_err(to_py)?;
                        PyRecordBatchReader::new(err_on_unwind(batch_reader))
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
/// bed_schema : str, list[tuple[str, str]], or dict[str, str]
///     The BED schema.
/// fields : list[str], optional
///     Names of the fields to project.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, bed_schema, region=None, index=None, fields=None, compressed=false))]
pub fn read_bed(
    py: Python,
    src: Py<PyAny>,
    bed_schema: Py<PyAny>,
    region: Option<String>,
    index: Option<Py<PyAny>>,
    fields: Option<Py<PyAny>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let bed_schema = resolve_bed_schema(py, &bed_schema)?;
    let scanner = BedScanner::new(
        bed_schema,
        resolve_fields(fields, py)?,
        CoordSystem::ZeroHalfOpen,
    )
    .map_err(to_py)?;

    let ipc = if let Some(region) = region {
        let region =
            oxbow::Region::parse(&region, oxbow::CoordSystem::ZeroHalfOpen).map_err(to_py)?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &src, index)?;
                let batches = scanner
                    .scan_query(fmt_reader, region, index.into_boxed(), None, None, None)
                    .map_err(to_py)?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bed::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &src, index)?;
                let batches = scanner
                    .scan_query(fmt_reader, region, index.into_boxed(), None, None, None)
                    .map_err(to_py)?;
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
        let batches = scanner.scan(fmt_reader, None, None, None).map_err(to_py)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

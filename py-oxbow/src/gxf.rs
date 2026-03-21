use std::io::Seek;
use std::sync::Arc;

use noodles::bgzf::io::Seek as _;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use pyo3_arrow::PyRecordBatchReader;
use pyo3_arrow::PySchema;

use crate::error::{err_on_unwind, to_py};
use crate::util::{
    pyobject_to_bufreader, resolve_coord_system, resolve_fields, resolve_index, PyVirtualPosition,
    Reader,
};
use oxbow::gxf::{GffScanner, GtfScanner};
use oxbow::util::batches_to_ipc;
use oxbow::util::index::IndexType;
use oxbow::CoordSystem;

/// A GTF file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the GTF file or a file-like object.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// attribute_defs : list[tuple[str, str]], optional [default: None]
///     Definitions for the ``"attributes"`` struct column. ``None`` omits the
///     attributes column. Use the ``attribute_defs()`` method to discover definitions.
/// coords : Literal["01", "11"], optional [default: "11"]
///    Coordinate system for returning positions and interpreting query ranges.
///    "01" for 0-based half-open, "11" for 1-based closed.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyGtfScanner {
    src: Py<PyAny>,
    reader: Reader,
    scanner: GtfScanner,
    compressed: bool,
}

#[pymethods]
impl PyGtfScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false, fields=None, attribute_defs=None, coords=None))]
    fn new(
        py: Python,
        src: Py<PyAny>,
        compressed: Option<bool>,
        fields: Option<Py<PyAny>>,
        attribute_defs: Option<Vec<(String, String)>>,
        coords: Option<String>,
    ) -> PyResult<Self> {
        let fields = resolve_fields(fields, py)?;
        let coord_system = resolve_coord_system(coords)?;
        let compressed = compressed.unwrap_or(false);
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let scanner = GtfScanner::new(
            None,
            fields,
            attribute_defs,
            coord_system.unwrap_or(CoordSystem::OneClosed),
        )
        .map_err(to_py)?;
        Ok(Self {
            src,
            reader,
            scanner,
            compressed,
        })
    }

    fn __getstate__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python<'_>) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
        let args = (self.src.clone_ref(py), self.compressed.into_py_any(py)?);
        let kwargs = PyDict::new(py);
        let model = self.scanner.model();
        kwargs.set_item("fields", model.field_names())?;
        if let Some(defs) = model.attr_defs() {
            let attr_defs: Vec<_> = defs.iter().map(|d| d.to_tuple()).collect();
            kwargs.set_item("attribute_defs", attr_defs)?;
        }
        kwargs.set_item("coords", model.coord_system().to_string())?;
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
    ///     of the attribute name and its type (always ``String`` for GTF).
    #[pyo3(signature = (scan_rows=1024))]
    fn attribute_defs(&mut self, scan_rows: Option<usize>) -> PyResult<Vec<(String, String)>> {
        let mut reader = self.reader.clone();
        match &mut reader {
            Reader::BgzfFile(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
                let defs = GtfScanner::attribute_defs(&mut fmt_reader, scan_rows).map_err(to_py)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
                let defs = GtfScanner::attribute_defs(&mut fmt_reader, scan_rows).map_err(to_py)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            _ => {
                let pos = reader.stream_position()?;
                let mut fmt_reader = noodles::gtf::io::Reader::new(reader);
                let defs = GtfScanner::attribute_defs(&mut fmt_reader, scan_rows).map_err(to_py)?;
                fmt_reader
                    .into_inner()
                    .seek(std::io::SeekFrom::Start(pos))?;
                Ok(defs)
            }
        }
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
        let fmt_reader = noodles::gtf::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, columns, batch_size, limit)
            .map_err(to_py)?;
        Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
    }

    /// Scan batches of records from specified byte ranges in the file.
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
    ///     The maximum number of records to scan.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    #[pyo3(signature = (byte_ranges, columns=None, batch_size=1024, limit=None))]
    fn scan_byte_ranges(
        &mut self,
        byte_ranges: Vec<(u64, u64)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let fmt_reader = noodles::gtf::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan_byte_ranges(fmt_reader, byte_ranges, columns, batch_size, limit)
            .map_err(to_py)?;
        Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
    }

    /// Scan batches of records from virtual position ranges in a BGZF file.
    ///
    /// Parameters
    /// ----------
    /// vpos_ranges : list[tuple[vpos, vpos]]
    ///     List of virtual position ranges as pairs.
    /// columns : list[str], optional
    ///     Names of the columns to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
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
                let fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
                let batch_reader = self
                    .scanner
                    .scan_virtual_ranges(fmt_reader, vpos_ranges, columns, batch_size, limit)
                    .map_err(to_py)?;
                Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
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
    /// Parameters
    /// ----------
    /// region : str
    ///     Genomic range string in the format "chr:start-end",
    ///     "chr:[start,end]" or "chr:[start,end)".
    /// index : path or file-like, optional
    ///     The index file to use for querying the region.
    /// columns : list[str], optional
    ///     Names of the columns to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
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
                let fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
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
                let fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
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
                "Scanning query ranges is only supported for indexed bgzf-compressed sources.",
            )),
        }
    }
}

/// A GFF file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the GFF file or a file-like object.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// attribute_defs : list[tuple[str, str]], optional [default: None]
///     Definitions for the ``"attributes"`` struct column. ``None`` omits the
///     attributes column. Use the ``attribute_defs()`` method to discover definitions.
/// coords : Literal["01", "11"], optional [default: "11"]
///    Coordinate system for returning positions and interpreting query ranges.
///    "01" for 0-based half-open, "11" for 1-based closed.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyGffScanner {
    src: Py<PyAny>,
    reader: Reader,
    scanner: GffScanner,
    compressed: bool,
}

#[pymethods]
impl PyGffScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false, fields=None, attribute_defs=None, coords=None))]
    fn new(
        py: Python,
        src: Py<PyAny>,
        compressed: Option<bool>,
        fields: Option<Py<PyAny>>,
        attribute_defs: Option<Vec<(String, String)>>,
        coords: Option<String>,
    ) -> PyResult<Self> {
        let fields = resolve_fields(fields, py)?;
        let coord_system = resolve_coord_system(coords)?;
        let compressed = compressed.unwrap_or(false);
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let scanner = GffScanner::new(
            None,
            fields,
            attribute_defs,
            coord_system.unwrap_or(CoordSystem::OneClosed),
        )
        .map_err(to_py)?;
        Ok(Self {
            src,
            reader,
            scanner,
            compressed,
        })
    }

    fn __getstate__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python<'_>) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
        let args = (self.src.clone_ref(py), self.compressed.into_py_any(py)?);
        let kwargs = PyDict::new(py);
        let model = self.scanner.model();
        kwargs.set_item("fields", model.field_names())?;
        if let Some(defs) = model.attr_defs() {
            let attr_defs: Vec<_> = defs.iter().map(|d| d.to_tuple()).collect();
            kwargs.set_item("attribute_defs", attr_defs)?;
        }
        kwargs.set_item("coords", model.coord_system().to_string())?;
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
                let defs = GffScanner::attribute_defs(&mut fmt_reader, scan_rows).map_err(to_py)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
                let defs = GffScanner::attribute_defs(&mut fmt_reader, scan_rows).map_err(to_py)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            _ => {
                let pos = reader.stream_position()?;
                let mut fmt_reader = noodles::gff::io::Reader::new(reader);
                let defs = GffScanner::attribute_defs(&mut fmt_reader, scan_rows).map_err(to_py)?;
                fmt_reader
                    .into_inner()
                    .seek(std::io::SeekFrom::Start(pos))?;
                Ok(defs)
            }
        }
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
    #[pyo3(signature = (columns=None, batch_size=1024, limit=None))]
    fn scan(
        &mut self,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let fmt_reader = noodles::gff::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, columns, batch_size, limit)
            .map_err(to_py)?;
        Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
    }

    /// Scan batches of records from specified byte ranges in the file.
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
    ///     The maximum number of records to scan.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    #[pyo3(signature = (byte_ranges, columns=None, batch_size=1024, limit=None))]
    fn scan_byte_ranges(
        &mut self,
        byte_ranges: Vec<(u64, u64)>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let fmt_reader = noodles::gff::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan_byte_ranges(fmt_reader, byte_ranges, columns, batch_size, limit)
            .map_err(to_py)?;
        Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
    }

    /// Scan batches of records from virtual position ranges in a BGZF file.
    ///
    /// Parameters
    /// ----------
    /// vpos_ranges : list[tuple[vpos, vpos]]
    ///     List of virtual position ranges as pairs.
    /// columns : list[str], optional
    ///     Names of the columns to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
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
                let fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
                let batch_reader = self
                    .scanner
                    .scan_virtual_ranges(fmt_reader, vpos_ranges, columns, batch_size, limit)
                    .map_err(to_py)?;
                Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
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
    /// Parameters
    /// ----------
    /// region : str
    ///     Genomic range string in the format "chr:start-end",
    ///     "chr:[start,end]" or "chr:[start,end)".
    /// index : path or file-like, optional
    ///     The index file to use for querying the region.
    /// columns : list[str], optional
    ///     Names of the columns to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
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
                let fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
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
                let fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
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
/// region : str
///     Genomic range string in the format "chr:start-end",
///     "chr:[start,end]" or "chr:[start,end)".
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
pub fn read_gtf(
    py: Python,
    src: Py<PyAny>,
    region: Option<String>,
    index: Option<Py<PyAny>>,
    fields: Option<Py<PyAny>>,
    attr_defs: Option<Vec<(String, String)>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let fields = resolve_fields(fields, py)?;
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let scanner =
        GtfScanner::new(None, fields, attr_defs, CoordSystem::OneClosed).map_err(to_py)?;

    let ipc = if let Some(region) = region {
        let region = oxbow::Region::parse(&region, oxbow::CoordSystem::OneClosed).map_err(to_py)?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &src, index)?;
                let batches = scanner
                    .scan_query(fmt_reader, region, index.into_boxed(), None, None, None)
                    .map_err(to_py)?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::gtf::io::Reader::new(bgzf_reader);
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
        let fmt_reader = noodles::gtf::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, None, None, None).map_err(to_py)?;
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
/// region : str
///     Genomic range string in the format "chr:start-end",
///     "chr:[start,end]" or "chr:[start,end)".
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
    src: Py<PyAny>,
    region: Option<String>,
    index: Option<Py<PyAny>>,
    fields: Option<Py<PyAny>>,
    attr_defs: Option<Vec<(String, String)>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let fields = resolve_fields(fields, py)?;
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let scanner =
        GffScanner::new(None, fields, attr_defs, CoordSystem::OneClosed).map_err(to_py)?;

    let ipc = if let Some(region) = region {
        let region = oxbow::Region::parse(&region, oxbow::CoordSystem::OneClosed).map_err(to_py)?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &src, index)?;
                let batches = scanner
                    .scan_query(fmt_reader, region, index.into_boxed(), None, None, None)
                    .map_err(to_py)?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::gff::io::Reader::new(bgzf_reader);
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
        let fmt_reader = noodles::gff::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, None, None, None).map_err(to_py)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

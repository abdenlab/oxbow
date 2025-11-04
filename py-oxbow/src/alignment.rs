use std::io::Seek;
use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use pyo3_arrow::PyRecordBatchReader;
use pyo3_arrow::PySchema;

use noodles::bgzf::io::Seek as _;
use noodles::core::Region;

use crate::util::{
    pyobject_to_bufreader, resolve_cram_index, resolve_fasta_repository, resolve_index, Reader,
};
use oxbow::alignment::{BamScanner, CramScanner, SamScanner};
use oxbow::util::batches_to_ipc;
use oxbow::util::index::IndexType;

/// A SAM file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the SAM file or a file-like object.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
#[pyclass(module = "oxbow.oxbow")]
pub struct PySamScanner {
    src: PyObject,
    reader: Reader,
    scanner: SamScanner,
    compressed: bool,
}

#[pymethods]
impl PySamScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false))]
    fn new(py: Python, src: PyObject, compressed: bool) -> PyResult<Self> {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let mut fmt_reader = noodles::sam::io::Reader::new(reader);
        let header = fmt_reader.read_header()?;
        let reader = fmt_reader.into_inner();
        let scanner = SamScanner::new(header);
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

    /// Discover tag definitions by sniffing `scan_rows` records.
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
    ///     A list of tag definitions, where each definition is a tuple of the
    ///     tag name and the SAM tag type code.
    #[pyo3(signature = (scan_rows=1024))]
    fn tag_defs(&mut self, scan_rows: Option<usize>) -> PyResult<Vec<(String, String)>> {
        let mut reader = self.reader.clone();
        match &mut reader {
            Reader::BgzfFile(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::sam::io::Reader::new(bgzf_reader);
                let defs = self.scanner.tag_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::sam::io::Reader::new(bgzf_reader);
                let defs = self.scanner.tag_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            _ => {
                let pos = reader.stream_position()?;
                let mut fmt_reader = noodles::sam::io::Reader::new(reader);
                let defs = self.scanner.tag_defs(&mut fmt_reader, scan_rows)?;
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
    /// tag_defs : list[tuple[str, str]], optional
    ///    Definitions of tag fields to project.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    #[pyo3(signature = (fields=None, tag_defs=None))]
    fn schema(
        &self,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
    ) -> PyResult<PySchema> {
        let schema = self.scanner.schema(fields, tag_defs)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Scan batches of records from the file.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// tag_defs : list[tuple[str, str]], optional
    ///     Definitions of tag fields to project.
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
    #[pyo3(signature = (fields=None, tag_defs=None, batch_size=1024, limit=None))]
    fn scan(
        &mut self,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let fmt_reader = noodles::sam::io::Reader::new(reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, fields, tag_defs, batch_size, limit)
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
    /// tag_defs : list[tuple[str, str]], optional
    ///     Definitions of tag fields to project.
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
    #[pyo3(signature = (region, index=None, fields=None, tag_defs=None, batch_size=1024, limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn scan_query(
        &mut self,
        py: Python,
        region: String,
        index: Option<PyObject>,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::sam::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &self.src, index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader, region, index, fields, tag_defs, batch_size, limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader, region, index, fields, tag_defs, batch_size, limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::sam::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &self.src, index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader, region, index, fields, tag_defs, batch_size, limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader, region, index, fields, tag_defs, batch_size, limit,
                            )
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

    /// Scan batches of records from the set of unaligned reads.
    ///
    /// This operation requires an index file.
    ///
    /// Parameters
    /// ----------
    /// index : path or file-like, optional
    ///     The index file to use for querying the region. If None and the
    ///     source was provided as a path, we will attempt to load the index
    ///     from the same path with an additional extension.
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// tag_defs : list[tuple[str, str]], optional
    ///     Definitions of tag fields to project.
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
    #[pyo3(signature = (index=None, fields=None, tag_defs=None, batch_size=1024, limit=None))]
    fn scan_unmapped(
        &mut self,
        py: Python,
        index: Option<PyObject>,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::sam::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &self.src, index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_unmapped(fmt_reader, index, fields, tag_defs, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_unmapped(fmt_reader, index, fields, tag_defs, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::sam::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &self.src, index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_unmapped(fmt_reader, index, fields, tag_defs, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_unmapped(fmt_reader, index, fields, tag_defs, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            _ => Err(PyErr::new::<PyValueError, _>(
                "Scanning unmapped reads is only supported for bgzf-compressed sources.",
            )),
        }
    }
}

/// A BAM file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the BAM file or a file-like object.
/// compressed : bool, optional [default: True]
///     Whether the source is BGZF-compressed.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyBamScanner {
    src: PyObject,
    reader: Reader,
    scanner: BamScanner,
    compressed: bool,
}

#[pymethods]
impl PyBamScanner {
    #[new]
    #[pyo3(signature = (src, compressed=true))]
    fn new(py: Python, src: PyObject, compressed: bool) -> PyResult<Self> {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let mut fmt_reader = noodles::bam::io::Reader::from(reader);
        let header = fmt_reader.read_header()?;
        let reader = fmt_reader.into_inner();
        let scanner = BamScanner::new(header);
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

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///    Names of the fixed fields to project.
    /// tag_defs : list[tuple[str, str]], optional
    ///    Definitions of tag fields to project.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    #[pyo3(signature = (fields=None, tag_defs=None))]
    fn schema(
        &self,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
    ) -> PyResult<PySchema> {
        let schema = self.scanner.schema(fields, tag_defs)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Discover tag definitions by sniffing `scan_rows` records.
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
    ///     A list of tag definitions, where each definition is a tuple of the
    ///     tag name and the SAM tag type code.
    #[pyo3(signature = (scan_rows=1024))]
    fn tag_defs(&mut self, scan_rows: Option<usize>) -> PyResult<Vec<(String, String)>> {
        let mut reader = self.reader.clone();
        match &mut reader {
            Reader::BgzfFile(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::bam::io::Reader::from(bgzf_reader);
                let defs = self.scanner.tag_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let pos = bgzf_reader.virtual_position();
                let mut fmt_reader = noodles::bam::io::Reader::from(bgzf_reader);
                let defs = self.scanner.tag_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader.into_inner().seek_to_virtual_position(pos)?;
                Ok(defs)
            }
            _ => {
                let pos = reader.stream_position()?;
                let mut fmt_reader = noodles::bam::io::Reader::from(reader);
                let defs = self.scanner.tag_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader
                    .into_inner()
                    .seek(std::io::SeekFrom::Start(pos))?;
                Ok(defs)
            }
        }
    }

    /// Scan batches of records from the file.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// tag_defs : list[tuple[str, str]], optional
    ///     Definitions of tag fields to project.
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
    #[pyo3(signature = (fields=None, tag_defs=None, batch_size=1024, limit=None))]
    fn scan(
        &mut self,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let fmt_reader = noodles::bam::io::Reader::from(reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, fields, tag_defs, batch_size, limit)
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
    /// tag_defs : list[tuple[str, str]], optional
    ///     Definitions of tag fields to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (region, index=None, fields=None, tag_defs=None, batch_size=1024, limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn scan_query(
        &mut self,
        py: Python,
        region: String,
        index: Option<PyObject>,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bam::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, &self.src, index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader, region, index, fields, tag_defs, batch_size, limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader, region, index, fields, tag_defs, batch_size, limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bam::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, &self.src, index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader, region, index, fields, tag_defs, batch_size, limit,
                            )
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_query(
                                fmt_reader, region, index, fields, tag_defs, batch_size, limit,
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

    /// Scan batches of records from the set of unaligned reads.
    ///
    /// This operation requires an index file.
    ///
    /// Parameters
    /// ----------
    /// index : path or file-like, optional
    ///     The index file to use for querying the region. If None and the
    ///     source was provided as a path, we will attempt to load the index
    ///     from the same path with an additional extension.
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// tag_defs : list[tuple[str, str]], optional
    ///     Definitions of tag fields to project.
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
    #[pyo3(signature = (index=None, fields=None, tag_defs=None, batch_size=1024, limit=None))]
    fn scan_unmapped(
        &mut self,
        py: Python,
        index: Option<PyObject>,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bam::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, &self.src, index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_unmapped(fmt_reader, index, fields, tag_defs, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_unmapped(fmt_reader, index, fields, tag_defs, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bam::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, &self.src, index)?;
                let py_batch_reader = match index {
                    IndexType::Linear(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_unmapped(fmt_reader, index, fields, tag_defs, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                    IndexType::Binned(index) => {
                        let batch_reader = self
                            .scanner
                            .scan_unmapped(fmt_reader, index, fields, tag_defs, batch_size, limit)
                            .map_err(PyErr::new::<PyValueError, _>)?;
                        PyRecordBatchReader::new(Box::new(batch_reader))
                    }
                };
                Ok(py_batch_reader)
            }
            _ => Err(PyErr::new::<PyValueError, _>(
                "Scanning unmapped reads is only supported for indexed bgzf-compressed sources.",
            )),
        }
    }
}

/// A CRAM file scanner.
///
/// Parameters
/// ----------
/// obj : str or file-like
///     The path to the BAM file or a file-like object.
#[pyclass]
pub struct PyCramScanner {
    src: PyObject,
    reader: Reader,
    scanner: CramScanner,
}

#[pymethods]
impl PyCramScanner {
    #[new]
    #[allow(unused_variables)]
    #[pyo3(signature = (src, compressed=None))]
    fn new(py: Python, src: PyObject, compressed: Option<bool>) -> PyResult<Self> {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), false)?;
        let mut fmt_reader = noodles::cram::io::Reader::new(reader);
        let header = fmt_reader.read_header()?;
        let reader = fmt_reader.into_inner();
        let scanner = CramScanner::new(header);
        Ok(Self {
            src,
            reader,
            scanner,
        })
    }

    fn __getstate__(&self, py: Python<'_>) -> PyResult<PyObject> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python<'_>) -> PyResult<(PyObject, PyObject)> {
        let args = (self.src.clone_ref(py),);
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

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///    Names of the fixed fields to project.
    /// tag_defs : list[tuple[str, str]], optional
    ///    Definitions of tag fields to project.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    #[pyo3(signature = (fields=None, tag_defs=None))]
    fn schema(
        &self,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
    ) -> PyResult<PySchema> {
        let schema = self.scanner.schema(fields, tag_defs)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Discover tag definitions by sniffing `scan_rows` records.
    ///
    /// The reader stream is reset to its original position after scanning.
    ///
    /// Parameters
    /// ----------
    /// scan_rows : int, optional [default: 1024]
    ///    The number of records to scan.
    ///
    /// Returns
    /// -------
    /// list[tuple[str, str]]
    ///     A list of tag definitions, where each definition is a tuple of the
    ///     tag name and the SAM tag type code.
    #[pyo3(signature = (scan_rows=1024))]
    fn tag_defs(&mut self, scan_rows: Option<usize>) -> PyResult<Vec<(String, String)>> {
        let mut reader = self.reader.clone();
        match &mut reader {
            Reader::File(_) | Reader::PyFileLike(_) => {
                let pos = reader.stream_position()?;
                let mut fmt_reader = noodles::cram::io::Reader::new(reader);
                let defs = self.scanner.tag_defs(&mut fmt_reader, scan_rows)?;
                fmt_reader
                    .into_inner()
                    .seek(std::io::SeekFrom::Start(pos))?;
                Ok(defs)
            }
            _ => {
                unreachable!()
            }
        }
    }

    /// Scan batches of records from the file.
    ///
    /// Parameters
    /// ----------
    /// reference : path or file-like, optional
    ///     The external reference FASTA file to use for decoding bases in the
    ///     CRAM records. If not provided, sequence bases or references must be
    ///     embedded in the CRAM file.
    /// reference_index : path or file-like, optional
    ///     The index file for the reference FASTA file.
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// tag_defs : list[tuple[str, str]], optional
    ///     Definitions of tag fields to project.
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
    #[pyo3(signature = (reference=None, reference_index=None, fields=None, tag_defs=None, batch_size=1024, limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn scan(
        &mut self,
        py: Python,
        reference: Option<PyObject>,
        reference_index: Option<PyObject>,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        let repo = resolve_fasta_repository(py, reference, reference_index)?;
        let fmt_reader = noodles::cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repo)
            .build_from_reader(reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, fields, tag_defs, batch_size, limit)?;
        let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
        Ok(py_batch_reader)
    }

    /// Scan batches of records from a genomic range.
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
    /// reference : path or file-like, optional
    ///     The external reference FASTA file to use for decoding bases in the
    ///     CRAM records. If not provided, sequence bases or references must be
    ///     embedded in the CRAM file.
    /// reference_index : path or file-like, optional
    ///     The index file for the reference FASTA file.
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// tag_defs : list[tuple[str, str]], optional
    ///     Definitions of tag fields to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (region, index=None, reference=None, reference_index=None, fields=None, tag_defs=None, batch_size=1024, limit=None))]
    #[allow(clippy::too_many_arguments)]
    fn scan_query(
        &mut self,
        py: Python,
        region: String,
        index: Option<PyObject>,
        reference: Option<PyObject>,
        reference_index: Option<PyObject>,
        fields: Option<Vec<String>>,
        tag_defs: Option<Vec<(String, String)>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;
        let index = resolve_cram_index(py, &self.src, index)?;
        let repo = resolve_fasta_repository(py, reference, reference_index)?;

        match self.reader.clone() {
            Reader::File(reader) => {
                let fmt_reader = noodles::cram::io::Reader::new(reader);
                let batch_reader = self.scanner.scan_query(
                    fmt_reader, repo, region, index, fields, tag_defs, batch_size, limit,
                )?;
                let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
                Ok(py_batch_reader)
            }
            Reader::PyFileLike(reader) => {
                let fmt_reader = noodles::cram::io::Reader::new(reader);
                let batch_reader = self.scanner.scan_query(
                    fmt_reader, repo, region, index, fields, tag_defs, batch_size, limit,
                )?;
                let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
                Ok(py_batch_reader)
            }
            _ => {
                unreachable!()
            }
        }
    }
}

/// Return Arrow IPC format from a SAM file.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the source file or a file-like object.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// tag_defs : list[tuple[str, str]], optional
///    Definitions of tag fields to project.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, region=None, index=None, fields=None, tag_defs=None, compressed=false))]
pub fn read_sam(
    py: Python,
    src: PyObject,
    region: Option<String>,
    index: Option<PyObject>,
    fields: Option<Vec<String>>,
    tag_defs: Option<Vec<(String, String)>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let mut fmt_reader = noodles::sam::io::Reader::new(reader);
    let header = fmt_reader.read_header()?;
    let scanner = SamScanner::new(header);
    let reader = fmt_reader.into_inner();

    let ipc = if let Some(region) = region {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::sam::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &src, index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    tag_defs,
                    None,
                    None,
                )?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::sam::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &src, index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    tag_defs,
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
        let fmt_reader = noodles::sam::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, fields, tag_defs, None, None)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

/// Return Arrow IPC format from a BAM file.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the source file or a file-like object.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// tag_defs : list[tuple[str, str]], optional
///    Definitions of tag fields to project.
/// compressed : bool, optional [default: True]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, region=None, index=None, fields=None, tag_defs=None, compressed=true))]
pub fn read_bam(
    py: Python,
    src: PyObject,
    region: Option<String>,
    index: Option<PyObject>,
    fields: Option<Vec<String>>,
    tag_defs: Option<Vec<(String, String)>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let mut fmt_reader = noodles::bam::io::Reader::from(reader);
    let header = fmt_reader.read_header()?;
    let scanner = BamScanner::new(header);
    let reader = fmt_reader.into_inner();

    let ipc = if let Some(region) = region {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bam::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, &src, index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    tag_defs,
                    None,
                    None,
                )?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bam::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, &src, index)?;
                let batches = scanner.scan_query(
                    fmt_reader,
                    region,
                    index.into_boxed(),
                    fields,
                    tag_defs,
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
        let fmt_reader = noodles::bam::io::Reader::from(reader);
        let batches = scanner.scan(fmt_reader, fields, tag_defs, None, None)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

/// Return Arrow IPC format from a CRAM file.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the source file or a file-like object.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// tag_defs : list[tuple[str, str]], optional
///    Definitions of tag fields to project.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, region=None, index=None, reference=None, reference_index=None, fields=None, tag_defs=None))]
#[allow(clippy::too_many_arguments)]
pub fn read_cram(
    py: Python,
    src: PyObject,
    region: Option<String>,
    index: Option<PyObject>,
    reference: Option<PyObject>,
    reference_index: Option<PyObject>,
    fields: Option<Vec<String>>,
    tag_defs: Option<Vec<(String, String)>>,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), false)?;
    let mut fmt_reader = noodles::cram::io::Reader::new(reader);
    let header = fmt_reader.read_header()?;
    let scanner = CramScanner::new(header);
    let reader = fmt_reader.into_inner();

    let repo = resolve_fasta_repository(py, reference, reference_index)?;

    let ipc = if let Some(region) = region {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        match reader {
            Reader::File(reader) => {
                let fmt_reader = noodles::cram::io::Reader::new(reader);
                let index = resolve_cram_index(py, &src, index)?;
                let batches = scanner.scan_query(
                    fmt_reader, repo, region, index, fields, tag_defs, None, None,
                )?;
                batches_to_ipc(batches)
            }
            Reader::PyFileLike(reader) => {
                let fmt_reader = noodles::cram::io::Reader::new(reader);
                let index = resolve_cram_index(py, &src, index)?;
                let batches = scanner.scan_query(
                    fmt_reader, repo, region, index, fields, tag_defs, None, None,
                )?;
                batches_to_ipc(batches)
            }
            _ => {
                unreachable!()
            }
        }
    } else {
        match reader {
            Reader::File(reader) => {
                let fmt_reader = noodles::cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(repo.clone())
                    .build_from_reader(reader);
                let batches = scanner.scan(fmt_reader, fields, tag_defs, None, None)?;
                batches_to_ipc(batches)
            }
            Reader::PyFileLike(reader) => {
                let fmt_reader = noodles::cram::io::reader::Builder::default()
                    .set_reference_sequence_repository(repo.clone())
                    .build_from_reader(reader);
                let batches = scanner.scan(fmt_reader, fields, tag_defs, None, None)?;
                batches_to_ipc(batches)
            }
            _ => {
                unreachable!()
            }
        }
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use pyo3_arrow::PyRecordBatchReader;
use pyo3_arrow::PySchema;

use flate2::read::MultiGzDecoder;
use noodles::bgzf::gzi::Reader as GziReader;
use noodles::bgzf::IndexedReader as IndexedBgzfReader;
use noodles::core::Region;

use crate::util::{pyobject_to_bufreader, resolve_faidx, Reader};
use oxbow::sequence::{FastaScanner, FastqScanner};
use oxbow::util::batches_to_ipc;

/// A FASTQ file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the FASTQ file or a file-like object.
/// compressed : bool, optional [default: False]
///     Whether the source is GZIP-compressed.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyFastqScanner {
    src: PyObject,
    reader: Reader,
    scanner: FastqScanner,
    compressed: bool,
}

#[pymethods]
impl PyFastqScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false))]
    fn new(py: Python, src: PyObject, compressed: bool) -> PyResult<Self> {
        let _src = src.clone_ref(py);
        let reader = pyobject_to_bufreader(py, src, false)?;
        let scanner = FastqScanner::new();
        Ok(Self {
            src: _src,
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

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///    Names of the fixed fields to project.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    #[pyo3(signature = (fields=None))]
    fn schema(&self, fields: Option<Vec<String>>) -> PyResult<PySchema> {
        let schema = self.scanner.schema(fields)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Scan the source as record batches.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
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
        if self.compressed {
            // Decode a fastq.gz using a conventional gzip stream decoder.
            let gz_reader = std::io::BufReader::new(MultiGzDecoder::new(reader));
            let fmt_reader = noodles::fastq::io::Reader::new(gz_reader);
            let batch_reader = self
                .scanner
                .scan(fmt_reader, fields, batch_size, limit)
                .map_err(PyErr::new::<PyValueError, _>)?;
            Ok(PyRecordBatchReader::new(Box::new(batch_reader)))
        } else {
            let fmt_reader = noodles::fastq::io::Reader::new(reader);
            let batch_reader = self
                .scanner
                .scan(fmt_reader, fields, batch_size, limit)
                .map_err(PyErr::new::<PyValueError, _>)?;
            Ok(PyRecordBatchReader::new(Box::new(batch_reader)))
        }
    }
}

/// A FASTA file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the FASTA file or a file-like object.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyFastaScanner {
    src: PyObject,
    reader: Reader,
    scanner: FastaScanner,
    compressed: bool,
}

#[pymethods]
impl PyFastaScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false))]
    fn new(py: Python, src: PyObject, compressed: bool) -> PyResult<Self> {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), false)?;
        let scanner = FastaScanner::new();
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

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///    Names of the fixed fields to project.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    #[pyo3(signature = (fields=None))]
    fn schema(&self, fields: Option<Vec<String>>) -> PyResult<PySchema> {
        let schema = self.scanner.schema(fields)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Scan the source as record batches.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// batch_size : int, optional [default: 1]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan. If None, records are scanned
    ///     until EOF.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    ///
    /// Notes
    /// -----
    /// Since reference sequences are often large, the default batch size is
    /// set to 1.
    #[pyo3(signature = (fields=None, batch_size=1, limit=None))]
    fn scan(
        &mut self,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        if self.compressed {
            let gz_reader = std::io::BufReader::new(MultiGzDecoder::new(reader));
            let fmt_reader = noodles::fasta::io::Reader::new(gz_reader);
            let batch_reader = self
                .scanner
                .scan(fmt_reader, fields, batch_size, limit)
                .map_err(PyErr::new::<PyValueError, _>)?;
            Ok(PyRecordBatchReader::new(Box::new(batch_reader)))
        } else {
            let fmt_reader = noodles::fasta::io::Reader::new(reader);
            let batch_reader = self
                .scanner
                .scan(fmt_reader, fields, batch_size, limit)
                .map_err(PyErr::new::<PyValueError, _>)?;
            Ok(PyRecordBatchReader::new(Box::new(batch_reader)))
        }
    }

    /// Scan sequence slices as record batches from a list of genomic ranges.
    ///
    /// This operation requires an index file.
    ///
    /// Parameters
    /// ----------
    /// regions : list[str]
    ///     Genomic ranges in the format "chr:start-end".
    /// index : path or file-like, optional
    ///     The FAI index file to use for slicing the reference sequences.
    ///     If None and the source was provided as a path, we will attempt to
    ///     load the index from the same path with an additional extension.
    /// gzi : path or file-like, optional
    ///     A GZI index file to use if the source is BGZF-encoded.
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
    ///
    /// Notes
    /// -----
    /// An FAI index is required to slice the reference sequences.
    /// If the source is BGZF-compressed, an additional GZI index is also
    /// required. The GZI index is used to translate uncompressed positions
    /// (from the FAI index) into compressed positions in the BGZF file.
    #[pyo3(signature = (regions, index=None, gzi=None, fields=None, batch_size=1024))]
    fn scan_query(
        &mut self,
        py: Python,
        regions: Vec<String>,
        index: Option<PyObject>,
        gzi: Option<PyObject>,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        // Load FAI index.
        let index = resolve_faidx(py, self.src.clone_ref(py), index)?;
        let reader = self.reader.clone();

        // Parse the genomic ranges.
        let regions: Vec<Region> = regions
            .into_iter()
            .map(|s| {
                s.parse::<Region>()
                    .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;

        if self.compressed {
            // A BGZF IndexedReader seeks within a compressed file using uncompressed positions.
            let gzi_source = match gzi {
                Some(gzi) => pyobject_to_bufreader(py, gzi.clone_ref(py), false)?,
                None => {
                    return Err(PyErr::new::<PyValueError, _>(
                        "GZI index required to query BGZF-encoded FASTA.",
                    ));
                }
            };
            let gzindex = GziReader::new(gzi_source).read_index()?;
            let bgzf_reader = IndexedBgzfReader::new(reader, gzindex);
            let fmt_reader = noodles::fasta::io::Reader::new(bgzf_reader);
            let batch_reader = self
                .scanner
                .scan_query(fmt_reader, regions, index, fields, batch_size)
                .map_err(PyErr::new::<PyValueError, _>)?;
            Ok(PyRecordBatchReader::new(Box::new(batch_reader)))
        } else {
            let fmt_reader = noodles::fasta::io::Reader::new(reader);
            let batch_reader = self
                .scanner
                .scan_query(fmt_reader, regions, index, fields, batch_size)
                .map_err(PyErr::new::<PyValueError, _>)?;
            Ok(PyRecordBatchReader::new(Box::new(batch_reader)))
        }
    }
}

/// Return Arrow IPC format from a FASTQ file.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the source file or a file-like object.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// compressed : bool, optional [default: False]
///     Whether the source is GZIP-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, fields=None, compressed=false))]
pub fn read_fastq(
    py: Python,
    src: PyObject,
    fields: Option<Vec<String>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src, false)?;
    let scanner = FastqScanner::new();

    let ipc = if compressed {
        let gz_reader = std::io::BufReader::new(MultiGzDecoder::new(reader));
        let fmt_reader = noodles::fastq::io::Reader::new(gz_reader);
        let batches = scanner.scan(fmt_reader, fields, None, None)?;
        batches_to_ipc(batches)
    } else {
        let fmt_reader = noodles::fastq::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, fields, None, None)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

/// Return Arrow IPC format from a FASTA file.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the source file or a file-like object.
/// regions : list[str], optional
///     Genomic ranges in the format "chr:start-end".
/// index : path or file-like, optional
///     The FAI index file to use for slicing the reference sequences.
///     If None and the source was provided as a path, we will attempt to
///     load the index from the same path with an additional extension.
/// gzi : path or file-like, optional
///     A GZI index file to use if the source is BGZF-encoded.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, regions=None, index=None, gzi=None, fields=None, compressed=false))]
pub fn read_fasta(
    py: Python,
    src: PyObject,
    regions: Option<Vec<String>>,
    index: Option<PyObject>,
    gzi: Option<PyObject>,
    fields: Option<Vec<String>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let scanner = FastaScanner::new();

    let ipc = if let Some(regions) = regions {
        let index = resolve_faidx(py, src, index)?;
        let regions: Vec<Region> = regions
            .into_iter()
            .map(|s| {
                s.parse::<Region>()
                    .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;
        if compressed {
            let gzi_source = match gzi {
                Some(gzi) => pyobject_to_bufreader(py, gzi.clone_ref(py), false)?,
                None => {
                    return Err(PyErr::new::<PyValueError, _>(
                        "GZI index required to query BGZF-encoded FASTA.",
                    ));
                }
            };
            let gzindex = GziReader::new(gzi_source).read_index()?;
            let bgzf_reader = IndexedBgzfReader::new(reader, gzindex);
            let fmt_reader = noodles::fasta::io::Reader::new(bgzf_reader);
            let batches = scanner
                .scan_query(fmt_reader, regions, index, fields, None)
                .map_err(PyErr::new::<PyValueError, _>)?;
            batches_to_ipc(batches)
        } else {
            let fmt_reader = noodles::fasta::io::Reader::new(reader);
            let batches = scanner
                .scan_query(fmt_reader, regions, index, fields, None)
                .map_err(PyErr::new::<PyValueError, _>)?;
            batches_to_ipc(batches)
        }
    } else {
        let fmt_reader = noodles::fasta::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, fields, None, None)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

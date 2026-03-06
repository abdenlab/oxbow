use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use pyo3_arrow::PyRecordBatchReader;
use pyo3_arrow::PySchema;

use flate2::read::MultiGzDecoder;
use noodles::bgzf::gzi::io::Reader as GziReader;
use noodles::bgzf::io::IndexedReader as IndexedBgzfReader;
use noodles::core::Region;

use crate::error::{err_on_unwind, to_py};
use crate::util::{pyobject_to_bufreader, resolve_faidx, PyVirtualPosition, Reader};
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
/// fields : list[str], optional
///     Names of the fixed fields to project.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyFastqScanner {
    src: Py<PyAny>,
    reader: Reader,
    scanner: FastqScanner,
    compressed: bool,
    fields: Option<Vec<String>>,
}

#[pymethods]
impl PyFastqScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false, fields=None))]
    fn new(
        py: Python,
        src: Py<PyAny>,
        compressed: bool,
        fields: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let _src = src.clone_ref(py);
        let reader = pyobject_to_bufreader(py, src, false)?;
        let scanner = FastqScanner::new(fields.clone()).map_err(to_py)?;
        Ok(Self {
            src: _src,
            reader,
            scanner,
            compressed,
            fields,
        })
    }

    fn __getstate__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python<'_>) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
        let args = (self.src.clone_ref(py), self.compressed.into_py_any(py)?);
        let kwargs = PyDict::new(py);
        if let Some(ref fields) = self.fields {
            kwargs.set_item("fields", fields)?;
        }
        Ok((args.into_py_any(py)?, kwargs.into_py_any(py)?))
    }

    /// Return the names of the fixed fields.
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

    /// Scan the source as record batches.
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
        if self.compressed {
            let gz_reader = std::io::BufReader::new(MultiGzDecoder::new(reader));
            let fmt_reader = noodles::fastq::io::Reader::new(gz_reader);
            let batch_reader = self
                .scanner
                .scan(fmt_reader, columns, batch_size, limit)
                .map_err(to_py)?;
            Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
        } else {
            let fmt_reader = noodles::fastq::io::Reader::new(reader);
            let batch_reader = self
                .scanner
                .scan(fmt_reader, columns, batch_size, limit)
                .map_err(to_py)?;
            Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
        }
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
        let fmt_reader = noodles::fastq::io::Reader::new(reader);
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
        if !self.compressed {
            return Err(PyErr::new::<PyValueError, _>(
                "scan_virtual_ranges requires a BGZF-compressed source.",
            ));
        }
        let vpos_ranges = vpos_ranges
            .into_iter()
            .map(|(start, end)| Ok((start.to_virtual_position()?, end.to_virtual_position()?)))
            .collect::<PyResult<Vec<_>>>()?;
        let reader = self.reader.clone();
        let bgzf_reader = noodles::bgzf::io::Reader::new(reader);
        let fmt_reader = noodles::fastq::io::Reader::new(bgzf_reader);
        let batch_reader = self
            .scanner
            .scan_virtual_ranges(fmt_reader, vpos_ranges, columns, batch_size, limit)
            .map_err(to_py)?;
        Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
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
/// fields : list[str], optional
///     Names of the fixed fields to project.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyFastaScanner {
    src: Py<PyAny>,
    reader: Reader,
    scanner: FastaScanner,
    compressed: bool,
    fields: Option<Vec<String>>,
}

#[pymethods]
impl PyFastaScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false, fields=None))]
    fn new(
        py: Python,
        src: Py<PyAny>,
        compressed: bool,
        fields: Option<Vec<String>>,
    ) -> PyResult<Self> {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), false)?;
        let scanner = FastaScanner::new(fields.clone()).map_err(to_py)?;
        Ok(Self {
            src,
            reader,
            scanner,
            compressed,
            fields,
        })
    }

    fn __getstate__(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python<'_>) -> PyResult<(Py<PyAny>, Py<PyAny>)> {
        let args = (self.src.clone_ref(py), self.compressed.into_py_any(py)?);
        let kwargs = PyDict::new(py);
        if let Some(ref fields) = self.fields {
            kwargs.set_item("fields", fields)?;
        }
        Ok((args.into_py_any(py)?, kwargs.into_py_any(py)?))
    }

    /// Return the names of the fixed fields.
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

    /// Scan the source as record batches.
    ///
    /// Parameters
    /// ----------
    /// columns : list[str], optional
    ///     Names of the columns to project.
    /// batch_size : int, optional [default: 1]
    ///     The number of records to include in each batch.
    /// limit : int, optional
    ///     The maximum number of records to scan. If None, records are scanned
    ///     until EOF.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///
    /// Notes
    /// -----
    /// Since reference sequences are often large, the default batch size is
    /// set to 1.
    #[pyo3(signature = (columns=None, batch_size=1, limit=None))]
    fn scan(
        &mut self,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let reader = self.reader.clone();
        if self.compressed {
            let gz_reader = std::io::BufReader::new(MultiGzDecoder::new(reader));
            let fmt_reader = noodles::fasta::io::Reader::new(gz_reader);
            let batch_reader = self
                .scanner
                .scan(fmt_reader, columns, batch_size, limit)
                .map_err(to_py)?;
            Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
        } else {
            let fmt_reader = noodles::fasta::io::Reader::new(reader);
            let batch_reader = self
                .scanner
                .scan(fmt_reader, columns, batch_size, limit)
                .map_err(to_py)?;
            Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
        }
    }

    /// Scan sequence slices as record batches from a list of genomic ranges.
    ///
    /// Parameters
    /// ----------
    /// regions : list[str]
    ///     Genomic ranges in the format "chr:start-end".
    /// index : path or file-like, optional
    ///     The FAI index file.
    /// gzi : path or file-like, optional
    ///     A GZI index file for BGZF-encoded sources.
    /// columns : list[str], optional
    ///     Names of the columns to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    #[pyo3(signature = (regions, index=None, gzi=None, columns=None, batch_size=1024))]
    fn scan_query(
        &mut self,
        py: Python,
        regions: Vec<String>,
        index: Option<Py<PyAny>>,
        gzi: Option<Py<PyAny>>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let index = resolve_faidx(py, &self.src, index)?;
        let reader = self.reader.clone();

        let regions: Vec<Region> = regions
            .into_iter()
            .map(|s| {
                s.parse::<Region>()
                    .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
            })
            .collect::<Result<Vec<_>, _>>()?;

        if self.compressed {
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
                .scan_query(fmt_reader, regions, index, columns, batch_size)
                .map_err(to_py)?;
            Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
        } else {
            let fmt_reader = noodles::fasta::io::Reader::new(reader);
            let batch_reader = self
                .scanner
                .scan_query(fmt_reader, regions, index, columns, batch_size)
                .map_err(to_py)?;
            Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
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
    src: Py<PyAny>,
    fields: Option<Vec<String>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src, false)?;
    let scanner = FastqScanner::new(fields).map_err(to_py)?;

    let ipc = if compressed {
        let gz_reader = std::io::BufReader::new(MultiGzDecoder::new(reader));
        let fmt_reader = noodles::fastq::io::Reader::new(gz_reader);
        let batches = scanner.scan(fmt_reader, None, None, None).map_err(to_py)?;
        batches_to_ipc(batches)
    } else {
        let fmt_reader = noodles::fastq::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, None, None, None).map_err(to_py)?;
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
///     The FAI index file.
/// gzi : path or file-like, optional
///     A GZI index file for BGZF-encoded sources.
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
    src: Py<PyAny>,
    regions: Option<Vec<String>>,
    index: Option<Py<PyAny>>,
    gzi: Option<Py<PyAny>>,
    fields: Option<Vec<String>>,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let scanner = FastaScanner::new(fields).map_err(to_py)?;

    let ipc = if let Some(regions) = regions {
        let index = resolve_faidx(py, &src, index)?;
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
                .scan_query(fmt_reader, regions, index, None, None)
                .map_err(to_py)?;
            batches_to_ipc(batches)
        } else {
            let fmt_reader = noodles::fasta::io::Reader::new(reader);
            let batches = scanner
                .scan_query(fmt_reader, regions, index, None, None)
                .map_err(to_py)?;
            batches_to_ipc(batches)
        }
    } else {
        let fmt_reader = noodles::fasta::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, None, None, None).map_err(to_py)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

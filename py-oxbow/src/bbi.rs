use std::io::Seek;
use std::sync::Arc;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use pyo3::IntoPyObjectExt;
use pyo3_arrow::PyRecordBatchReader;
use pyo3_arrow::PySchema;

use bigtools::bed::autosql::parse::parse_autosql;
use noodles::core::Region;

use crate::util::{pyobject_to_bufreader, Reader};
use oxbow::bbi::model::base::field::FieldDef;
use oxbow::bbi::{BBIReader, BBIZoomScanner, BedSchema, BigBedScanner, BigWigScanner};
use oxbow::util::batches_to_ipc;

#[pyclass(eq, eq_int, module = "oxbow.oxbow")]
#[derive(Clone, PartialEq)]
pub enum PyBBIFileType {
    BigWig,
    BigBed,
}

/// A BigWig file scanner.
///
/// Parameters
/// ----------
/// obj : str or file-like
///     The path to the BigWig file or a file-like object.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyBigWigScanner {
    _src: PyObject,
    reader: Reader,
    scanner: BigWigScanner,
}

#[pymethods]
impl PyBigWigScanner {
    #[new]
    #[pyo3(signature = (src))]
    fn new(py: Python, src: PyObject) -> PyResult<Self> {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), false)?;
        let fmt_reader = bigtools::BigWigRead::open(reader).unwrap();
        let info = fmt_reader.info().clone();
        let reader = fmt_reader.into_inner();
        let scanner = BigWigScanner::new(info);
        Ok(Self {
            _src: src,
            reader,
            scanner,
        })
    }

    fn __getstate__(&self, py: Python) -> PyResult<PyObject> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python) -> PyResult<(PyObject, PyObject)> {
        let args = (self._src.clone_ref(py),);
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

    /// Return the names of the fields.
    fn field_names(&self) -> Vec<String> {
        self.scanner.field_names()
    }

    /// Return the zoom/reduction level resolutions.
    fn zoom_levels(&self) -> Vec<u32> {
        self.scanner.zoom_levels()
    }

    /// Return a scanner for a specific zoom level.
    ///
    /// Parameters
    /// ---------
    /// zoom_level : int
    ///     The resolution (in bp) of the zoom level to scan.
    ///
    /// Returns
    /// -------
    /// PyBBIZoomScanner
    ///     A scanner for the specified zoom level.
    fn get_zoom(&mut self, zoom_level: u32) -> PyResult<PyBBIZoomScanner> {
        Python::with_gil(|py| {
            let py_zoom = PyBBIZoomScanner::new(
                py,
                self._src.clone_ref(py),
                PyBBIFileType::BigWig,
                zoom_level,
            );
            Ok(py_zoom)
        })
    }

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fields to project.
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
        let info = self.scanner.info().clone();
        let fmt_reader = bigtools::BigWigRead::with_info(info, reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, fields, batch_size, limit)
            .map_err(PyErr::new::<PyValueError, _>)?;
        let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
        Ok(py_batch_reader)
    }

    /// Scan batches of records from a genomic range query.
    ///
    /// Parameters
    /// ----------
    /// region : str
    ///     Genomic region in the format "chr:start-end".
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (region, fields=None, batch_size=1024, limit=None))]
    fn scan_query(
        &mut self,
        region: String,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        let reader = self.reader.clone();
        let info = self.scanner.info().clone();
        let fmt_reader = bigtools::BigWigRead::with_info(info, reader);
        let batch_reader = self
            .scanner
            .scan_query(fmt_reader, region, fields, batch_size, limit)
            .map_err(PyErr::new::<PyValueError, _>)?;
        let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
        Ok(py_batch_reader)
    }
}

/// A BigBed file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the BigBed file or a file-like object.
/// schema : str, optional
///     The BED schema to use for parsing BigBed records. May be a
///     ``bed[m[+[n]]]`` string, ``"bedgraph"``, or ``"autosql"``. If not
///     specified, the BigBed is interpreted as BED3+, where all ancillary
///     fields are returned as a single lumped string field named "rest".
///     If "autosql", the file's AutoSql definition is used to parse the
///     records, if it exists.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyBigBedScanner {
    _src: PyObject,
    _schema: Option<String>,
    reader: Reader,
    scanner: BigBedScanner,
}

#[pymethods]
impl PyBigBedScanner {
    #[new]
    #[pyo3(signature = (src, schema="bed3+"))]
    fn new(py: Python, src: PyObject, schema: Option<&str>) -> PyResult<Self> {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), false)?;
        let mut fmt_reader = bigtools::BigBedRead::open(reader).unwrap();
        let bed_schema = match schema {
            Some("autosql") => {
                let autosql_str = fmt_reader.autosql().unwrap().unwrap();
                let mut declarations = parse_autosql(&autosql_str).unwrap();
                if declarations.len() > 1 {
                    return Err(PyErr::new::<PyValueError, _>(
                        "Autosql schema must contain exactly one declaration.",
                    ));
                }
                let declaration = declarations.pop().unwrap();
                let fields: Vec<_> = declaration
                    .fields
                    .iter()
                    .skip(3)
                    .map(FieldDef::try_from)
                    .collect::<Result<_, _>>()?;
                BedSchema::new(3, Some(fields))
            }
            Some(schema) => schema.parse(),
            None => BedSchema::new_from_nm(3, None),
        }?;
        let info = fmt_reader.info().clone();
        let reader = fmt_reader.into_inner();
        let scanner = BigBedScanner::new(bed_schema, info);
        let _schema: Option<String> = schema.map(|schema| schema.to_string());
        // let schema = schema.map(|s| s.to_owned());
        Ok(Self {
            _src: src,
            _schema,
            reader,
            scanner,
        })
    }

    fn __getstate__(&self, py: Python) -> PyResult<PyObject> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python) -> PyResult<(PyObject, PyObject)> {
        let args = (
            self._src.clone_ref(py),
            self._schema.clone().into_py_any(py)?,
        );
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

    /// Return the names of the fields.
    fn field_names(&self) -> Vec<String> {
        self.scanner.field_names()
    }

    /// Return the raw autosql schema definition.
    ///
    /// Returns
    /// -------
    /// str
    fn read_autosql(&mut self) -> PyResult<String> {
        self.reader.seek(std::io::SeekFrom::Start(0)).unwrap();
        let reader = self.reader.clone();
        let mut fmt_reader = bigtools::BigBedRead::open(reader).unwrap();
        let autosql_str = fmt_reader.autosql().unwrap().unwrap();
        Ok(autosql_str)
    }

    /// Return the zoom/reduction level resolutions.
    fn zoom_levels(&self) -> Vec<u32> {
        self.scanner.zoom_levels()
    }

    /// Return a scanner for a specific zoom level.
    ///
    /// Parameters
    /// ---------
    /// zoom_level : int
    ///     The resolution (in bp) of the zoom level to scan.
    ///
    /// Returns
    /// -------
    /// PyBBIZoomScanner
    ///     A scanner for the specified zoom level.
    fn get_zoom(&mut self, zoom_level: u32) -> PyResult<PyBBIZoomScanner> {
        Python::with_gil(|py| {
            let py_zoom = PyBBIZoomScanner::new(
                py,
                self._src.clone_ref(py),
                PyBBIFileType::BigBed,
                zoom_level,
            );
            Ok(py_zoom)
        })
    }

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fields to project.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    #[pyo3(signature = (fields=None))]
    fn schema(&self, fields: Option<Vec<String>>) -> PyResult<PySchema> {
        let schema = self.scanner.schema(fields)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Scan batches of records from the file's base-level records.
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
        let info = self.scanner.info().clone();
        let fmt_reader = bigtools::BigBedRead::with_info(info, reader);
        let batch_reader = self
            .scanner
            .scan(fmt_reader, fields, batch_size, limit)
            .map_err(PyErr::new::<PyValueError, _>)?;
        let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
        Ok(py_batch_reader)
    }

    /// Scan batches of records from a genomic range query.
    ///
    /// Parameters
    /// ----------
    /// region : str
    ///     Genomic region in the format "chr:start-end".
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (region, fields=None, batch_size=1024, limit=None))]
    fn scan_query(
        &mut self,
        region: String,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        let reader = self.reader.clone();
        let info = self.scanner.info().clone();
        let fmt_reader = bigtools::BigBedRead::with_info(info, reader);
        let batch_reader = self
            .scanner
            .scan_query(fmt_reader, region, fields, batch_size, limit)
            .map_err(PyErr::new::<PyValueError, _>)?;
        let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
        Ok(py_batch_reader)
    }
}

/// A BBI file zoom level scanner.
///
/// Can only be initialized from a BigBed or BigWig scanner.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyBBIZoomScanner {
    src: PyObject,
    reader: Reader,
    bbi_type: PyBBIFileType,
    zoom_level: u32,
    scanner: BBIZoomScanner,
}

#[pymethods]
impl PyBBIZoomScanner {
    #[new]
    pub fn new(py: Python, src: PyObject, bbi_type: PyBBIFileType, zoom_level: u32) -> Self {
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), false)
            .expect("Failed to convert PyObject to BufReader");
        match bbi_type {
            PyBBIFileType::BigBed => {
                let fmt_reader = bigtools::BigBedRead::open(reader).unwrap();
                let ref_names = fmt_reader
                    .chroms()
                    .iter()
                    .map(|info| info.name.to_string())
                    .collect();
                let zoom_levels: Vec<u32> = fmt_reader
                    .info()
                    .zoom_headers
                    .iter()
                    .map(|header| header.reduction_level)
                    .collect();
                if !zoom_levels.contains(&zoom_level) {
                    panic!(
                        "Zoom resolution {} not found in BBI file. Available reduction levels: {:?}.",
                        zoom_level, zoom_levels
                    );
                }
                let reader = fmt_reader.into_inner();
                let scanner = BBIZoomScanner::new(ref_names, zoom_level);
                Self {
                    src,
                    reader,
                    bbi_type,
                    zoom_level,
                    scanner,
                }
            }
            PyBBIFileType::BigWig => {
                let fmt_reader = bigtools::BigWigRead::open(reader).unwrap();
                let ref_names = fmt_reader
                    .chroms()
                    .iter()
                    .map(|info| info.name.to_string())
                    .collect();
                let zoom_levels: Vec<u32> = fmt_reader
                    .info()
                    .zoom_headers
                    .iter()
                    .map(|header| header.reduction_level)
                    .collect();
                if !zoom_levels.contains(&zoom_level) {
                    panic!(
                        "Zoom resolution {} not found in BBI file. Available reduction levels: {:?}.",
                        zoom_level, zoom_levels
                    );
                }
                let reader = fmt_reader.into_inner();
                let scanner = BBIZoomScanner::new(ref_names, zoom_level);
                Self {
                    src,
                    reader,
                    bbi_type,
                    zoom_level,
                    scanner,
                }
            }
        }
    }

    fn __getstate__(&self, py: Python) -> PyResult<PyObject> {
        Ok(py.None())
    }

    fn __getnewargs_ex__(&self, py: Python) -> PyResult<(PyObject, PyObject)> {
        let args = (
            self.src.clone_ref(py),
            self.bbi_type.clone().into_py_any(py)?,
            self.zoom_level.into_py_any(py)?,
        );
        let kwargs = PyDict::new(py);
        Ok((args.into_py_any(py)?, kwargs.into_py_any(py)?))
    }

    /// Return the names of the reference sequences.
    fn field_names(&self) -> Vec<String> {
        self.scanner.field_names()
    }

    /// Return the Arrow schema.
    ///
    /// Parameters
    /// ----------
    /// fields : list[str], optional
    ///     Names of the fields to project.
    ///
    /// Returns
    /// -------
    /// arro3 Schema (pycapsule)
    #[pyo3(signature = (fields=None))]
    fn schema(&self, fields: Option<Vec<String>>) -> PyResult<PySchema> {
        let schema = self.scanner.schema(fields)?;
        Ok(PySchema::new(Arc::new(schema)))
    }

    /// Scan batches of records from the zoom level.
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
        self.reader.seek(std::io::SeekFrom::Start(0)).unwrap();
        let reader = self.reader.clone();
        match self.bbi_type {
            PyBBIFileType::BigBed => {
                let fmt_reader = bigtools::BigBedRead::open(reader).unwrap();
                let reader = BBIReader::BigBed(fmt_reader);
                let batch_reader = self
                    .scanner
                    .scan(reader, fields, batch_size, limit)
                    .map_err(PyErr::new::<PyValueError, _>)?;
                let py_batch_reader = PyRecordBatchReader::new(batch_reader);
                Ok(py_batch_reader)
            }
            PyBBIFileType::BigWig => {
                let fmt_reader = bigtools::BigWigRead::open(reader).unwrap();
                let reader = BBIReader::BigWig(fmt_reader);
                let batch_reader = self
                    .scanner
                    .scan(reader, fields, batch_size, limit)
                    .map_err(PyErr::new::<PyValueError, _>)?;
                let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
                Ok(py_batch_reader)
            }
        }
    }

    /// Scan batches of records from a genomic range query.
    ///
    /// Parameters
    /// ----------
    /// region : str
    ///     Genomic region in the format "chr:start-end".
    /// fields : list[str], optional
    ///     Names of the fixed fields to project.
    /// batch_size : int, optional [default: 1024]
    ///     The number of records to include in each batch.
    ///
    /// Returns
    /// -------
    /// arro3 RecordBatchReader (pycapsule)
    ///     An iterator yielding Arrow record batches.
    #[pyo3(signature = (region, fields=None, batch_size=1024, limit=None))]
    fn scan_query(
        &mut self,
        region: String,
        fields: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        self.reader.seek(std::io::SeekFrom::Start(0)).unwrap();
        let reader = self.reader.clone();
        match self.bbi_type {
            PyBBIFileType::BigBed => {
                let fmt_reader = bigtools::BigBedRead::open(reader).unwrap();
                let reader = BBIReader::BigBed(fmt_reader);
                let batch_reader = self
                    .scanner
                    .scan_query(reader, region, fields, batch_size, limit)
                    .map_err(PyErr::new::<PyValueError, _>)?;
                let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
                Ok(py_batch_reader)
            }
            PyBBIFileType::BigWig => {
                let fmt_reader = bigtools::BigWigRead::open(reader).unwrap();
                let reader = BBIReader::BigWig(fmt_reader);
                let batch_reader = self
                    .scanner
                    .scan_query(reader, region, fields, batch_size, limit)
                    .map_err(PyErr::new::<PyValueError, _>)?;
                let py_batch_reader = PyRecordBatchReader::new(Box::new(batch_reader));
                Ok(py_batch_reader)
            }
        }
    }
}

/// Return Arrow IPC format from a BigWig file.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the source file or a file-like object.
/// fields : list[str], optional
///     Names of the fixed fields to project.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, region=None, fields=None))]
pub fn read_bigwig(
    py: Python,
    src: PyObject,
    region: Option<String>,
    fields: Option<Vec<String>>,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), false)?;

    let ipc = if let Some(region) = region {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        let fmt_reader = bigtools::BigWigRead::open(reader)
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;
        let info = fmt_reader.info().clone();
        let scanner = BigWigScanner::new(info);
        let batches = scanner
            .scan_query(fmt_reader, region, fields, None, None)
            .map_err(PyErr::new::<PyValueError, _>)?;
        batches_to_ipc(batches)
    } else {
        let fmt_reader = bigtools::BigWigRead::open(reader)
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;
        let info = fmt_reader.info().clone();
        let scanner = BigWigScanner::new(info);
        let batches = scanner
            .scan(fmt_reader, fields, None, None)
            .map_err(PyErr::new::<PyValueError, _>)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

/// Return Arrow IPC format from a BigBed file.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the source file or a file-like object.
/// bed_schema : str
///     The BED schema to use for parsing BigBed records.
/// fields : list[str], optional
///     Names of the fixed fields to project.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, bed_schema="bed3+", region=None, fields=None))]
pub fn read_bigbed(
    py: Python,
    src: PyObject,
    bed_schema: &str,
    region: Option<String>,
    fields: Option<Vec<String>>,
) -> PyResult<Vec<u8>> {
    let bed_schema = bed_schema.parse()?;
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), false)?;

    let ipc = if let Some(region) = region {
        let region = region
            .parse::<Region>()
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

        let fmt_reader = bigtools::BigBedRead::open(reader)
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;
        let info = fmt_reader.info().clone();
        let scanner = BigBedScanner::new(bed_schema, info);
        let batches = scanner
            .scan_query(fmt_reader, region, fields, None, None)
            .map_err(PyErr::new::<PyValueError, _>)?;
        batches_to_ipc(batches)
    } else {
        let fmt_reader = bigtools::BigBedRead::open(reader)
            .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;
        let info = fmt_reader.info().clone();
        let scanner = BigBedScanner::new(bed_schema, info);
        let batches = scanner
            .scan(fmt_reader, fields, None, None)
            .map_err(PyErr::new::<PyValueError, _>)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

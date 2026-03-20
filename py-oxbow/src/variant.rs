use std::sync::Arc;

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
use oxbow::util::batches_to_ipc;
use oxbow::util::index::IndexType;
use oxbow::variant::{BcfScanner, GenotypeBy, VcfScanner};
use oxbow::CoordSystem;

/// A VCF file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///    The path to the VCF file or a file-like object.
/// compressed : bool, optional [default: False]
///    Whether the source is BGZF-compressed.
/// fields : list[str], optional
///    Names of the fixed fields to project.
/// info_fields : list[str], optional
///    Names of the INFO fields to project.
/// genotype_fields : list[str], optional
///    Names of the sample-specific genotype fields to project.
/// genotype_by : Literal["sample", "field"], optional [default: "sample"]
///    How to project the genotype fields. If "sample", the columns
///    correspond to the samples. If "field", the columns correspond to
///    the genotype fields.
/// samples : list[str] or None, optional [default: None]
///    Names of the samples to include in the genotype output. ``"*"`` for
///    all samples, a list to select specific samples, or ``None`` to omit
///    all sample genotype data.
/// samples_nested : bool, optional [default: False]
///   Whether to nest sample genotype data under a single ``"samples"`` struct
///   column.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyVcfScanner {
    src: Py<PyAny>,
    reader: Reader,
    scanner: VcfScanner,
    compressed: bool,
}

#[pymethods]
impl PyVcfScanner {
    #[new]
    #[pyo3(signature = (src, compressed=false, fields=None, info_fields=None, genotype_fields=None, genotype_by=None, samples=None, samples_nested=false, coords=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python,
        src: Py<PyAny>,
        compressed: bool,
        fields: Option<Py<PyAny>>,
        info_fields: Option<Py<PyAny>>,
        genotype_fields: Option<Py<PyAny>>,
        genotype_by: Option<String>,
        samples: Option<Py<PyAny>>,
        samples_nested: bool,
        coords: Option<String>,
    ) -> PyResult<Self> {
        let coord_system = resolve_coord_system(coords)?;
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let mut fmt_reader = noodles::vcf::io::Reader::new(reader);
        let header = fmt_reader.read_header()?;
        let reader = fmt_reader.into_inner();
        let gt_by = resolve_genotype_by(genotype_by)?;
        let scanner = VcfScanner::new(
            header,
            resolve_fields(fields, py)?,
            resolve_fields(info_fields, py)?,
            resolve_fields(genotype_fields, py)?,
            gt_by,
            resolve_fields(samples, py)?,
            Some(samples_nested),
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
        if let Some(defs) = model.info_defs() {
            let names: Vec<&str> = defs.iter().map(|d| d.name.as_str()).collect();
            kwargs.set_item("info_fields", names)?;
        }
        if let Some(defs) = model.genotype_defs() {
            let names: Vec<&str> = defs.iter().map(|d| d.name.as_str()).collect();
            kwargs.set_item("genotype_fields", names)?;
        }
        if let Some(samples) = model.samples() {
            kwargs.set_item("samples", samples)?;
        }
        let gt_by = match model.genotype_by() {
            GenotypeBy::Sample => "sample",
            GenotypeBy::Field => "field",
        };
        kwargs.set_item("genotype_by", gt_by)?;
        kwargs.set_item("samples_nested", model.samples_nested())?;
        kwargs.set_item("coords", model.coord_system().to_string())?;
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
        let fmt_reader = noodles::vcf::io::Reader::new(reader);
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
        let fmt_reader = noodles::vcf::io::Reader::new(reader);
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
                let fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
                let batch_reader = self
                    .scanner
                    .scan_virtual_ranges(fmt_reader, vpos_ranges, columns, batch_size, limit)
                    .map_err(to_py)?;
                Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
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
        region: String,
        index: Option<Py<PyAny>>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region =
            oxbow::Region::parse(&region, self.scanner.model().coord_system()).map_err(to_py)?;
        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
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
                let fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
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

/// A BCF file scanner.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the BCF file or a file-like object.
/// compressed : bool, optional [default: True]
///     Whether the source is BGZF-compressed. If None, it is assumed to be
///     compressed.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// info_fields : list[str], optional
///     Names of the INFO fields to project.
/// genotype_fields : list[str], optional
///     Names of the sample-specific genotype fields to project.
/// genotype_by : Literal["sample", "field"], optional [default: "sample"]
///     How to project the genotype fields. If "sample", the columns
///     correspond to the samples. If "field", the columns correspond to
///     the genotype fields.
/// samples : list[str] or None, optional [default: None]
///    Names of the samples to include in the genotype output. ``"*"`` for
///    all samples, a list to select specific samples, or ``None`` to omit
///    all sample genotype data.
/// samples_nested : bool, optional [default: False]
///   Whether to nest sample genotype data under a single ``"samples"`` struct
///   column.
#[pyclass(module = "oxbow.oxbow")]
pub struct PyBcfScanner {
    src: Py<PyAny>,
    reader: Reader,
    scanner: BcfScanner,
    compressed: bool,
}

#[pymethods]
impl PyBcfScanner {
    #[new]
    #[pyo3(signature = (src, compressed=true, fields=None, info_fields=None, genotype_fields=None, genotype_by=None, samples=None, samples_nested=false, coords=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python,
        src: Py<PyAny>,
        compressed: bool,
        fields: Option<Py<PyAny>>,
        info_fields: Option<Py<PyAny>>,
        genotype_fields: Option<Py<PyAny>>,
        genotype_by: Option<String>,
        samples: Option<Py<PyAny>>,
        samples_nested: bool,
        coords: Option<String>,
    ) -> PyResult<Self> {
        let coord_system = resolve_coord_system(coords)?;
        let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
        let mut fmt_reader = noodles::bcf::io::Reader::from(reader);
        let header = fmt_reader.read_header()?;
        let reader = fmt_reader.into_inner();
        let gt_by = resolve_genotype_by(genotype_by)?;
        let scanner = BcfScanner::new(
            header,
            resolve_fields(fields, py)?,
            resolve_fields(info_fields, py)?,
            resolve_fields(genotype_fields, py)?,
            gt_by,
            resolve_fields(samples, py)?,
            Some(samples_nested),
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
        if let Some(defs) = model.info_defs() {
            let names: Vec<&str> = defs.iter().map(|d| d.name.as_str()).collect();
            kwargs.set_item("info_fields", names)?;
        }
        if let Some(defs) = model.genotype_defs() {
            let names: Vec<&str> = defs.iter().map(|d| d.name.as_str()).collect();
            kwargs.set_item("genotype_fields", names)?;
        }
        if let Some(samples) = model.samples() {
            kwargs.set_item("samples", samples)?;
        }
        let gt_by = match model.genotype_by() {
            GenotypeBy::Sample => "sample",
            GenotypeBy::Field => "field",
        };
        kwargs.set_item("genotype_by", gt_by)?;
        kwargs.set_item("samples_nested", model.samples_nested())?;
        kwargs.set_item("coords", model.coord_system().to_string())?;
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
        let fmt_reader = noodles::bcf::io::Reader::from(reader);
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
        let fmt_reader = noodles::bcf::io::Reader::from(reader);
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
                let fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
                let batch_reader = self
                    .scanner
                    .scan_virtual_ranges(fmt_reader, vpos_ranges, columns, batch_size, limit)
                    .map_err(to_py)?;
                Ok(PyRecordBatchReader::new(err_on_unwind(batch_reader)))
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
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
        region: String,
        index: Option<Py<PyAny>>,
        columns: Option<Vec<String>>,
        batch_size: Option<usize>,
        limit: Option<usize>,
    ) -> PyResult<PyRecordBatchReader> {
        let region =
            oxbow::Region::parse(&region, self.scanner.model().coord_system()).map_err(to_py)?;

        match self.reader.clone() {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
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
                let fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
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
/// region : str, optional
///     Genomic region in the format "chr:start-end".
/// index : path or file-like, optional
///     The index file to use for querying the region.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// info_fields : list[str], optional
///     Names of the INFO fields to project.
/// genotype_fields : list[str], optional
///     Names of the sample-specific genotype fields to project.
/// genotype_by : Literal["sample", "field"], optional [default: "sample"]
///     How to project the genotype fields. If "sample", the columns
///     correspond to the samples. If "field", the columns correspond to
///     the genotype fields.
/// samples : list[str], optional
///     Names of the samples to include in the genotype fields.
/// samples_nested : bool, optional [default: False]
///     Whether to nest the sample-specific genotype fields under a "samples" struct column.
/// compressed : bool, optional [default: False]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, region=None, index=None, fields=None, info_fields=None, genotype_fields=None, genotype_by=None, samples=None, samples_nested=false, compressed=false))]
#[allow(clippy::too_many_arguments)]
pub fn read_vcf(
    py: Python,
    src: Py<PyAny>,
    region: Option<String>,
    index: Option<Py<PyAny>>,
    fields: Option<Py<PyAny>>,
    info_fields: Option<Py<PyAny>>,
    genotype_fields: Option<Py<PyAny>>,
    genotype_by: Option<String>,
    samples: Option<Py<PyAny>>,
    samples_nested: bool,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let mut fmt_reader = noodles::vcf::io::Reader::new(reader);
    let header = fmt_reader.read_header()?;
    let reader = fmt_reader.into_inner();
    let genotype_by = resolve_genotype_by(genotype_by)?;
    let scanner = VcfScanner::new(
        header,
        resolve_fields(fields, py)?,
        resolve_fields(info_fields, py)?,
        resolve_fields(genotype_fields, py)?,
        genotype_by,
        resolve_fields(samples, py)?,
        Some(samples_nested),
        CoordSystem::OneClosed,
    )
    .map_err(to_py)?;

    let ipc = if let Some(region) = region {
        let region = oxbow::Region::parse(&region, oxbow::CoordSystem::OneClosed).map_err(to_py)?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
                let index = resolve_index(py, &src, index)?;
                let batches = scanner
                    .scan_query(fmt_reader, region, index.into_boxed(), None, None, None)
                    .map_err(to_py)?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::vcf::io::Reader::new(bgzf_reader);
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
        let fmt_reader = noodles::vcf::io::Reader::new(reader);
        let batches = scanner.scan(fmt_reader, None, None, None).map_err(to_py)?;
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
/// region : str, optional
///     Genomic region in the format "chr:start-end".
/// index : path or file-like, optional
///     The index file to use for querying the region.
/// fields : list[str], optional
///     Names of the fixed fields to project.
/// info_fields : list[str], optional
///     Names of the INFO fields to project.
/// genotype_fields : list[str], optional
///     Names of the sample-specific genotype fields to project.
/// genotype_by : Literal["sample", "field"], optional [default: "sample"]
///     How to project the genotype fields. If "sample", the columns
///     correspond to the samples. If "field", the columns correspond to
///     the genotype fields.
/// samples : list[str], optional
///     Names of the samples to include in the genotype fields.
/// samples_nested : bool, optional [default: False]
///     Whether to nest the sample-specific genotype fields under a "samples" struct column.
/// compressed : bool, optional [default: True]
///     Whether the source is BGZF-compressed.
///
/// Returns
/// -------
/// bytes
///     Arrow IPC
#[pyfunction]
#[pyo3(signature = (src, region=None, index=None, fields=None, info_fields=None, genotype_fields=None, genotype_by=None, samples=None, samples_nested=false, compressed=true))]
#[allow(clippy::too_many_arguments)]
pub fn read_bcf(
    py: Python,
    src: Py<PyAny>,
    region: Option<String>,
    index: Option<Py<PyAny>>,
    fields: Option<Py<PyAny>>,
    info_fields: Option<Py<PyAny>>,
    genotype_fields: Option<Py<PyAny>>,
    genotype_by: Option<String>,
    samples: Option<Py<PyAny>>,
    samples_nested: bool,
    compressed: bool,
) -> PyResult<Vec<u8>> {
    let reader = pyobject_to_bufreader(py, src.clone_ref(py), compressed)?;
    let mut fmt_reader = noodles::bcf::io::Reader::from(reader);
    let header = fmt_reader.read_header()?;
    let reader = fmt_reader.into_inner();
    let genotype_by = resolve_genotype_by(genotype_by)?;
    let scanner = BcfScanner::new(
        header,
        resolve_fields(fields, py)?,
        resolve_fields(info_fields, py)?,
        resolve_fields(genotype_fields, py)?,
        genotype_by,
        resolve_fields(samples, py)?,
        Some(samples_nested),
        CoordSystem::OneClosed,
    )
    .map_err(to_py)?;

    let ipc = if let Some(region) = region {
        let region = oxbow::Region::parse(&region, oxbow::CoordSystem::OneClosed).map_err(to_py)?;

        match reader {
            Reader::BgzfFile(bgzf_reader) => {
                let fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
                let index = resolve_index(py, &src, index)?;
                let batches = scanner
                    .scan_query(fmt_reader, region, index.into_boxed(), None, None, None)
                    .map_err(to_py)?;
                batches_to_ipc(batches)
            }
            Reader::BgzfPyFileLike(bgzf_reader) => {
                let fmt_reader = noodles::bcf::io::Reader::from(bgzf_reader);
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
        let fmt_reader = noodles::bcf::io::Reader::from(reader);
        let batches = scanner.scan(fmt_reader, None, None, None).map_err(to_py)?;
        batches_to_ipc(batches)
    };

    ipc.map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

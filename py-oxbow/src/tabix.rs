use std::sync::Arc;

use arrow::array::{ArrayRef, Int32Array, RecordBatchIterator, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use noodles::core::Region;
use noodles::csi;
use noodles::csi::io::IndexedRecord;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyString;

use crate::filelike::PyFileLikeObject;
use crate::util::resolve_index;
use oxbow::util::batches_to_ipc;

/// Read records from a tabix-indexed file in the given genomic region.
///
/// Parameters
/// ----------
/// src : str or file-like
///     The path to the BGZF-compressed tabix-indexed file, or a seekable
///     file-like object.
/// region : str
///     Genomic region to query, e.g. ``"chr1:1000-2000"``.
/// index : str or file-like, optional
///     Path to the tabix index (.tbi or .csi), or a file-like object.
///     Required when ``src`` is a file-like object; inferred from ``src``
///     when ``src`` is a file path (appends ``.tbi``).
///
/// Returns
/// -------
/// bytes
///     Arrow IPC bytes with "chrom" (Utf8), "start" (Int32), "end" (Int32),
///     and "raw" (Utf8) columns.  Coordinates are 1-based, matching the
///     tabix index convention.
///
/// Raises
/// ------
/// ValueError
///     If ``src`` is a file-like and ``index`` is not provided, or if the
///     region chromosome is not present in the index
///     ("missing reference sequence name").
#[pyfunction]
#[pyo3(signature = (src, region, index=None))]
pub fn read_tabix(
    py: Python,
    src: Py<PyAny>,
    region: String,
    index: Option<Py<PyAny>>,
) -> PyResult<Vec<u8>> {
    let region = region
        .parse::<Region>()
        .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

    let csi_index = resolve_index(py, &src, index)?;

    let records: Vec<(String, i32, i32, String)> =
        if let Ok(py_string) = src.cast_bound::<PyString>(py) {
            let path = py_string.to_string();
            let file = std::fs::File::open(&path)?;
            collect_query(file, csi_index.into_boxed(), &region)
                .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?
        } else {
            let file_like = PyFileLikeObject::new(src, true, false, true)?;
            collect_query(file_like, csi_index.into_boxed(), &region)
                .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?
        };

    let schema = Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::Int32, false),
        Field::new("end", DataType::Int32, false),
        Field::new("raw", DataType::Utf8, false),
    ]));

    let chroms: Vec<&str> = records.iter().map(|(c, _, _, _)| c.as_str()).collect();
    let starts: Vec<i32> = records.iter().map(|(_, s, _, _)| *s).collect();
    let ends: Vec<i32> = records.iter().map(|(_, _, e, _)| *e).collect();
    let raws: Vec<&str> = records.iter().map(|(_, _, _, r)| r.as_str()).collect();

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(chroms)) as ArrayRef,
            Arc::new(Int32Array::from(starts)) as ArrayRef,
            Arc::new(Int32Array::from(ends)) as ArrayRef,
            Arc::new(StringArray::from(raws)) as ArrayRef,
        ],
    )
    .map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))?;

    let reader = RecordBatchIterator::new(vec![Ok(batch)], schema);
    batches_to_ipc(reader).map_err(|e| PyErr::new::<PyValueError, _>(e.to_string()))
}

fn collect_query<R: std::io::Read + std::io::Seek>(
    reader: R,
    index: Box<dyn csi::BinningIndex>,
    region: &Region,
) -> std::io::Result<Vec<(String, i32, i32, String)>> {
    let mut indexed_reader = csi::io::IndexedReader::new(reader, index);
    let query = indexed_reader.query(region)?;
    query
        .map(|result| {
            let record = result?;
            let chrom = record.indexed_reference_sequence_name().to_string();
            let start = record.indexed_start_position().get() as i32;
            let end = record.indexed_end_position().get() as i32;
            let raw = record.as_ref().to_string();
            Ok((chrom, start, end, raw))
        })
        .collect()
}

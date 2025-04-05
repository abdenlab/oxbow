use std::io::BufReader;
use std::io::{BufRead, Read, Seek};

use pyo3::prelude::*;
use pyo3::types::PyString;

use noodles::bgzf::gzi::Index as GzIndex;
use noodles::bgzf::io::Seek as BgzfSeek;
use noodles::bgzf::VirtualPosition;
use noodles::fasta::fai::Index as FaIndex;

use crate::filelike::PyFileLikeObject;
use oxbow::util::index::IndexType;
use oxbow::util::index::{index_from_reader, index_from_source_path};

pub const BUFFER_SIZE_BYTES: usize = const { 1024 * 1024 };

/// A cloneable buffered reader that can read from files or file-like objects.
pub enum Reader {
    File(BufReader<std::fs::File>),
    PyFileLike(BufReader<PyFileLikeObject>),
    BgzfFile(noodles::bgzf::Reader<BufReader<std::fs::File>>),
    BgzfPyFileLike(noodles::bgzf::Reader<BufReader<PyFileLikeObject>>),
}

impl Read for Reader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::File(reader) => reader.read(buf),
            Self::PyFileLike(reader) => reader.read(buf),
            Self::BgzfFile(reader) => reader.read(buf),
            Self::BgzfPyFileLike(reader) => reader.read(buf),
        }
    }
}

impl BufRead for Reader {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self {
            Self::File(reader) => reader.fill_buf(),
            Self::PyFileLike(reader) => reader.fill_buf(),
            Self::BgzfFile(reader) => reader.fill_buf(),
            Self::BgzfPyFileLike(reader) => reader.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            Self::File(reader) => reader.consume(amt),
            Self::PyFileLike(reader) => reader.consume(amt),
            Self::BgzfFile(reader) => reader.consume(amt),
            Self::BgzfPyFileLike(reader) => reader.consume(amt),
        }
    }
}

impl Seek for Reader {
    fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        match self {
            Self::File(reader) => reader.seek(pos),
            Self::PyFileLike(reader) => reader.seek(pos),
            Self::BgzfFile(_) => panic!("Seek not supported for BGZF files."),
            Self::BgzfPyFileLike(_) => panic!("Seek not supported for BGZF files."),
        }
    }
}

impl BgzfSeek for Reader {
    fn seek_to_virtual_position(
        &mut self,
        pos: VirtualPosition,
    ) -> std::io::Result<VirtualPosition> {
        match self {
            Self::File(_) => panic!("Seek not supported for non-BGZF files."),
            Self::PyFileLike(_) => panic!("Seek not supported for non-BGZF files."),
            Self::BgzfFile(reader) => reader.seek(pos),
            Self::BgzfPyFileLike(reader) => reader.seek(pos),
        }
    }

    fn seek_with_index(&mut self, index: &GzIndex, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        match self {
            Self::File(_) => panic!("Seek not supported for non-BGZF files."),
            Self::PyFileLike(_) => panic!("Seek not supported for non-BGZF files."),
            Self::BgzfFile(reader) => reader.seek_with_index(index, pos),
            Self::BgzfPyFileLike(reader) => reader.seek_with_index(index, pos),
        }
    }
}

impl Reader {
    #[allow(clippy::seek_from_current)]
    pub fn clone(&mut self) -> Self {
        match self {
            Self::File(source) => {
                let _ = source.seek(std::io::SeekFrom::Current(0)).unwrap(); // discards buffer
                let file = source.get_ref().try_clone().unwrap();
                let reader = BufReader::with_capacity(BUFFER_SIZE_BYTES, file);
                Self::File(reader)
            }
            Self::PyFileLike(source) => {
                let _ = source.seek(std::io::SeekFrom::Current(0)).unwrap(); // discards buffer
                let file_like = source.get_ref().clone();
                let reader = BufReader::with_capacity(BUFFER_SIZE_BYTES, file_like);
                Self::PyFileLike(reader)
            }
            Self::BgzfFile(source) => {
                let vpos = source.virtual_position();
                let file = source.get_ref().get_ref().try_clone().unwrap();
                let mut reader =
                    noodles::bgzf::Reader::new(BufReader::with_capacity(BUFFER_SIZE_BYTES, file));
                reader.seek_to_virtual_position(vpos).unwrap();
                Self::BgzfFile(reader)
            }
            Self::BgzfPyFileLike(source) => {
                let vpos = source.virtual_position();
                let file_like = source.get_ref().get_ref().clone();
                let mut reader = noodles::bgzf::Reader::new(BufReader::with_capacity(
                    BUFFER_SIZE_BYTES,
                    file_like,
                ));
                reader.seek_to_virtual_position(vpos).unwrap();
                Self::BgzfPyFileLike(reader)
            }
        }
    }
}

pub fn pyobject_to_bufreader(
    py: Python,
    obj: PyObject,
    compressed: bool,
) -> std::io::Result<Reader> {
    if let Ok(string_ref) = obj.downcast_bound::<PyString>(py) {
        let path = string_ref.to_string_lossy().into_owned();
        let file = std::fs::File::open(path)?;
        let reader = BufReader::with_capacity(BUFFER_SIZE_BYTES, file);
        if compressed {
            let reader = noodles::bgzf::Reader::new(reader);
            Ok(Reader::BgzfFile(reader))
        } else {
            Ok(Reader::File(reader))
        }
    } else {
        let file_like = PyFileLikeObject::new(obj, true, false, true)?;
        let reader = BufReader::with_capacity(BUFFER_SIZE_BYTES, file_like);
        if compressed {
            let reader = noodles::bgzf::Reader::new(reader);
            Ok(Reader::BgzfPyFileLike(reader))
        } else {
            Ok(Reader::PyFileLike(reader))
        }
    }
}

#[allow(dead_code)]
pub fn resolve_index(
    py: Python,
    source: PyObject,
    index: Option<PyObject>,
) -> std::io::Result<IndexType> {
    match index {
        // Index file not provided
        None => match source.downcast_bound::<PyString>(py) {
            // If source is a path, try to find companion index file path
            Ok(py_string) => {
                let source_path = py_string.to_string();
                index_from_source_path(&source_path)
            }
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "No index file found.",
            )),
        },
        // Index file explicitly provided
        Some(index) => {
            if let Ok(py_string) = index.downcast_bound::<PyString>(py) {
                let path = py_string.to_string();
                let index_file = std::fs::File::open(&path)?;
                index_from_reader(index_file)
            } else {
                let index_pyfile = PyFileLikeObject::new(index, true, false, true)?;
                index_from_reader(index_pyfile)
            }
        }
    }
}

pub fn resolve_faidx(
    py: Python,
    source: PyObject,
    faidx: Option<PyObject>,
) -> std::io::Result<FaIndex> {
    match faidx {
        // Index file not provided
        None => match source.downcast_bound::<PyString>(py) {
            // If source is a path, try to find companion index file path
            Ok(py_string) => {
                let source_path = py_string.to_string();
                let index_path = format!("{}.fai", source_path);
                let reader = BufReader::new(std::fs::File::open(&index_path)?);
                let mut index_reader = noodles::fasta::fai::io::Reader::new(reader);
                index_reader.read_index()
            }
            _ => Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "No FAI index file found.",
            )),
        },
        // Index file explicitly provided
        Some(faidx) => {
            if let Ok(py_string) = faidx.downcast_bound::<PyString>(py) {
                let index_path = py_string.to_string();
                let reader = BufReader::new(std::fs::File::open(&index_path)?);
                let mut index_reader = noodles::fasta::fai::io::Reader::new(reader);
                index_reader.read_index()
            } else {
                let pyfile = PyFileLikeObject::new(faidx, true, false, true)?;
                let reader = BufReader::new(pyfile);
                let mut index_reader = noodles::fasta::fai::io::Reader::new(reader);
                index_reader.read_index()
            }
        }
    }
}

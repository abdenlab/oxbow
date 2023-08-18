use std::io::{self, Read, Seek, SeekFrom, Write};

use pyo3::{
    exceptions::PyTypeError, types::PyBytes, IntoPy, PyErr, PyObject, PyResult, Python, ToPyObject,
};

/// Represents a file-like object in python. This simply wraps the Rust io
/// traits, calling into python io methods.
#[derive(Clone)]
pub struct PyFileLikeObject {
    inner: PyObject,
}

impl PyFileLikeObject {
    /// Creates a `PyFileLikeObject`, validating that it conforms to some
    /// required set of methods (one or more of `read`, `write`, or `seek`).
    /// Will return a `TypeError` if object does not have the required methods.
    pub fn new(object: PyObject, read: bool, write: bool, seek: bool) -> PyResult<Self> {
        Python::with_gil(|py| {
            if read && object.getattr(py, "read").is_err() {
                return Err(PyErr::new::<PyTypeError, _>(
                    "Object does not have a .read() method.",
                ));
            }

            if seek && object.getattr(py, "seek").is_err() {
                return Err(PyErr::new::<PyTypeError, _>(
                    "Object does not have a .seek() method.",
                ));
            }

            if write && object.getattr(py, "write").is_err() {
                return Err(PyErr::new::<PyTypeError, _>(
                    "Object does not have a .write() method.",
                ));
            }

            Ok(PyFileLikeObject { inner: object })
        })
    }
}

/// Extracts a string repr from, and returns an IO error to send back to rust.
fn to_io_error(py: Python<'_>, e: PyErr) -> io::Error {
    let pyobj: PyObject = e.into_py(py);

    match pyobj.call_method(py, "__str__", (), None) {
        Ok(repr) => match repr.extract::<String>(py) {
            Ok(s) => io::Error::new(io::ErrorKind::Other, s),
            Err(_e) => io::Error::new(io::ErrorKind::Other, "An unknown error has occurred"),
        },
        Err(_) => io::Error::new(io::ErrorKind::Other, "An unknown error has occurred"),
    }
}

impl Read for PyFileLikeObject {
    fn read(&mut self, mut buf: &mut [u8]) -> Result<usize, io::Error> {
        Python::with_gil(|py| {
            let res = self
                .inner
                .call_method(py, "read", (buf.len(),), None)
                .map_err(|e| to_io_error(py, e))?;
            let pybytes: &PyBytes = res.downcast(py).map_err(|e| to_io_error(py, e.into()))?;
            let bytes = pybytes.as_bytes();
            buf.write_all(bytes)?;
            Ok(bytes.len())
        })
    }
}

impl Write for PyFileLikeObject {
    fn write(&mut self, buf: &[u8]) -> Result<usize, io::Error> {
        Python::with_gil(|py| {
            let arg = PyBytes::new(py, buf).to_object(py);

            let number_bytes_written = self
                .inner
                .call_method(py, "write", (arg,), None)
                .map_err(|e| to_io_error(py, e))?;

            if number_bytes_written.is_none(py) {
                return Err(io::Error::new(io::ErrorKind::Other, "no bytes written"));
            }

            number_bytes_written
                .extract(py)
                .map_err(|e| to_io_error(py, e))
        })
    }

    fn flush(&mut self) -> Result<(), io::Error> {
        Python::with_gil(|py| {
            self.inner
                .call_method(py, "flush", (), None)
                .map_err(|e| to_io_error(py, e))?;

            Ok(())
        })
    }
}

impl Seek for PyFileLikeObject {
    fn seek(&mut self, pos: SeekFrom) -> Result<u64, io::Error> {
        Python::with_gil(|py| {
            let (whence, offset) = match pos {
                SeekFrom::Start(i) => (0, i as i64),
                SeekFrom::Current(i) => (1, i),
                SeekFrom::End(i) => (2, i),
            };

            let new_position = self
                .inner
                .call_method(py, "seek", (offset, whence), None)
                .map_err(|e| to_io_error(py, e))?;

            new_position.extract(py).map_err(|e| to_io_error(py, e))
        })
    }
}

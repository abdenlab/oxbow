use std::backtrace::Backtrace;
use std::fmt;

use arrow::error::ArrowError;

/// A captured backtrace, wrapped to avoid triggering thiserror's unstable
/// `Error::provide()` codegen (which requires `error_generic_member_access`).
pub struct CapturedBacktrace(Backtrace);

impl CapturedBacktrace {
    fn capture() -> Self {
        Self(Backtrace::capture())
    }

    /// Returns the backtrace status.
    pub fn status(&self) -> std::backtrace::BacktraceStatus {
        self.0.status()
    }
}

impl fmt::Debug for CapturedBacktrace {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

impl fmt::Display for CapturedBacktrace {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// The unified error type for the oxbow crate.
#[derive(Debug, thiserror::Error)]
pub enum OxbowError {
    /// Invalid user-provided configuration: unknown column names, invalid field
    /// definitions, invalid tag type codes, bad schema parameters, etc.
    #[error("{message}")]
    InvalidInput {
        message: String,
        trace: CapturedBacktrace,
    },

    /// Data does not conform to the expected format during parsing: type
    /// mismatches in tag values, INFO fields, genotype fields, BED fields,
    /// malformed records, etc.
    #[error("{message}")]
    InvalidData {
        message: String,
        trace: CapturedBacktrace,
    },

    /// A required resource was not found: missing reference sequences, missing
    /// index headers, missing samples.
    #[error("{message}")]
    NotFound {
        message: String,
        trace: CapturedBacktrace,
    },

    /// An I/O error from the underlying reader/writer.
    #[error(transparent)]
    Io(std::io::Error),

    /// An Arrow error from batch construction or schema operations.
    #[error(transparent)]
    Arrow(ArrowError),

    /// An error from an upstream library (noodles, bigtools, etc.) that
    /// doesn't fit the other categories.
    #[error(transparent)]
    External(Box<dyn std::error::Error + Send + Sync>),
}

/// A specialized `Result` type for oxbow operations.
pub type Result<T> = std::result::Result<T, OxbowError>;

impl OxbowError {
    /// Creates an `InvalidInput` error with the given message.
    pub fn invalid_input(message: impl Into<String>) -> Self {
        Self::InvalidInput {
            message: message.into(),
            trace: CapturedBacktrace::capture(),
        }
    }

    /// Creates an `InvalidData` error with the given message.
    pub fn invalid_data(message: impl Into<String>) -> Self {
        Self::InvalidData {
            message: message.into(),
            trace: CapturedBacktrace::capture(),
        }
    }

    /// Creates a `NotFound` error with the given message.
    pub fn not_found(message: impl Into<String>) -> Self {
        Self::NotFound {
            message: message.into(),
            trace: CapturedBacktrace::capture(),
        }
    }

    /// Returns the captured backtrace, if any.
    pub fn backtrace(&self) -> Option<&CapturedBacktrace> {
        match self {
            Self::InvalidInput { trace, .. }
            | Self::InvalidData { trace, .. }
            | Self::NotFound { trace, .. } => Some(trace),
            Self::Io(_) | Self::Arrow(_) | Self::External(_) => None,
        }
    }

    /// Formats the error message with the Rust backtrace appended, if one was
    /// captured. Useful for forwarding to Python exception messages.
    pub fn display_with_backtrace(&self) -> String {
        if let Some(bt) = self.backtrace() {
            if bt.status() == std::backtrace::BacktraceStatus::Captured {
                return format!("{self}\n\nRust backtrace:\n{bt}");
            }
        }
        self.to_string()
    }
}

impl From<std::io::Error> for OxbowError {
    fn from(err: std::io::Error) -> Self {
        Self::Io(err)
    }
}

impl From<ArrowError> for OxbowError {
    fn from(err: ArrowError) -> Self {
        Self::Arrow(err)
    }
}

impl From<OxbowError> for ArrowError {
    fn from(err: OxbowError) -> Self {
        // Bake the backtrace into the error message because only the string
        // survives the Arrow C Stream Interface boundary to Python/pyarrow.
        let msg = err.display_with_backtrace();
        ArrowError::ExternalError(Box::new(std::io::Error::other(msg)))
    }
}

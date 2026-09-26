use thiserror::Error;

#[derive(Error, Debug)]
pub enum Error {
    #[error("invalid zero-based scan index: {0}")]
    InvalidScanIndex(u32),

    #[error("invalid scan number: {0}")]
    InvalidScanNumber(u32),

    #[error("corrupt data: {0}")]
    CorruptData(&'static str),

    #[error("unsupported layout: {0}")]
    UnsupportedLayout(&'static str),

    #[error("{context}: {source}")]
    Context {
        context: ErrorContext,
        #[source]
        source: Box<Error>,
    },

    #[error("invalid scan event: {0}")]
    InvalidScanEvent(String),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("invalid magic number: expected 0xa101, got 0x{0:04x}")]
    BadMagic(u16),

    #[error("invalid signature: expected \"Finnigan\", got {0:?}")]
    BadSignature(String),

    #[error("unsupported version: {0}")]
    UnsupportedVersion(u32),

    #[error("unexpected end of file at offset {offset:#x} (needed {needed} bytes)")]
    UnexpectedEof { offset: u64, needed: usize },

    #[error("invalid UTF-16 data at offset {0:#x}")]
    InvalidUtf16(u64),

    #[error("address out of range: {0:#x}")]
    AddressOutOfRange(u64),

    #[error("invalid generic data type code: {0}")]
    InvalidTypeCode(u32),

    #[error("operation not supported for this scan-data format: {0}")]
    UnsupportedOperation(&'static str),

    #[error(
        "refusing to allocate for {requested} bytes at offset {offset:#x}: \
         only {available} bytes remain in the input"
    )]
    AllocationTooLarge {
        offset: u64,
        requested: u64,
        available: u64,
    },
}

pub type Result<T> = std::result::Result<T, Error>;

/// Stable categories for logging and recovery without parsing error strings.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ErrorKind {
    Io,
    Corrupt,
    Unsupported,
    InvalidInput,
    /// Reserved for an explicit configured resource budget. Declared data
    /// extending beyond the input is corruption, not a resource budget failure.
    ResourceLimit,
}

/// The operation that failed. Nested contexts retain lower-level operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Operation {
    Open,
    ReadScan,
    DecodeSpectrum,
    ValidateExport,
    WriteMzml,
    WriteIndexedMzml,
}

/// Machine-readable location of a failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ErrorContext {
    pub operation: Operation,
    /// Input path for errors from `RawFileReader::open_path`.
    pub path: Option<std::path::PathBuf>,
    /// Zero-based offset from the first scan in the run.
    pub scan_index: Option<u32>,
    pub scan_number: Option<u32>,
    /// Absolute packet start in the file, when known. An inner decoder offset
    /// may instead be relative to this packet. This is not an offending-byte claim.
    pub byte_offset: Option<u64>,
}

impl ErrorContext {
    pub fn new(operation: Operation) -> Self {
        Self {
            operation,
            path: None,
            scan_index: None,
            scan_number: None,
            byte_offset: None,
        }
    }
}

impl std::fmt::Display for ErrorContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.operation)?;
        if let Some(path) = &self.path {
            write!(f, " for {}", path.display())?;
        }
        if let Some(index) = self.scan_index {
            write!(f, " at scan index {index}")?;
        }
        if let Some(number) = self.scan_number {
            write!(f, " (scan {number})")?;
        }
        if let Some(offset) = self.byte_offset {
            write!(f, " at packet offset {offset:#x}")?;
        }
        Ok(())
    }
}

impl Error {
    pub fn with_context(self, context: ErrorContext) -> Self {
        Self::Context {
            context,
            source: Box::new(self),
        }
    }

    pub fn context(&self) -> Option<&ErrorContext> {
        match self {
            Self::Context { context, .. } => Some(context),
            _ => None,
        }
    }

    /// Return the underlying reader error, retaining I/O sources on that error.
    pub fn root_cause(&self) -> &Self {
        match self {
            Self::Context { source, .. } => source.root_cause(),
            _ => self,
        }
    }

    pub fn kind(&self) -> ErrorKind {
        match self.root_cause() {
            Self::Io(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => {
                ErrorKind::Corrupt
            }
            Self::Io(_) => ErrorKind::Io,
            Self::InvalidScanIndex(_) | Self::InvalidScanNumber(_) => ErrorKind::InvalidInput,
            Self::UnsupportedVersion(_)
            | Self::UnsupportedOperation(_)
            | Self::UnsupportedLayout(_) => ErrorKind::Unsupported,
            Self::AllocationTooLarge { .. } => ErrorKind::Corrupt,
            _ => ErrorKind::Corrupt,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error as _;

    #[test]
    fn context_retains_typed_cause_and_io_source() {
        let error = Error::Io(std::io::Error::from(std::io::ErrorKind::PermissionDenied))
            .with_context(ErrorContext::new(Operation::ReadScan))
            .with_context(ErrorContext::new(Operation::DecodeSpectrum));
        assert_eq!(error.kind(), ErrorKind::Io);
        assert_eq!(
            error.context().unwrap().operation,
            Operation::DecodeSpectrum
        );
        assert!(matches!(error.root_cause(), Error::Io(_)));
        assert!(error.source().unwrap().source().unwrap().source().is_some());
    }

    #[test]
    fn categories_separate_input_corruption_and_unsupported_formats() {
        assert_eq!(
            Error::AllocationTooLarge {
                offset: 0,
                requested: 100,
                available: 20
            }
            .kind(),
            ErrorKind::Corrupt
        );
        assert_eq!(Error::InvalidScanIndex(2).kind(), ErrorKind::InvalidInput);
        assert_eq!(
            Error::UnsupportedVersion(100).kind(),
            ErrorKind::Unsupported
        );
        assert_eq!(
            Error::CorruptData("nonfinite array").kind(),
            ErrorKind::Corrupt
        );
        assert_eq!(
            Error::Io(std::io::Error::from(std::io::ErrorKind::UnexpectedEof)).kind(),
            ErrorKind::Corrupt
        );
    }
}

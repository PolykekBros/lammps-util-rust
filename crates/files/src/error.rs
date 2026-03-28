use std::{error, fmt, io, num};

use crate::snapshot::{HEADER_ATOMS, HEADER_NUM_OF_ATOMS, HEADER_SYM_BOX, HEADER_TIMESTEP};

pub struct Error {
    err: Box<ErrorImpl>,
}

struct ErrorImpl {
    kind: ErrorKind,
    content: String,
    line: usize,
}

#[derive(Debug)]
pub enum ErrorKind {
    ExpectedTimestepHeader,
    MissingTimestep,
    InvalidTimestep(num::ParseIntError),
    ExpectedAtomCountHeader,
    MissingAtomCount,
    InvalidAtomCount(num::ParseIntError),
    ExpectedSymboxHeader,
    MissingSymboxField,
    InvalidSymboxField(num::ParseFloatError),
    ExpectedAtomsHeader,
    MissingAtomKeys,
    DuplicateAtomKeys(String),
    DuplicateTimestep(u64),
    MissingAtomRowField,
    InvalidAtomRowField(num::ParseFloatError),
    Io(io::Error),
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self {
            ErrorKind::ExpectedTimestepHeader => write!(f, "expected {HEADER_TIMESTEP}"),
            ErrorKind::MissingTimestep => write!(f, "missing timestep value"),
            ErrorKind::InvalidTimestep(e) => write!(f, "invalid timestep: {e}"),
            ErrorKind::ExpectedAtomCountHeader => write!(f, "expected {HEADER_NUM_OF_ATOMS}"),
            ErrorKind::MissingAtomCount => write!(f, "missing atom count"),
            ErrorKind::InvalidAtomCount(e) => write!(f, "invalid atom count: {e}"),
            ErrorKind::ExpectedSymboxHeader => write!(f, "expected {HEADER_SYM_BOX}"),
            ErrorKind::MissingSymboxField => write!(f, "missing box dimension field"),
            ErrorKind::InvalidSymboxField(e) => write!(f, "invalid box dimension: {e}"),
            ErrorKind::ExpectedAtomsHeader => write!(f, "expected {HEADER_ATOMS}"),
            ErrorKind::MissingAtomKeys => write!(f, "missing atom attribute keys"),
            ErrorKind::DuplicateAtomKeys(key) => write!(f, "duplicate atom key: '{key}'"),
            ErrorKind::DuplicateTimestep(timestep) => write!(f, "duplicate timestep: '{timestep}'"),
            ErrorKind::MissingAtomRowField => write!(f, "missing field in atom row"),
            ErrorKind::InvalidAtomRowField(e) => write!(f, "invalid atom field: {e}"),
            ErrorKind::Io(e) => write!(f, "IO error: {e}"),
        }
    }
}

impl error::Error for ErrorKind {}

impl Error {
    pub(crate) fn new(kind: ErrorKind, content: String, line: usize) -> Self {
        Self {
            err: Box::new(ErrorImpl {
                kind,
                content,
                line,
            }),
        }
    }

    #[must_use]
    pub fn kind(&self) -> &ErrorKind {
        &self.err.kind
    }

    #[must_use]
    pub fn content(&self) -> &str {
        &self.err.content
    }

    #[must_use]
    pub fn line(&self) -> usize {
        self.err.line
    }
}

impl fmt::Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Error")
            .field("kind", &self.err.kind)
            .field("content", &self.err.content)
            .field("line", &self.err.line)
            .finish()
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "error at line {}: ", self.err.line)?;
        write!(f, "{}", self.kind())?;
        if !self.err.content.is_empty() {
            write!(f, " (found: {:?})", self.err.content)?;
        }
        Ok(())
    }
}

impl error::Error for Error {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match &self.err.kind {
            ErrorKind::Io(e) => Some(e),
            ErrorKind::InvalidTimestep(e) | ErrorKind::InvalidAtomCount(e) => Some(e),
            ErrorKind::InvalidSymboxField(e) | ErrorKind::InvalidAtomRowField(e) => Some(e),
            _ => None,
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

use std::{error, fmt, io, num};

pub struct Error {
    err: Box<ErrorImpl>,
}

struct ErrorImpl {
    kind: ErrorKind,
    content: String,
    line: usize,
    column: usize,
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
    MissingAtomRowField,
    InvalidAtomRowField(num::ParseFloatError),
    Io(io::Error),
}

impl Error {
    /// Internal helper to create a new error
    pub(crate) fn new(kind: ErrorKind, content: String, line: usize, column: usize) -> Self {
        Self {
            err: Box::new(ErrorImpl {
                kind,
                content,
                line,
                column,
            }),
        }
    }

    pub fn kind(&self) -> &ErrorKind {
        &self.err.kind
    }

    pub fn content(&self) -> &str {
        &self.err.content
    }

    pub fn line(&self) -> usize {
        self.err.line
    }

    pub fn column(&self) -> usize {
        self.err.column
    }
}

impl fmt::Debug for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // For debugging, we show the full internal state
        f.debug_struct("Error")
            .field("kind", &self.err.kind)
            .field("content", &self.err.content)
            .field("line", &self.err.line)
            .field("column", &self.err.column)
            .finish()
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "error at {}:{}: ", self.err.line, self.err.column)?;

        match &self.err.kind {
            ErrorKind::ExpectedTimestepHeader => write!(f, "expected ITEM: TIMESTEP"),
            ErrorKind::MissingTimestep => write!(f, "missing timestep value"),
            ErrorKind::InvalidTimestep(e) => write!(f, "invalid timestep: {e}"),
            ErrorKind::ExpectedAtomCountHeader => write!(f, "expected ITEM: NUMBER OF ATOMS"),
            ErrorKind::MissingAtomCount => write!(f, "missing atom count"),
            ErrorKind::InvalidAtomCount(e) => write!(f, "invalid atom count: {e}"),
            ErrorKind::ExpectedSymboxHeader => write!(f, "expected ITEM: BOX BOUNDS"),
            ErrorKind::MissingSymboxField => write!(f, "missing box dimension field"),
            ErrorKind::InvalidSymboxField(e) => write!(f, "invalid box dimension: {e}"),
            ErrorKind::ExpectedAtomsHeader => write!(f, "expected ITEM: ATOMS"),
            ErrorKind::MissingAtomKeys => write!(f, "missing atom attribute keys"),
            ErrorKind::DuplicateAtomKeys(key) => write!(f, "duplicate atom key: '{key}'"),
            ErrorKind::MissingAtomRowField => write!(f, "missing field in atom row"),
            ErrorKind::InvalidAtomRowField(e) => write!(f, "invalid atom field: {e}"),
            ErrorKind::Io(e) => return write!(f, "IO error: {e}"), // Exit early, no content to show
        }?;

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
            ErrorKind::InvalidTimestep(e) => Some(e),
            ErrorKind::InvalidAtomCount(e) => Some(e),
            ErrorKind::InvalidSymboxField(e) => Some(e),
            ErrorKind::InvalidAtomRowField(e) => Some(e),
            _ => None,
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Self {
        Self::new(ErrorKind::Io(err), String::new(), 0, 0)
    }
}

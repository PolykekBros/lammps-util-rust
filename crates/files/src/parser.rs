use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use crate::error::{Error, ErrorKind, Result};

#[derive(Debug)]
pub struct Parser<B> {
    br: B,
    buffer: String,
    is_err: bool,
    current: usize,
}

impl Parser<BufReader<File>> {
    pub fn open<P: AsRef<Path>>(p: P) -> Result<Self> {
        let f = File::open(p).map_err(|e| Error::new(ErrorKind::Io(e), String::new(), 0))?;
        let br = BufReader::new(f);
        Ok(Self::new(br))
    }
}

impl<B: BufRead> Parser<B> {
    pub fn new(br: B) -> Self {
        Self {
            br,
            buffer: String::with_capacity(1024),
            current: 0,
            is_err: false,
        }
    }

    pub fn next(&mut self) -> Option<Result<(&str, usize)>> {
        if self.is_err {
            return None;
        }
        self.buffer.clear();
        let line = match self.br.read_line(&mut self.buffer) {
            Ok(0) => return None,
            Ok(_) => Ok((self.buffer.trim(), self.current + 1)),
            Err(e) => {
                self.is_err = true;
                Err(Error::new(
                    ErrorKind::Io(e),
                    String::new(),
                    self.current + 1,
                ))
            }
        };
        self.current += 1;
        Some(line)
    }

    pub fn parse_line<F, T>(&mut self, f: F) -> Result<T>
    where
        F: FnOnce(&str) -> std::result::Result<T, ErrorKind>,
    {
        let (line, current) = self.next().transpose()?.unwrap_or_default();
        f(line).map_err(|kind| Error::new(kind, line.to_string(), current))
    }

    pub fn skip(&mut self, n: usize) -> Result<()> {
        for _ in 0..n {
            if self.next().transpose()?.is_none() {
                break;
            }
        }
        Ok(())
    }

    pub fn get_current(&self) -> usize {
        self.current
    }
}

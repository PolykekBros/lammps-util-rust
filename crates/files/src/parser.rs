use std::{
    fs::File,
    path::Path,
    io::{Lines, BufRead, BufReader},
};

use crate::error::{Result, Error};

#[derive(Debug)]
pub struct Parser<B> {
    lines: Lines<B>,
    is_err: bool,
    current: usize,
}

impl Parser<BufReader<File>> {
    pub fn open<P: AsRef<Path>>(p: P) -> Result<Self> {
        let f = File::open(p)?;
        let br = BufReader::new(f); 
        Ok(Self::new(br))
    }
}

impl<B: BufRead> Parser<B> {
    pub fn new(reader: B) -> Self {
        Self {
            lines: reader.lines(),
            current: 0,
            is_err: false,
        }
    }

    pub fn next(&mut self) -> Option<Result<String>> {
        let next = self.lines.next()?
            .map_err(Error::from);
        if !self.is_err {
            self.current += 1;        
        }
        if next.is_err() {
            self.is_err = true;
        }
        Some(next)
    }

    pub fn current(&self) -> usize {
        self.current
    }
}

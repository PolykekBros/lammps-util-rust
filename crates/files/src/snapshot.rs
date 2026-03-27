use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

// Assuming geomutil_util is available as per previous implementation
use geomutil_util::{BoundingBox3, Point3};

use crate::{
    error::{Error, ErrorKind, Result},
    parser::Parser,
};

pub const HEADER_TIMESTEP: &str = "ITEM: TIMESTEP";
pub const HEADER_NUM_OF_ATOMS: &str = "ITEM: NUMBER OF ATOMS";
pub const HEADER_SYM_BOX: &str = "ITEM: BOX BOUNDS";
pub const HEADER_ATOMS: &str = "ITEM: ATOMS";

#[derive(Debug, Clone)]
pub struct SymBox {
    pub boundaries: String,
    pub bbox: BoundingBox3<f64>,
}

#[derive(Debug, Clone)]
pub struct Snapshot {
    pub timestep: u64,
    pub atoms_count: usize,
    pub sym_box: SymBox,
    pub keys: Vec<String>,
    pub atoms: Vec<Vec<f64>>,
}

impl Snapshot {
    pub fn parse<B: BufRead>(parser: &mut Parser<B>) -> Result<Self> {
        let timestep = parse_timestep(parser)?;
        let atoms_count = parse_atom_count(parser)?;
        let sym_box = parse_sym_box(parser)?;
        let (keys, atoms) = parse_atoms(parser, atoms_count)?;

        Ok(Self {
            timestep,
            atoms_count,
            sym_box,
            keys,
            atoms,
        })
    }

    pub fn get_atom_attr(&self, atom_idx: usize, key: &str) -> Option<f64> {
        let key_idx = self.keys.iter().position(|k| k == key)?;
        self.atoms.get(atom_idx)?.get(key_idx).copied()
    }
}

pub struct Snapshots<B> {
    parser: Parser<B>,
}

impl Snapshots<BufReader<File>> {
    pub fn open<P: AsRef<Path>>(p: P) -> Result<Self> {
        let parser = Parser::open(p)?;
        Ok(Self::new(parser))
    }
}

impl<B: BufRead> Snapshots<B> {
    pub fn new(parser: Parser<B>) -> Self {
        Self { parser }
    }
}

impl<B: BufRead> Iterator for Snapshots<B> {
    type Item = Result<Snapshot>;

    fn next(&mut self) -> Option<Self::Item> {
        match Snapshot::parse(&mut self.parser) {
            Ok(s) => Some(Ok(s)),
            Err(e) => {
                // If we are at the very end of the file and ExpectedTimestepHeader happened with empty content,
                // it means parser.next() returned None, which indicates EOF before any new snapshot.
                match e.kind() {
                    ErrorKind::ExpectedTimestepHeader if e.content().is_empty() => None,
                    _ => Some(Err(e)),
                }
            }
        }
    }
}

fn parse_timestep<B: BufRead>(parser: &mut Parser<B>) -> Result<u64> {
    let line = parser.next().ok_or_else(|| {
        Error::new(
            ErrorKind::ExpectedTimestepHeader,
            String::new(),
            parser.current(),
            1,
        )
    })??;
    if line.trim() != HEADER_TIMESTEP {
        return Err(Error::new(
            ErrorKind::ExpectedTimestepHeader,
            line,
            parser.current(),
            1,
        ));
    }
    let timestep = parser.next().ok_or_else(|| {
        Error::new(
            ErrorKind::MissingTimestep,
            String::new(),
            parser.current(),
            1,
        )
    })??;
    timestep.parse::<u64>().map_err(|e| {
        Error::new(
            ErrorKind::InvalidTimestep(e),
            val_line,
            parser.current(),
            1,
        )
    })
}

fn parse_atom_count<B: BufRead>(parser: &mut Parser<B>) -> Result<usize> {
    let line = parser.next().ok_or_else(|| {
        Error::new(
            ErrorKind::ExpectedAtomCountHeader,
            String::new(),
            parser.current(),
            1,
        )
    })??;

    if line.trim() != HEADER_NUM_OF_ATOMS {
        let col = line.find(line.trim()).unwrap_or(0) + 1;
        return Err(Error::new(
            ErrorKind::ExpectedAtomCountHeader,
            line,
            parser.current(),
            col,
        ));
    }

    let val_line = parser.next().ok_or_else(|| {
        Error::new(
            ErrorKind::MissingAtomCount,
            String::new(),
            parser.current(),
            1,
        )
    })??;

    let trimmed = val_line.trim();
    if trimmed.is_empty() {
        return Err(Error::new(
            ErrorKind::MissingAtomCount,
            val_line,
            parser.current(),
            1,
        ));
    }

    let col = val_line.find(trimmed).unwrap_or(0) + 1;
    trimmed.parse::<usize>().map_err(|e| {
        Error::new(
            ErrorKind::InvalidAtomCount(e),
            val_line,
            parser.current(),
            col,
        )
    })
}

fn parse_sym_box<B: BufRead>(parser: &mut Parser<B>) -> Result<SymBox> {
    let line = parser.next().ok_or_else(|| {
        Error::new(
            ErrorKind::ExpectedSymboxHeader,
            String::new(),
            parser.current(),
            1,
        )
    })??;

    if !line.trim().starts_with(HEADER_SYM_BOX) {
        let col = line.find(line.trim()).unwrap_or(0) + 1;
        return Err(Error::new(
            ErrorKind::ExpectedSymboxHeader,
            line,
            parser.current(),
            col,
        ));
    }

    let header_end = line.find(HEADER_SYM_BOX).unwrap() + HEADER_SYM_BOX.len();
    let boundaries = line[header_end..].trim().to_string();

    let mut bounds = [[0.0f64; 2]; 3];
    for i in 0..3 {
        let b_line = parser.next().ok_or_else(|| {
            Error::new(
                ErrorKind::MissingSymboxField,
                String::new(),
                parser.current(),
                1,
            )
        })??;

        let tokens = get_tokens(&b_line);
        if tokens.len() < 2 {
            return Err(Error::new(
                ErrorKind::MissingSymboxField,
                b_line,
                parser.current(),
                1,
            ));
        }

        for j in 0..2 {
            let (token, col) = tokens[j];
            bounds[i][j] = token.parse::<f64>().map_err(|e| {
                Error::new(
                    ErrorKind::InvalidSymboxField(e),
                    b_line.clone(),
                    parser.current(),
                    col,
                )
            })?;
        }
    }

    Ok(SymBox {
        boundaries,
        bbox: BoundingBox3::new(
            Point3 {
                coords: [bounds[0][0], bounds[1][0], bounds[2][0]],
            },
            Point3 {
                coords: [bounds[0][1], bounds[1][1], bounds[2][1]],
            },
        ),
    })
}

fn parse_atoms<B: BufRead>(
    parser: &mut Parser<B>,
    count: usize,
) -> Result<(Vec<String>, Vec<Vec<f64>>)> {
    let line = parser.next().ok_or_else(|| {
        Error::new(
            ErrorKind::ExpectedAtomsHeader,
            String::new(),
            parser.current(),
            1,
        )
    })??;

    if !line.trim().starts_with(HEADER_ATOMS) {
        let col = line.find(line.trim()).unwrap_or(0) + 1;
        return Err(Error::new(
            ErrorKind::ExpectedAtomsHeader,
            line,
            parser.current(),
            col,
        ));
    }

    let header_pos = line.find(HEADER_ATOMS).unwrap();
    let header_end = header_pos + HEADER_ATOMS.len();
    let key_tokens = get_tokens(&line[header_end..]);

    if key_tokens.is_empty() {
        return Err(Error::new(
            ErrorKind::MissingAtomKeys,
            line,
            parser.current(),
            header_end + 1,
        ));
    }

    let mut keys = Vec::with_capacity(key_tokens.len());
    let mut seen = HashMap::new();
    for (key, col_in_part) in key_tokens {
        let absolute_col = header_end + col_in_part;
        if seen.insert(key, ()).is_some() {
            return Err(Error::new(
                ErrorKind::DuplicateAtomKeys(key.to_string()),
                line,
                parser.current(),
                absolute_col,
            ));
        }
        keys.push(key.to_string());
    }

    let mut atoms = Vec::with_capacity(count);
    for _ in 0..count {
        let a_line = parser.next().ok_or_else(|| {
            Error::new(
                ErrorKind::MissingAtomRowField,
                String::new(),
                parser.current(),
                1,
            )
        })??;

        let tokens = get_tokens(&a_line);
        if tokens.len() != keys.len() {
            return Err(Error::new(
                ErrorKind::MissingAtomRowField,
                a_line,
                parser.current(),
                1,
            ));
        }

        let mut row = Vec::with_capacity(keys.len());
        for (token, col) in tokens {
            let val = token.parse::<f64>().map_err(|e| {
                Error::new(
                    ErrorKind::InvalidAtomRowField(e),
                    a_line.clone(),
                    parser.current(),
                    col,
                )
            })?;
            row.push(val);
        }
        atoms.push(row);
    }

    Ok((keys, atoms))
}

/// Helper to get tokens with their 1-based column positions
fn get_tokens(line: &str) -> Vec<(&str, usize)> {
    let mut tokens = Vec::new();
    let mut last_pos = 0;
    for token in line.split_whitespace() {
        let pos = line[last_pos..].find(token).unwrap() + last_pos;
        tokens.push((token, pos + 1));
        last_pos = pos + token.len();
    }
    tokens
}

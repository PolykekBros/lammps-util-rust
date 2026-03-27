use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

// Assuming geomutil_util is available as per previous implementation
use geomutil_util::BoundingBox3;

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
            0,
        )
    })??;

    if line.trim() != HEADER_TIMESTEP {
        return Err(Error::new(
            ErrorKind::ExpectedTimestepHeader,
            line,
            parser.current(),
            0,
        ));
    }

    let val_line = parser.next().ok_or_else(|| {
        Error::new(
            ErrorKind::MissingTimestep,
            String::new(),
            parser.current(),
            0,
        )
    })??;

    val_line
        .trim()
        .parse::<u64>()
        .map_err(|e| Error::new(ErrorKind::InvalidTimestep(e), val_line, parser.current(), 0))
}

fn parse_atom_count<B: BufRead>(parser: &mut Parser<B>) -> Result<usize> {
    let line = parser.next().ok_or_else(|| {
        Error::new(
            ErrorKind::ExpectedAtomCountHeader,
            String::new(),
            parser.current(),
            0,
        )
    })??;

    if line.trim() != HEADER_NUM_OF_ATOMS {
        return Err(Error::new(
            ErrorKind::ExpectedAtomCountHeader,
            line,
            parser.current(),
            0,
        ));
    }

    let val_line = parser.next().ok_or_else(|| {
        Error::new(
            ErrorKind::MissingAtomCount,
            String::new(),
            parser.current(),
            0,
        )
    })??;

    val_line.trim().parse::<usize>().map_err(|e| {
        Error::new(
            ErrorKind::InvalidAtomCount(e),
            val_line,
            parser.current(),
            0,
        )
    })
}

fn parse_sym_box<B: BufRead>(parser: &mut Parser<B>) -> Result<SymBox> {
    let line = parser.next().ok_or_else(|| {
        Error::new(
            ErrorKind::ExpectedSymboxHeader,
            String::new(),
            parser.current(),
            0,
        )
    })??;

    if !line.trim().starts_with(HEADER_SYM_BOX) {
        return Err(Error::new(
            ErrorKind::ExpectedSymboxHeader,
            line,
            parser.current(),
            0,
        ));
    }

    let boundaries = line.trim()[HEADER_SYM_BOX.len()..].trim().to_string();

    let mut bounds = [[0.0f64; 2]; 3];
    for i in 0..3 {
        let b_line = parser.next().ok_or_else(|| {
            Error::new(
                ErrorKind::MissingSymboxField,
                String::new(),
                parser.current(),
                0,
            )
        })??;
        let fields: Vec<&str> = b_line.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(Error::new(
                ErrorKind::MissingSymboxField,
                b_line,
                parser.current(),
                0,
            ));
        }
        bounds[i][0] = fields[0].parse::<f64>().map_err(|e| {
            Error::new(
                ErrorKind::InvalidSymboxField(e),
                b_line.clone(),
                parser.current(),
                0,
            )
        })?;
        bounds[i][1] = fields[1].parse::<f64>().map_err(|e| {
            Error::new(
                ErrorKind::InvalidSymboxField(e),
                b_line.clone(),
                parser.current(),
                0,
            )
        })?;
    }

    Ok(SymBox {
        boundaries,
        bbox: BoundingBox3::new(
            [bounds[0][0], bounds[1][0], bounds[2][0]],
            [bounds[0][1], bounds[1][1], bounds[2][1]],
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
            0,
        )
    })??;

    if !line.trim().starts_with(HEADER_ATOMS) {
        return Err(Error::new(
            ErrorKind::ExpectedAtomsHeader,
            line,
            parser.current(),
            0,
        ));
    }

    let keys: Vec<String> = line.trim()[HEADER_ATOMS.len()..]
        .split_whitespace()
        .map(|s| s.to_string())
        .collect();

    if keys.is_empty() {
        return Err(Error::new(
            ErrorKind::MissingAtomKeys,
            line,
            parser.current(),
            0,
        ));
    }

    // Check for duplicate keys
    let mut seen = HashMap::new();
    for key in &keys {
        if seen.insert(key, ()).is_some() {
            return Err(Error::new(
                ErrorKind::DuplicateAtomKeys(key.clone()),
                line,
                parser.current(),
                0,
            ));
        }
    }

    let mut atoms = Vec::with_capacity(count);
    for _ in 0..count {
        let a_line = parser.next().ok_or_else(|| {
            Error::new(
                ErrorKind::MissingAtomRowField,
                String::new(),
                parser.current(),
                0,
            )
        })??;
        let fields: Vec<f64> = a_line
            .split_whitespace()
            .map(|s| s.parse::<f64>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|e| {
                Error::new(
                    ErrorKind::InvalidAtomRowField(e),
                    a_line.clone(),
                    parser.current(),
                    0,
                )
            })?;

        if fields.len() != keys.len() {
            return Err(Error::new(
                ErrorKind::MissingAtomRowField,
                a_line,
                parser.current(),
                0,
            ));
        }
        atoms.push(fields);
    }

    Ok((keys, atoms))
}

use std::{
    collections::HashMap,
    fmt,
    io::{self, BufRead},
    iter,
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
pub struct SnapshotMeta {
    pub timestep: u64,
    pub atoms_count: usize,
    pub sym_box: SymBox,
    pub keys: HashMap<String, usize>,
}

impl SnapshotMeta {
    /// Parses a `SnapshotMeta` from a `Parser`.
    ///
    /// # Errors
    ///
    /// Returns an `Error` if the snapshot meta cannot be parsed.
    pub fn parse<B: BufRead>(parser: &mut Parser<B>) -> Result<Self> {
        let timestep = parse_timestep(parser)?;
        let atoms_count = parse_atom_count(parser)?;
        let sym_box = parse_sym_box(parser)?;
        let keys = parse_keys(parser)?;
        Ok(Self {
            timestep,
            atoms_count,
            sym_box,
            keys,
        })
    }
}

fn parse_timestep<B: BufRead>(parser: &mut Parser<B>) -> Result<u64> {
    let line = parser.next().transpose()?.unwrap_or_default();
    if !line.trim().starts_with(HEADER_TIMESTEP) {
        return Err(Error::new(
            ErrorKind::ExpectedSymboxHeader,
            line,
            parser.current(),
        ));
    }
    let timestep = parser
        .next()
        .ok_or_else(|| Error::new(ErrorKind::MissingTimestep, String::new(), parser.current()))??;
    timestep
        .parse()
        .map_err(|e| Error::new(ErrorKind::InvalidTimestep(e), timestep, parser.current()))
}

fn parse_atom_count<B: BufRead>(parser: &mut Parser<B>) -> Result<usize> {
    let line = parser.next().transpose()?.unwrap_or_default();
    if !line.trim().starts_with(HEADER_NUM_OF_ATOMS) {
        return Err(Error::new(
            ErrorKind::ExpectedSymboxHeader,
            line,
            parser.current(),
        ));
    }
    let atom_count = parser.next().ok_or_else(|| {
        Error::new(ErrorKind::MissingAtomCount, String::new(), parser.current())
    })??;
    atom_count
        .parse()
        .map_err(|e| Error::new(ErrorKind::InvalidAtomCount(e), atom_count, parser.current()))
}

fn parse_sym_box<B: BufRead>(parser: &mut Parser<B>) -> Result<SymBox> {
    let line = parser.next().transpose()?.unwrap_or_default();
    if !line.trim().starts_with(HEADER_SYM_BOX) {
        return Err(Error::new(
            ErrorKind::ExpectedSymboxHeader,
            line,
            parser.current(),
        ));
    }
    let boundaries = line[HEADER_SYM_BOX.len()..].trim().to_owned();
    let mut lo = Point3::default();
    let mut hi = Point3::default();
    for i in 0..3 {
        let line = parser.next().transpose()?.unwrap_or_default();
        let pair = line
            .split_whitespace()
            .map(|s| {
                s.parse::<f64>().map_err(|e| {
                    let line = line.clone();
                    Error::new(ErrorKind::InvalidSymboxField(e), line, parser.current())
                })
            })
            .collect::<Result<Vec<_>>>()?;
        if pair.len() < 2 {
            return Err(Error::new(
                ErrorKind::MissingSymboxField,
                line,
                parser.current(),
            ));
        }
        lo[i] = pair[0];
        hi[i] = pair[1];
    }
    Ok(SymBox {
        boundaries,
        bbox: BoundingBox3::new(lo, hi),
    })
}

fn parse_keys<B: BufRead>(parser: &mut Parser<B>) -> Result<HashMap<String, usize>> {
    let line = parser.next().transpose()?.unwrap_or_default();
    if !line.trim().starts_with(HEADER_ATOMS) {
        return Err(Error::new(
            ErrorKind::ExpectedAtomsHeader,
            line,
            parser.current(),
        ));
    }
    let tokens = line[HEADER_ATOMS.len()..].trim();
    let mut keys = HashMap::new();
    for (idx, key) in tokens.split_whitespace().enumerate() {
        if keys.contains_key(key) {
            return Err(Error::new(
                ErrorKind::DuplicateAtomKeys(key.to_owned()),
                line.clone(),
                parser.current(),
            ));
        }
        keys.insert(key.to_owned(), idx);
    }
    if keys.is_empty() {
        return Err(Error::new(
            ErrorKind::MissingAtomKeys,
            line.clone(),
            parser.current(),
        ));
    }
    Ok(keys)
}

#[derive(Clone)]
pub struct Snapshot {
    meta: SnapshotMeta,
    keys: Vec<String>,
    atoms: Vec<f64>,
}

impl Snapshot {
    #[must_use]
    pub fn new(meta: SnapshotMeta) -> Self {
        let atoms = vec![0.0; meta.atoms_count * meta.keys.len()];
        let mut keys: Vec<(&String, &usize)> = meta.keys.iter().collect();
        keys.sort_by(|a, b| a.1.cmp(b.1));
        let keys = keys.into_iter().map(|i| i.0.to_owned()).collect();
        Self { meta, keys, atoms }
    }

    /// Parses a `Snapshot` from a `Parser` given a `SnapshotMeta`.
    ///
    /// # Errors
    ///
    /// Returns an `Error` if the snapshot cannot be parsed.
    pub fn parse<B: BufRead>(parser: &mut Parser<B>, meta: SnapshotMeta) -> Result<Self> {
        let mut snapshot = Self::new(meta);
        for i in 0..snapshot.get_atoms_count() {
            let line = parser.next().transpose()?.unwrap_or_default();
            let values = line
                .split_whitespace()
                .map(|s| {
                    s.parse::<f64>().map_err(|e| {
                        let line = line.clone();
                        Error::new(
                            ErrorKind::InvalidAtomRowField(e),
                            line.clone(),
                            parser.current(),
                        )
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            if values.len() != snapshot.get_keys().len() {
                return Err(Error::new(
                    ErrorKind::MissingAtomRowField,
                    line.clone(),
                    parser.current(),
                ));
            }
            for (j, val) in values.into_iter().enumerate() {
                snapshot.set_atom_value(j, i, val);
            }
        }
        Ok(snapshot)
    }

    /// Writes the `Snapshot` to a `Write`.
    ///
    /// # Errors
    ///
    /// Returns an `io::Error` if the write operation fails.
    pub fn write<W>(&self, w: &mut W) -> io::Result<()>
    where
        W: io::Write,
    {
        writeln!(w, "{HEADER_TIMESTEP}")?;
        writeln!(w, "{}", self.get_timestep())?;
        writeln!(w, "{HEADER_NUM_OF_ATOMS}")?;
        writeln!(w, "{}", self.get_atoms_count())?;
        writeln!(w, "{HEADER_SYM_BOX} {}", self.get_symbox().boundaries)?;
        for (lo, hi) in iter::zip(self.get_symbox().bbox.lower, self.get_symbox().bbox.upper) {
            writeln!(w, "{lo} {hi}",)?;
        }
        writeln!(w, "{HEADER_ATOMS} {}", self.get_keys().join(" "))?;
        for i in 0..self.get_atoms_count() {
            write!(w, "{}", self.atoms[i])?;
            for j in 1..self.keys.len() {
                write!(w, " {}", self.get_atom_value(j, i))?;
            }
            writeln!(w)?;
        }
        Ok(())
    }

    #[must_use]
    pub fn get_atoms_count(&self) -> usize {
        self.meta.atoms_count
    }

    #[must_use]
    pub fn get_timestep(&self) -> u64 {
        self.meta.timestep
    }

    #[must_use]
    pub fn get_symbox(&self) -> &SymBox {
        &self.meta.sym_box
    }

    #[must_use]
    pub const fn get_keys_map(&self) -> &HashMap<String, usize> {
        &self.meta.keys
    }

    #[must_use]
    pub fn get_keys(&self) -> &[String] {
        &self.keys
    }

    #[must_use]
    pub fn get_property_index(&self, key: &str) -> usize {
        self.meta.keys[key]
    }

    #[must_use]
    pub fn get_property(&self, key: &str) -> &[f64] {
        let start = self.meta.keys[key] * self.get_atoms_count();
        let end = start + self.get_atoms_count();
        &self.atoms[start..end]
    }

    #[must_use]
    pub fn get_property_mut(&mut self, key: &str) -> &mut [f64] {
        let start = self.meta.keys[key] * self.get_atoms_count();
        let end = start + self.get_atoms_count();
        &mut self.atoms[start..end]
    }

    #[must_use]
    pub fn get_atom_value(&self, property_index: usize, atom_index: usize) -> f64 {
        let idx = self.get_atoms_count() * property_index + atom_index;
        self.atoms[idx]
    }

    pub fn set_atom_value(&mut self, property_index: usize, atom_index: usize, value: f64) {
        let idx = self.get_atoms_count() * property_index + atom_index;
        self.atoms[idx] = value;
    }

    pub fn get_zero_lvl(&self) -> f64 {
        self.get_property("z")
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    // #[must_use]
    // pub fn get_coordinates(&self) -> Vec<XYZ> {
    //     izip!(
    //         self.get_property("x").iter().copied(),
    //         self.get_property("y").iter().copied(),
    //         self.get_property("z").iter().copied(),
    //     )
    //     .enumerate()
    //     .map(|(i, xyz)| XYZ::new(Into::<[f64; 3]>::into(xyz).into(), i, false))
    //     .collect()
    // }
}

impl fmt::Debug for Snapshot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DumpSnapshot")
            .field("timestep", &self.get_timestep())
            .field("atoms_count", &self.get_atoms_count())
            .field("keys", &self.get_keys())
            .finish_non_exhaustive()
    }
}

#[must_use]
pub fn copy_snapshot(input_snapshot: &Snapshot) -> Snapshot {
    copy_snapshot_with_indices_with_keys(
        input_snapshot,
        std::iter::empty(),
        0..input_snapshot.get_atoms_count(),
    )
}

#[must_use]
pub fn copy_snapshot_with_indices(
    input_snapshot: &Snapshot,
    indices: impl Iterator<Item = usize>,
) -> Snapshot {
    copy_snapshot_with_indices_with_keys(input_snapshot, std::iter::empty(), indices)
}

#[must_use]
pub fn copy_snapshot_with_keys<'a>(
    input_snapshot: &Snapshot,
    additional_keys: impl Iterator<Item = &'a str>,
) -> Snapshot {
    copy_snapshot_with_indices_with_keys(
        input_snapshot,
        additional_keys,
        0..input_snapshot.get_atoms_count(),
    )
}

#[must_use]
pub fn copy_snapshot_with_indices_with_keys<'a, 'b>(
    input_snapshot: &Snapshot,
    additional_keys: impl Iterator<Item = &'a str>,
    indices: impl Iterator<Item = usize>,
) -> Snapshot {
    let mut input_meta = input_snapshot.meta.clone();
    for key in additional_keys {
        input_meta
            .keys
            .insert(key.to_string(), input_meta.keys.len());
    }
    let indices = indices.collect::<Vec<_>>();
    let mut snapshot = Snapshot::new(input_meta);
    for (new_i, i) in indices.into_iter().enumerate() {
        for (j, _) in input_snapshot.get_keys().iter().enumerate() {
            snapshot.set_atom_value(j, new_i, input_snapshot.get_atom_value(j, i));
        }
    }
    snapshot
}

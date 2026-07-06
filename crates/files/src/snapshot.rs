use std::{
    collections::HashMap,
    fmt,
    io::{self, BufRead},
    iter,
};

// Assuming geomutil_util is available as per previous implementation
use geomutil_util::{BoundingBox3, Point3};

use crate::{
    error::{ErrorKind, Result},
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
    parser.parse_line(|line| {
        line.starts_with(HEADER_TIMESTEP)
            .then_some(())
            .ok_or(ErrorKind::ExpectedTimestepHeader)
    })?;
    parser.parse_line(|line| line.parse().map_err(ErrorKind::InvalidTimestep))
}

fn parse_atom_count<B: BufRead>(parser: &mut Parser<B>) -> Result<usize> {
    parser.parse_line(|line| {
        line.starts_with(HEADER_NUM_OF_ATOMS)
            .then_some(())
            .ok_or(ErrorKind::ExpectedAtomCountHeader)
    })?;
    parser.parse_line(|line| line.parse().map_err(ErrorKind::InvalidAtomCount))
}

fn parse_sym_box<B: BufRead>(parser: &mut Parser<B>) -> Result<SymBox> {
    let boundaries = parser.parse_line(|line| {
        line.starts_with(HEADER_SYM_BOX)
            .then(|| line[HEADER_SYM_BOX.len()..].trim().to_owned())
            .ok_or(ErrorKind::ExpectedSymboxHeader)
    })?;
    let mut lo = Point3::default();
    let mut hi = Point3::default();
    for i in 0..3 {
        let pair: Vec<f64> = parser.parse_line(|line| {
            let pair = line
                .split_whitespace()
                .map(|s| s.parse().map_err(ErrorKind::InvalidSymboxField))
                .collect::<std::result::Result<Vec<_>, _>>()?;
            (pair.len() == 2)
                .then_some(pair)
                .ok_or(ErrorKind::MissingSymboxField)
        })?;
        lo[i] = pair[0];
        hi[i] = pair[1];
    }
    Ok(SymBox {
        boundaries,
        bbox: BoundingBox3::new(lo, hi),
    })
}

fn parse_keys<B: BufRead>(parser: &mut Parser<B>) -> Result<HashMap<String, usize>> {
    parser.parse_line(|line| {
        line.starts_with(HEADER_ATOMS)
            .then_some(())
            .ok_or(ErrorKind::ExpectedAtomsHeader)?;
        let tokens = line[HEADER_ATOMS.len()..].trim();
        let keys = tokens.split_whitespace().enumerate().try_fold(
            HashMap::new(),
            |mut map, (idx, key)| {
                map.insert(key.to_owned(), idx)
                    .is_none()
                    .then_some(map)
                    .ok_or_else(|| ErrorKind::DuplicateAtomKeys(key.to_owned()))
            },
        )?;
        (!keys.is_empty())
            .then_some(keys)
            .ok_or(ErrorKind::MissingAtomKeys)
    })
}

fn parse_atom_row<B: BufRead>(parser: &mut Parser<B>, n: usize) -> Result<Vec<f64>> {
    parser.parse_line(|line| {
        let values = line
            .split_whitespace()
            .map(|s| s.parse::<f64>().map_err(ErrorKind::InvalidAtomRowField))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        (values.len() == n)
            .then_some(values)
            .ok_or(ErrorKind::MissingAtomRowField)
    })
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
            let values = parse_atom_row(parser, snapshot.get_keys().len())?;
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
            writeln!(w, "{lo} {hi}")?;
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

    fn get_property_atoms_range(&self, key: &str) -> (usize, usize) {
        let idx = self.get_property_index(key);
        let start = idx * self.get_atoms_count();
        let end = start + self.get_atoms_count();
        (start, end)
    }

    #[must_use]
    pub fn get_property(&self, key: &str) -> &[f64] {
        let (start, end) = self.get_property_atoms_range(key);
        &self.atoms[start..end]
    }

    #[must_use]
    pub fn get_property_mut(&mut self, key: &str) -> &mut [f64] {
        let (start, end) = self.get_property_atoms_range(key);
        &mut self.atoms[start..end]
    }

    #[must_use]
    pub fn get_property_index_checked(&self, key: &str) -> Option<usize> {
        self.meta.keys.get(key).copied()
    }

    fn get_property_atoms_range_checked(&self, key: &str) -> Option<(usize, usize)> {
        let idx = self.get_property_index_checked(key)?;
        let start = idx * self.get_atoms_count();
        let end = start + self.get_atoms_count();
        Some((start, end))
    }

    #[must_use]
    pub fn get_property_checked(&self, key: &str) -> Option<&[f64]> {
        let (start, end) = self.get_property_atoms_range_checked(key)?;
        Some(&self.atoms[start..end])
    }

    #[must_use]
    pub fn get_property_mut_checked(&mut self, key: &str) -> Option<&mut [f64]> {
        let (start, end) = self.get_property_atoms_range_checked(key)?;
        Some(&mut self.atoms[start..end])
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

    #[must_use]
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
pub fn copy_snapshot_with_indices_with_keys<'a>(
    input_snapshot: &Snapshot,
    additional_keys: impl Iterator<Item = &'a str>,
    indices: impl Iterator<Item = usize>,
) -> Snapshot {
    let indices = indices.collect::<Vec<_>>();
    let mut input_meta = input_snapshot.meta.clone();
    input_meta.atoms_count = indices.len();
    for key in additional_keys {
        input_meta
            .keys
            .insert(key.to_string(), input_meta.keys.len());
    }
    let mut snapshot = Snapshot::new(input_meta);
    for (new_i, i) in indices.into_iter().enumerate() {
        for (j, _) in input_snapshot.get_keys().iter().enumerate() {
            snapshot.set_atom_value(j, new_i, input_snapshot.get_atom_value(j, i));
        }
    }
    snapshot
}

#[cfg(test)]
mod tests {
    use super::*;
    // Note: Adjust the ErrorKind imports depending on your exact error module structure
    use crate::error::ErrorKind;

    // A realistic mock of the snapshot format
    const VALID_DUMP: &str = "\
ITEM: TIMESTEP
1000
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 20.0
-5.0 5.0
ITEM: ATOMS id type x y z
1 1 0.5 1.5 2.5
2 2 8.5 9.5 -2.5
";

    #[test]
    fn test_parse_valid_snapshot() {
        // `as_bytes()` implements `BufRead`, making it perfect for our parser
        let mut parser = Parser::new(VALID_DUMP.as_bytes());

        // 1. Test Metadata Parsing
        let meta = SnapshotMeta::parse(&mut parser);
        assert!(meta.is_ok());
        let meta = meta.unwrap();
        assert_eq!(meta.timestep, 1000);
        assert_eq!(meta.atoms_count, 2);
        assert_eq!(meta.sym_box.boundaries, "pp pp pp");

        // 2. Test Snapshot Parsing
        let snapshot = Snapshot::parse(&mut parser, meta);
        assert!(snapshot.is_ok());
        let snapshot = snapshot.unwrap();
        assert_eq!(snapshot.get_atoms_count(), 2);
        assert_eq!(snapshot.get_keys(), &["id", "type", "x", "y", "z"]);

        // 3. Test flattened array access (Struct-of-Arrays math)
        let id_idx = snapshot.get_property_index("id");
        let x_idx = snapshot.get_property_index("x");

        // Atom 0: id = 1.0, x = 0.5
        assert_eq!(snapshot.get_atom_value(id_idx, 0), 1.0);
        assert_eq!(snapshot.get_atom_value(x_idx, 0), 0.5);

        // Atom 1: id = 2.0, x = 8.5
        assert_eq!(snapshot.get_atom_value(id_idx, 1), 2.0);
        assert_eq!(snapshot.get_atom_value(x_idx, 1), 8.5);
    }

    #[test]
    fn test_parse_missing_atom_columns() {
        const INVALID_DUMP: &str = "\
ITEM: TIMESTEP
1000
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 20.0
-5.0 5.0
ITEM: ATOMS id type x y z
1 1 0.5 1.5
2 2 8.5 9.5 -2.5
";
        let mut parser = Parser::new(INVALID_DUMP.as_bytes());
        let meta = SnapshotMeta::parse(&mut parser).unwrap();

        let result = Snapshot::parse(&mut parser, meta);
        assert!(result.is_err());

        let err = result.unwrap_err();
        assert!(matches!(err.kind(), ErrorKind::MissingAtomRowField));
        assert_eq!(err.line(), 10); // Should fail exactly on line 9
    }

    #[test]
    fn test_copy_snapshot_subset() {
        let mut parser = Parser::new(VALID_DUMP.as_bytes());
        let meta = SnapshotMeta::parse(&mut parser).unwrap();
        let original_snapshot = Snapshot::parse(&mut parser, meta).unwrap();

        // Copy ONLY the second atom (index 1), and add a new property "vx"
        let new_keys = vec!["vx"];
        let indices = vec![1];

        let copied_snapshot = copy_snapshot_with_indices_with_keys(
            &original_snapshot,
            new_keys.into_iter(),
            indices.into_iter(),
        );

        // Verify the atoms count was correctly updated
        assert_eq!(copied_snapshot.get_atoms_count(), 1);

        // Verify keys were appended
        assert_eq!(
            copied_snapshot.get_keys(),
            &["id", "type", "x", "y", "z", "vx"]
        );

        // Verify the data was copied correctly to the new index 0
        let id_idx = copied_snapshot.get_property_index("id");
        let x_idx = copied_snapshot.get_property_index("x");
        let vx_idx = copied_snapshot.get_property_index("vx");

        assert_eq!(copied_snapshot.get_atom_value(id_idx, 0), 2.0); // original id was 2
        assert_eq!(copied_snapshot.get_atom_value(x_idx, 0), 8.5); // original x was 8.5

        // Verify new keys are initialized to 0.0
        assert_eq!(copied_snapshot.get_atom_value(vx_idx, 0), 0.0);
    }
}

use core::fmt;
use std::collections::HashMap;
use std::fs;
use std::io;
use std::path::Path;

// use memmap2::Mmap;
const HEADER_TIMESTEP: &str = "ITEM: TIMESTEP";
const HEADER_NUM_OF_ATOMS: &str = "ITEM: NUMBER OF ATOMS";
const HEADER_SYM_BOX: &str = "ITEM: BOX";
const HEADER_ATOMS: &str = "ITEM: ATOMS";

pub struct DumpSnapshot {
    pub step: u64,
    pub atoms_count: usize,
    keys: HashMap<String, usize>,
    atoms: Vec<f64>,
}

pub struct DumpFile {
    snapshots: HashMap<u64, DumpSnapshot>,
}

#[derive(Debug)]
pub enum DumpParsingError {
    InvalidOrMissingTimestep,
    InvalidOrMissingNumberOfAtoms,
    MissingSymBox,
    MissingAtomKeys,
    DuplicateAtomKeys,
    DuplicateSnapshots,
    InvalidOrMissingAtomRow,
    IO(io::Error),
}

impl std::fmt::Display for DumpParsingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl std::error::Error for DumpParsingError {}

impl DumpFile {
    pub fn new(path: &Path, timesteps: &[u64]) -> Result<Self, DumpParsingError> {
        let lines = fs::read_to_string(path).map_err(DumpParsingError::IO)?;
        let mut lines = lines.trim().split("\n");

        let mut timesteps = timesteps.to_vec();
        timesteps.sort();

        let mut dump = Self {
            snapshots: HashMap::new(),
        };

        'read_all: loop {
            let timestep = match (lines.next(), lines.next().map(str::parse::<u64>)) {
                (Some(HEADER_TIMESTEP), Some(Ok(n))) => n,
                (None, _) => break Ok(dump),
                (_, _) => break Err(DumpParsingError::InvalidOrMissingTimestep),
            };
            let number_of_atoms = match (lines.next(), lines.next().map(str::parse::<usize>)) {
                (Some(HEADER_NUM_OF_ATOMS), Some(Ok(n))) => n,
                (_, _) => break Err(DumpParsingError::InvalidOrMissingNumberOfAtoms),
            };
            if !timesteps.is_empty() {
                if &timestep > timesteps.last().unwrap_or(&u64::MAX) {
                    break Ok(dump);
                } else if !timesteps.contains(&timestep) {
                    for _ in 0..number_of_atoms + 4 + 1 {
                        if lines.next().is_none() {
                            break;
                        }
                    }
                    continue;
                }
            }
            if dump.snapshots.contains_key(&timestep) {
                break Err(DumpParsingError::DuplicateSnapshots);
            }
            match lines
                .next()
                .and_then(|l| l.split_at_checked(HEADER_SYM_BOX.len()))
            {
                Some((HEADER_SYM_BOX, _)) => {
                    // TODO: read box parameters
                    lines.next();
                    lines.next();
                    lines.next();
                }
                _ => break Err(DumpParsingError::MissingSymBox),
            }
            let mut snapshot_keys = HashMap::new();
            match lines
                .next()
                .and_then(|l| l.split_at_checked(HEADER_ATOMS.len()))
            {
                Some((HEADER_ATOMS, keys)) => {
                    for key in keys.split_whitespace() {
                        if snapshot_keys
                            .insert(key.to_string(), snapshot_keys.len())
                            .is_some()
                        {
                            break 'read_all Err(DumpParsingError::DuplicateAtomKeys);
                        }
                    }
                }
                _ => break Err(DumpParsingError::MissingAtomKeys),
            }
            let mut snapshot = DumpSnapshot {
                step: timestep,
                atoms_count: number_of_atoms,
                atoms: vec![0.0; number_of_atoms * snapshot_keys.len()],
                keys: snapshot_keys,
            };
            for i in 0..number_of_atoms {
                let Some::<Vec<f64>>(values) = lines
                    .next()
                    .map(str::split_whitespace)
                    .and_then(|tokens| tokens.map(str::parse::<f64>).map(Result::ok).collect())
                else {
                    break 'read_all Err(DumpParsingError::InvalidOrMissingAtomRow);
                };
                for (j, val) in values.iter().enumerate() {
                    snapshot.atoms[number_of_atoms * j + i] = *val;
                }
            }
            dump.snapshots.insert(snapshot.step, snapshot);
        }
    }

    pub fn get_snapshots(&self) -> Vec<&DumpSnapshot> {
        let mut entries: Vec<(&u64, &DumpSnapshot)> = self.snapshots.iter().collect();
        entries.sort_by(|a, b| a.0.cmp(b.0));
        entries.into_iter().map(|i| i.1).collect()
    }

    pub fn get_property(&self, timestep: u64, key: &str) -> &[f64] {
        self.snapshots[&timestep].get_property(key)
    }
}

impl DumpSnapshot {
    pub fn get_keys(&self) -> Vec<&String> {
        let mut entries: Vec<(&String, &usize)> = self.keys.iter().collect();
        entries.sort_by(|a, b| a.1.cmp(b.1));
        entries.into_iter().map(|i| i.0).collect()
    }

    pub fn get_property(&self, key: &str) -> &[f64] {
        let start = self.keys[key] * self.atoms_count;
        let end = start + self.atoms_count;
        &self.atoms[start..end]
    }
}

impl fmt::Debug for DumpSnapshot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DumpSnapshot")
            .field("step", &self.step)
            .field("atoms_count", &self.atoms_count)
            .field("keys", &self.keys)
            .finish()
    }
}

pub fn get_zero_lvl(snap: &DumpSnapshot) -> f64 {
    snap.get_property("z")
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max)
}

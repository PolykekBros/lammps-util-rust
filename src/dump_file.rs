use std::collections::HashMap;
use std::fs;
use std::io;
use std::path::Path;

use crate::dump_snapshot::{DumpSnapshot, HEADER_NUM_OF_ATOMS, HEADER_TIMESTEP};

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

        loop {
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
            let snapshot = DumpSnapshot::read(&mut lines, timestep, number_of_atoms)?;
            dump.snapshots.insert(snapshot.step, snapshot);
        }
    }

    pub fn save(&self, path: &Path) -> io::Result<()> {
        let f = fs::File::create(path)?;
        let mut w = io::BufWriter::new(f);
        for snapshot in self.get_snapshots() {
            snapshot.write(&mut w)?;
        }
        Ok(())
    }

    pub fn get_snapshots(&self) -> Vec<&DumpSnapshot> {
        let mut entries: Vec<(&u64, &DumpSnapshot)> = self.snapshots.iter().collect();
        entries.sort_by(|a, b| a.0.cmp(b.0));
        entries.into_iter().map(|i| i.1).collect()
    }

    #[inline]
    pub fn get_property(&self, timestep: u64, key: &str) -> &[f64] {
        self.snapshots[&timestep].get_property(key)
    }
}

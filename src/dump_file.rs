use anyhow::{anyhow, Context, Result};
use std::fs;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::{collections::HashMap, io::BufReader};

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
        write!(f, "{self:?}")
    }
}

impl std::error::Error for DumpParsingError {}

impl DumpFile {
    #[must_use] pub fn new(snapshots: Vec<DumpSnapshot>) -> Self {
        let mut snapshots_map = HashMap::new();
        for snapshot in snapshots {
            snapshots_map.entry(snapshot.step).insert_entry(snapshot);
        }
        Self {
            snapshots: snapshots_map,
        }
    }

    pub fn read(path: &Path, timesteps: &[u64]) -> Result<Self> {
        let mut lines = BufReader::new(
            File::open(path).context(format!("Reading {}", path.to_string_lossy()))?,
        )
        .lines()
        .map_while(Result::ok);
        let mut timesteps = timesteps.to_vec();
        timesteps.sort_unstable();

        let mut dump = Self {
            snapshots: HashMap::new(),
        };

        loop {
            let timestep = match (
                lines.next().filter(|s| s == HEADER_TIMESTEP),
                lines.next().map(|s| s.as_str().parse::<u64>()),
            ) {
                (Some(_), Some(Ok(n))) => n,
                (None, _) => break Ok(dump),
                (_, _) => break Err(anyhow!(DumpParsingError::InvalidOrMissingTimestep)),
            };
            let number_of_atoms = match lines
                .next()
                .filter(|s| s == HEADER_NUM_OF_ATOMS)
                .zip(lines.next().map(|s| s.as_str().parse::<usize>()))
            {
                Some((_, Ok(n))) => n,
                _ => break Err(anyhow!(DumpParsingError::InvalidOrMissingNumberOfAtoms)),
            };
            if !timesteps.is_empty() {
                if &timestep > timesteps.last().unwrap_or(&u64::MAX) {
                    break Ok(dump);
                } else if !timesteps.contains(&timestep) {
                    for _ in 0..=(number_of_atoms + 4) {
                        if lines.next().is_none() {
                            break;
                        }
                    }
                    continue;
                }
            }
            if dump.snapshots.contains_key(&timestep) {
                break Err(anyhow!(DumpParsingError::DuplicateSnapshots));
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

    #[must_use] pub fn get_snapshots(&self) -> Vec<&DumpSnapshot> {
        let mut entries: Vec<(&u64, &DumpSnapshot)> = self.snapshots.iter().collect();
        entries.sort_by(|a, b| a.0.cmp(b.0));
        entries.into_iter().map(|i| i.1).collect()
    }

    #[inline]
    #[must_use] pub fn get_property(&self, timestep: u64, key: &str) -> &[f64] {
        self.snapshots[&timestep].get_property(key)
    }
}

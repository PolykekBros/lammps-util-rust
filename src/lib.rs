use std::collections::HashMap;
use std::fs;
use std::io;
use std::path::Path;

// use memmap2::Mmap;

pub struct DumpTimestep {
    start: usize,
    pub atom_count: usize,
}

pub struct Dump {
    atoms: Vec<f64>,
    pub timesteps: HashMap<u64, DumpTimestep>,
    keys: HashMap<String, usize>,
}

#[derive(Debug)]
pub enum DumpParsingError {
    NoTimestep,
    NoNumberOfAtoms,
    DuplicateKeys,
    DuplicateTimesteps,
    IO(io::Error),
    InvalidNumber,
    MissingAtomRow,
    InvalidAtomRow,
}

impl Dump {
    pub fn new(path: &Path, timesteps: &[u64]) -> Result<Self, DumpParsingError> {
        let lines = fs::read_to_string(path).map_err(DumpParsingError::IO)?;
        let mut lines = lines.split("\n");

        let mut dump = Self {
            atoms: Vec::new(),
            timesteps: HashMap::new(),
            keys: HashMap::new(),
        };

        let mut timestep: Option<u64> = None;
        let mut atom_count: Option<usize> = None;
        while let Some(line) = lines.next() {
            if line == "ITEM: TIMESTEP" {
                timestep = Some(
                    lines
                        .next()
                        .ok_or(DumpParsingError::NoTimestep)?
                        .parse::<u64>()
                        .map_err(|_| DumpParsingError::InvalidNumber)?,
                );
            } else if line == "ITEM: NUMBER OF ATOMS" {
                atom_count = Some(
                    lines
                        .next()
                        .ok_or(DumpParsingError::NoNumberOfAtoms)?
                        .parse::<usize>()
                        .map_err(|_| DumpParsingError::InvalidNumber)?,
                );
            } else if line.contains("ITEM: ATOMS") {
                if dump.keys.is_empty() {
                    let keys = line.split_whitespace().skip(2);
                    for key in keys {
                        if dump.keys.contains_key(key) {
                            return Err(DumpParsingError::DuplicateKeys);
                        }
                        dump.keys.insert(key.to_string(), dump.keys.len());
                    }
                }
                let Some(timestep) = timestep else {
                    return Err(DumpParsingError::NoTimestep);
                };
                let Some(atom_count) = atom_count else {
                    return Err(DumpParsingError::NoNumberOfAtoms);
                };
                if &timestep > timesteps.last().unwrap_or(&u64::MAX) {
                    break;
                }
                if timesteps.is_empty() || timesteps.contains(&timestep) {
                    if dump.timesteps.contains_key(&timestep) {
                        return Err(DumpParsingError::DuplicateTimesteps);
                    }
                    dump.timesteps.insert(
                        timestep,
                        DumpTimestep {
                            start: atom_count * dump.timesteps.len() * dump.keys.len(),
                            atom_count,
                        },
                    );
                    let mut atoms = vec![0.0; atom_count * dump.keys.len()];
                    for i in 0..atom_count {
                        let Some(tokens) = lines.next() else {
                            return Err(DumpParsingError::MissingAtomRow);
                        };
                        let tokens: Vec<&str> = tokens.split_whitespace().collect();
                        if tokens.len() != dump.keys.len() {
                            return Err(DumpParsingError::InvalidAtomRow);
                        }
                        for j in dump.keys.values() {
                            atoms[atom_count * j + i] = tokens[*j]
                                .parse::<f64>()
                                .map_err(|_| DumpParsingError::InvalidNumber)?
                        }
                    }
                    dump.atoms.extend_from_slice(&atoms);
                }
            }
        }

        Ok(dump)
    }

    pub fn get_property(&self, timestep: u64, key: &str) -> &[f64] {
        let tstep = &self.timesteps[&timestep];
        let start = tstep.start + self.keys[key] * tstep.atom_count;
        let end = start + tstep.atom_count;
        &self.atoms[start..end]
    }

    pub fn get_keys(&self) -> Vec<&String> {
        let mut entries: Vec<(&String, &usize)> = self.keys.iter().collect();
        entries.sort_by(|a, b| a.1.cmp(b.1));
        entries.into_iter().map(|i| i.0).collect()
    }
}

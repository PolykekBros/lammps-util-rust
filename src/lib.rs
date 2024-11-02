use std::collections::HashMap;
use std::fs;
use std::io;
use std::path::Path;

// use memmap2::Mmap;
const HEADER_TIMESTEP = "ITEM: TIMESTEP";
const HEADER_NUM_OF_ATOMS = "ITEM: NUMBER OF ATOMS";
const HEADER_SYM_BOX = "ITEM: BOX";
const HEADER_ATOMS = "ITEM: ATOMS";

pub struct DumpTimestep {
    start: usize,
    pub number: u64,
    pub atom_count: usize,
}

pub struct Dump {
    atoms: Vec<f64>,
    pub timesteps: HashMap<u64, DumpTimestep>,
    keys: HashMap<String, usize>,
}

#[derive(Debug)]
pub enum DumpParsingError {
    InvalidOrMissingTimestep,
    InvalidOrMissingNumberOfAtoms,
    MissingSymBox,
    MissingAtomKeys,
    DuplicateKeys,
    DuplicateTimesteps,
    MissingAtomRow,
    InvalidAtomRow,
    IO(io::Error),
}

impl Dump {
    pub fn new(path: &Path, timesteps: &[u64]) -> Result<Self, DumpParsingError> {
        let lines = fs::read_to_string(path).map_err(DumpParsingError::IO)?;
        let mut lines = lines.trim().split("\n");

        let mut timesteps = timesteps.to_vec();
        timesteps.sort();

        let mut dump = Self {
            atoms: Vec::new(),
            timesteps: HashMap::new(),
            keys: HashMap::new(),
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
            match lines.next().map(|l| l.split_at_checked(HEADER_SYM_BOX.len())).flatten() {
                Some((HEADER_SYM_BOX, _)) => {
                    // TODO: read box parameters
                    lines.next();
                    lines.next();
                    lines.next();
                }
                _ => break Err(DumpParsingError::MissingSymBox),
            }
            match lines.next().map(|l| l.split_at_checked(HEADER_ATOMS.len())).flatten() {
                Some((HEADER_ATOMS, keys)) => {
                    if dump.keys.is_empty() {
                        for key in keys.split_whitespace() {
                            if dump.keys.insert(key.to_string(), dump.keys.len()).is_some() {
                                break 'read_all Err(DumpParsingError::DuplicateKeys);
                            }
                        }
                    }
                }
                _ => break Err(DumpParsingError::MissingAtomKeys),
            }
            if dump
                .timesteps
                .insert(
                    timestep,
                    DumpTimestep {
                        start: number_of_atoms * dump.timesteps.len() * dump.keys.len(),
                        number: timestep,
                        atom_count: number_of_atoms,
                    },
                )
                .is_some()
            {
                break Err(DumpParsingError::DuplicateKeys);
            }

            let mut atoms = vec![0.0; number_of_atoms * dump.keys.len()];
            for i in 0..number_of_atoms {
                let Some(tokens) = lines.next() else {
                    break 'read_all Err(DumpParsingError::MissingAtomRow);
                };
                let tokens: Vec<&str> = tokens.split_whitespace().collect();
                if tokens.len() != dump.keys.len() {
                    return Err(DumpParsingError::InvalidAtomRow);
                }
                for j in dump.keys.values() {
                    atoms[number_of_atoms * j + i] = tokens[*j]
                        .parse::<f64>()
                        .map_err(|_| DumpParsingError::InvalidAtomRow)?
                }
            }
            dump.atoms.extend_from_slice(&atoms);
        }
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

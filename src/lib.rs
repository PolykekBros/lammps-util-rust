use std::collections::HashMap;
use std::fs;
use std::path::Path;

// use memmap2::Mmap;

struct DumpTimestep {
    start: usize,
    atom_count: usize,
}

struct Dump {
    atoms: Vec<f64>,
    timesteps: HashMap<u64, DumpTimestep>,
    keys: HashMap<String, usize>,
}

enum DumpParsingError {
    NoTimestep,
    NoNumberOfAtoms,
    DuplicateKeys,
    DuplicateTimesteps,
}

impl Dump {
    pub fn new(path: &Path, timesteps: &Vec<u64>) -> Result<Self, DumpParsingError> {
        let mut lines = fs::read_to_string(path)?.split("\n");

        let mut dump = Self {
            atoms: Vec::new(),
            timesteps: HashMap::new(),
            keys: HashMap::new(),
        };

        let mut timestep: Option<u64> = None;
        let mut atom_count: Option<usize> = None;
        while let Some(line) = lines.next() {
            if line == "ITEM: TIMESTEP" {
                let timestep = Some(lines.next().ok_or(DumpParsingError::NoTimestep)?.parse()?);
            } else if line == "ITEM: NUMBER OF ATOMS" {
                let atom_count = Some(
                    lines
                        .next()
                        .ok_or(DumpParsingError::NoNumberOfAtoms)?
                        .parse()?,
                );
            } else if line.contains("ITEM: ATOMS") {
                if dump.keys.is_empty() {
                    let keys = line.split_whitespace().skip(2);
                    for key in keys {
                        if dump.keys.contains_key(key) {
                            return Err(DumpParsingError::DuplicateKeys);
                        }
                        dump.keys[key] = dump.keys.len();
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
                if timesteps.len() == 0 || timesteps.contains(&timestep) {
                    if dump.timesteps.contains_key(&timestep) {
                        return Err(DumpParsingError::DuplicateTimesteps);
                    }
                    dump.timesteps[&timestep] = DumpTimestep {
                        start: atom_count * dump.timesteps.len() * dump.keys.len(),
                        atom_count,
                    };
                    ...
                }
            }
        }

        Ok(dump)
    }
}

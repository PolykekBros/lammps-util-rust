use std::collections::HashMap;
use std::fmt;

use crate::dump_file::DumpParsingError;

const HEADER_SYM_BOX: &str = "ITEM: BOX";
const HEADER_ATOMS: &str = "ITEM: ATOMS";

pub struct SymBox {
    pub boundaries: String,
    pub xlo: f64,
    pub xhi: f64,
    pub ylo: f64,
    pub yhi: f64,
    pub zlo: f64,
    pub zhi: f64,
}

pub struct DumpSnapshot {
    pub step: u64,
    pub atoms_count: usize,
    pub sym_box: SymBox,
    keys: HashMap<String, usize>,
    atoms: Vec<f64>,
}

impl DumpSnapshot {
    pub fn new<'a, I>(lines: I, step: u64, atoms_count: usize) -> Result<Self, DumpParsingError>
    where
        I: Iterator<Item = &'a str>,
    {
        let sym_box = match lines
            .next()
            .and_then(|l| l.split_at_checked(HEADER_SYM_BOX.len()))
        {
            Some((HEADER_SYM_BOX, boundaries)) => {
                // TODO: read box parameters
                lines.next();
                lines.next();
                lines.next();
                Ok(SymBox {
                    boundaries: boundaries.to_string(),
                    xlo: 0.0,
                    xhi: 0.0,
                    ylo: 0.0,
                    yhi: 0.0,
                    zlo: 0.0,
                    zhi: 0.0,
                })
            }
            _ => Err(DumpParsingError::MissingSymBox),
        }?;
        let mut keys = HashMap::new();
        match lines
            .next()
            .and_then(|l| l.split_at_checked(HEADER_ATOMS.len()))
        {
            Some((HEADER_ATOMS, keys)) => {
                for key in keys.split_whitespace() {
                    if keys.insert(key.to_string(), keys.len()).is_some() {
                        return Err(DumpParsingError::DuplicateAtomKeys);
                    }
                }
            }
            _ => return Err(DumpParsingError::MissingAtomKeys),
        }
        let mut snapshot = Self {
            step,
            atoms_count,
            atoms: vec![0.0; atoms_count * keys.len()],
            keys,
            sym_box,
        };
        for i in 0..atoms_count {
            let values: Vec<f64> = lines
                .next()
                .map(str::split_whitespace)
                .map(str::parse::<f64>)
                .collect();
            for (j, val) in values.iter().enumerate() {
                snapshot.atoms[atoms_count * j + i] = *val;
            }
        }
        Ok(snapshot)
    }

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

    #[inline]
    pub fn get_zero_lvl(&self) -> f64 {
        self.get_property("z")
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
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

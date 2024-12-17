use std::collections::HashMap;
use std::fmt;
use std::io;

use crate::dump_file::DumpParsingError;

pub const HEADER_TIMESTEP: &str = "ITEM: TIMESTEP";
pub const HEADER_NUM_OF_ATOMS: &str = "ITEM: NUMBER OF ATOMS";
pub const HEADER_SYM_BOX: &str = "ITEM: BOX BOUNDS";
pub const HEADER_ATOMS: &str = "ITEM: ATOMS";

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
    keys_map: HashMap<String, usize>,
    keys: Vec<String>,
    atoms: Vec<f64>,
}

impl DumpSnapshot {
    pub fn new(
        keys_map: HashMap<String, usize>,
        step: u64,
        atoms_count: usize,
        sym_box: SymBox,
    ) -> Self {
        let mut keys: Vec<(&String, &usize)> = keys_map.iter().collect();
        keys.sort_by(|a, b| a.1.cmp(b.1));
        let keys = keys.into_iter().map(|k| k.0.clone()).collect();
        Self {
            step,
            atoms_count,
            atoms: vec![0.0; atoms_count * keys_map.len()],
            keys,
            keys_map,
            sym_box,
        }
    }

    pub fn read<'a, I>(
        lines: &mut I,
        step: u64,
        atoms_count: usize,
    ) -> Result<Self, DumpParsingError>
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
                    boundaries: boundaries[1..].to_string(),
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
        let mut keys_map = HashMap::new();
        match lines
            .next()
            .and_then(|l| l.split_at_checked(HEADER_ATOMS.len()))
        {
            Some((HEADER_ATOMS, tokens)) => {
                for token in tokens.split_whitespace() {
                    if keys_map.insert(token.to_string(), keys_map.len()).is_some() {
                        return Err(DumpParsingError::DuplicateAtomKeys);
                    }
                }
            }
            _ => return Err(DumpParsingError::MissingAtomKeys),
        }
        let mut snapshot = DumpSnapshot::new(keys_map, step, atoms_count, sym_box);
        for i in 0..atoms_count {
            let values: Vec<f64> = lines
                .next()
                .map(|l| l.split_whitespace().flat_map(str::parse::<f64>).collect())
                .ok_or(DumpParsingError::InvalidOrMissingAtomRow)?;
            for (j, val) in values.iter().enumerate() {
                snapshot.atoms[atoms_count * j + i] = *val;
            }
        }
        Ok(snapshot)
    }

    pub fn write<W>(&self, w: &mut W) -> io::Result<()>
    where
        W: io::Write,
    {
        writeln!(w, "{HEADER_TIMESTEP}")?;
        writeln!(w, "{}", self.step)?;
        writeln!(w, "{HEADER_NUM_OF_ATOMS}")?;
        writeln!(w, "{}", self.atoms_count)?;
        writeln!(w, "{HEADER_SYM_BOX} {}", self.sym_box.boundaries)?;
        writeln!(w, "{} {}", self.sym_box.xlo, self.sym_box.xhi)?;
        writeln!(w, "{} {}", self.sym_box.ylo, self.sym_box.yhi)?;
        writeln!(w, "{} {}", self.sym_box.zlo, self.sym_box.zhi)?;
        writeln!(w, "{HEADER_ATOMS} {}", self.keys.join(" "))?;
        for i in 0..self.atoms_count {
            write!(w, "{}", self.atoms[i])?;
            for j in 1..self.keys.len() {
                write!(w, " {}", self.atoms[j * self.atoms_count + i])?;
            }
            writeln!(w)?;
        }
        Ok(())
    }

    #[inline]
    pub fn get_keys(&self) -> &[String] {
        &self.keys
    }

    pub fn get_property(&self, key: &str) -> &[f64] {
        let start = self.keys_map[key] * self.atoms_count;
        let end = start + self.atoms_count;
        &self.atoms[start..end]
    }

    pub fn get_property_mut(&mut self, key: &str) -> &mut [f64] {
        let start = self.keys_map[key] * self.atoms_count;
        let end = start + self.atoms_count;
        &mut self.atoms[start..end]
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

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead},
    path::Path,
};

use crate::{
    Snapshot,
    error::{Error, ErrorKind, Result},
    parser::Parser,
    snapshot::SnapshotMeta,
};

pub struct Dump {
    timesteps: HashMap<u64, usize>,
    snapshots: Vec<Snapshot>,
}

impl Dump {
    #[must_use]
    pub fn new(snapshots: Vec<Snapshot>) -> Self {
        let mut snapshots = snapshots;
        snapshots.sort_by_key(Snapshot::get_timestep);
        let timesteps = snapshots
            .iter()
            .enumerate()
            .map(|(idx, s)| (s.get_timestep(), idx))
            .collect();
        Self {
            timesteps,
            snapshots,
        }
    }

    /// Opens a `Dump` from a file path given a list of timesteps.
    ///
    /// # Errors
    ///
    /// Returns an `Error` if the dump cannot be opened or parsed.
    pub fn open<P: AsRef<Path>>(path: P, timesteps: &[u64]) -> Result<Self> {
        let mut parser = Parser::open(path)?;
        let mut snapshots = HashMap::<u64, Snapshot>::new();
        while let Some(meta) = next_snapshot_meta(&mut parser) {
            let meta = meta?;
            if !timesteps.contains(&meta.timestep) {
                continue;
            }
            if snapshots.contains_key(&meta.timestep) {
                return Err(Error::new(
                    ErrorKind::DuplicateTimestep(meta.timestep),
                    meta.timestep.to_string(),
                    parser.current(),
                ));
            }
            let snapshot = Snapshot::parse(&mut parser, meta)?;
            snapshots.insert(snapshot.get_timestep(), snapshot);
        }
        let snapshots = snapshots.into_values().collect();
        Ok(Self::new(snapshots))
    }

    /// Saves the `Dump` to a file path.
    ///
    /// # Errors
    ///
    /// Returns an `io::Error` if the dump cannot be saved.
    pub fn save(&self, path: &Path) -> io::Result<()> {
        let f = File::create(path)?;
        let mut w = io::BufWriter::new(f);
        for snapshot in self.get_snapshots() {
            snapshot.write(&mut w)?;
        }
        Ok(())
    }

    #[must_use]
    pub fn get_snapshots(&self) -> &[Snapshot] {
        &self.snapshots
    }

    #[inline]
    #[must_use]
    pub fn get_property(&self, timestep: u64, key: &str) -> &[f64] {
        let idx = self.timesteps[&timestep];
        self.snapshots[idx].get_property(key)
    }

    #[inline]
    #[must_use]
    pub fn get_property_checked(&self, timestep: u64, key: &str) -> Option<&[f64]> {
        let idx = self.timesteps.get(&timestep)?;
        let property = self.snapshots[*idx].get_property(key);
        Some(property)
    }
}

fn next_snapshot_meta<B: BufRead>(parser: &mut Parser<B>) -> Option<Result<SnapshotMeta>> {
    match SnapshotMeta::parse(parser) {
        Ok(meta) => Some(Ok(meta)),
        Err(err) => {
            if matches!(err.kind(), ErrorKind::ExpectedTimestepHeader) {
                None
            } else {
                Some(Err(err))
            }
        }
    }
}

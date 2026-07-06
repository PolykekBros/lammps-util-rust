use std::{
    collections::{HashMap, HashSet},
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

#[derive(Clone, Debug)]
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
        let filter_timesteps = !timesteps.is_empty();
        let timesteps = timesteps.iter().copied().collect::<HashSet<_>>();
        let mut parser = Parser::open(path)?;
        let mut seen = HashSet::new();
        let mut snapshots = Vec::new();
        let mut total_timesteps = 0;
        while let Some(meta) = next_snapshot_meta(&mut parser) {
            let meta = meta?;
            total_timesteps += 1;
            if filter_timesteps && !timesteps.contains(&meta.timestep) {
                parser.skip(meta.atoms_count)?;
                continue;
            }
            if seen.contains(&meta.timestep) {
                return Err(Error::new(
                    ErrorKind::DuplicateTimestep(meta.timestep),
                    meta.timestep.to_string(),
                    parser.get_current(),
                ));
            }
            let snapshot = Snapshot::parse(&mut parser, meta)?;
            seen.insert(snapshot.get_timestep());
            snapshots.push(snapshot);
        }
        if total_timesteps == 0 {
            return Err(Error::new(
                ErrorKind::EmptyFile,
                String::new(),
                parser.get_current(),
            ));
        }
        Ok(Self::new(snapshots))
    }

    /// Saves the `Dump` to a file path.
    ///
    /// # Errors
    ///
    /// Returns an `io::Error` if the dump cannot be saved.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
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
        self.snapshots[*idx].get_property_checked(key)
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

#[cfg(test)]
mod tests {
    use super::*;

    const TWO_TIMESTEP_DUMP: &str = "\
ITEM: TIMESTEP
1000
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z
1 1 0.5 0.5 0.5
ITEM: TIMESTEP
2000
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z
1 1 0.5 0.5 0.5
";

    #[test]
    fn test_open_empty_file() {
        let temp_dir = std::env::temp_dir();
        let path = temp_dir.join("empty.dump");
        std::fs::write(&path, "").unwrap();

        let result = Dump::open(&path, &[0]);
        let _ = std::fs::remove_file(&path);

        assert!(result.is_err());
        let err = result.err().unwrap();
        assert!(matches!(err.kind(), ErrorKind::EmptyFile));
        assert_eq!(err.line(), 0);
    }

    #[test]
    fn test_open_non_empty_no_timesteps() {
        let temp_dir = std::env::temp_dir();
        let path = temp_dir.join("random.dump");
        std::fs::write(&path, "random garbage text\n").unwrap();

        let result = Dump::open(&path, &[0]);
        let _ = std::fs::remove_file(&path);

        assert!(result.is_err());
        let err = result.err().unwrap();
        assert!(matches!(err.kind(), ErrorKind::EmptyFile));
        assert_eq!(err.line(), 1);
    }

    #[test]
    fn test_open_read_all_timesteps() {
        let temp_dir = std::env::temp_dir();
        let path = temp_dir.join("two_timesteps.dump");
        std::fs::write(&path, TWO_TIMESTEP_DUMP).unwrap();

        let result = Dump::open(&path, &[]);
        let _ = std::fs::remove_file(&path);

        assert!(result.is_ok());
        let dump = result.unwrap();
        assert_eq!(dump.get_snapshots().len(), 2);
        assert_eq!(dump.get_snapshots()[0].get_timestep(), 1000);
        assert_eq!(dump.get_snapshots()[1].get_timestep(), 2000);
    }

    #[test]
    fn test_open_nonexistent_timesteps() {
        let temp_dir = std::env::temp_dir();
        let path = temp_dir.join("two_timesteps_filtered.dump");
        std::fs::write(&path, TWO_TIMESTEP_DUMP).unwrap();

        let result = Dump::open(&path, &[999]);
        let _ = std::fs::remove_file(&path);

        assert!(result.is_ok());
        let dump = result.unwrap();
        assert_eq!(dump.get_snapshots().len(), 0);
    }
}

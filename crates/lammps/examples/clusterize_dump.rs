use anyhow::Result;
use lammps_util_rust::{clusterize_snapshot, DumpFile};
use std::path::Path;

fn main() -> Result<()> {
    let dump = DumpFile::read(Path::new("examples/dump.simple"), &Vec::new())?;
    let snapshot = dump.get_snapshots()[0];
    let snapshot_cluster = clusterize_snapshot(snapshot, 3.0);
    let dump_cluster = DumpFile::new(vec![snapshot_cluster]);
    dump_cluster.save(Path::new("examples/dump.simple_clusterized"))?;
    Ok(())
}

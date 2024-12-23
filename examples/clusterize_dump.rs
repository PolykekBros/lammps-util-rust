use lammps_util_rust::{Clusterizer, DumpFile};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dump = DumpFile::read(Path::new("examples/dump.simple"), &Vec::new())?;
    let dump = Clusterizer::new(dump.get_snapshots()[0], 3.0).clusterize();
    dump.save(Path::new("examples/dump.simple_clusterized"))?;
    Ok(())
}

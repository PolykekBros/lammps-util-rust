use lammps_util_rust::DumpFile;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dump = DumpFile::read(Path::new("examples/dump.simple"), &Vec::new())?;
    dump.save(Path::new("examples/dump.simple_saved"))?;
    Ok(())
}

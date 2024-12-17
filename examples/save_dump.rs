use lammps_util_rust::dump_file::DumpFile;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dump = DumpFile::new(Path::new("examples/dump.simple"), &Vec::new())?;
    dump.save(Path::new("examples/dump.simple_saved"))?;
    Ok(())
}

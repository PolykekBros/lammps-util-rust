use lammps_util_rust::{DumpFile, DumpParsingError};
use std::path::Path;

fn main() -> Result<(), DumpParsingError> {
    let dump = DumpFile::new(Path::new("examples/dump.simple"), &Vec::new())?;
    for snapshot in dump.get_snapshots() {
        println!("\ntimestep: {}", snapshot.step);
        let keys = snapshot.get_keys();
        print!("keys:");
        for key in keys.iter() {
            print!(" {key}")
        }
        println!();
        for i in 0..snapshot.atoms_count {
            for key in keys.iter() {
                print!("{}\t", snapshot.get_property(key)[i]);
            }
            println!();
        }
    }
    Ok(())
}

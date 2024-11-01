use lammps_util_rust::{Dump, DumpParsingError};
use std::path::Path;

fn main() -> Result<(), DumpParsingError> {
    let dump = Dump::new(Path::new("examples/dump.simple"), &Vec::new())?;
    for tstep in dump.timesteps.keys() {
        println!("\ntimestep: {tstep}");
        print!("keys:");
        for key in dump.get_keys() {
            print!(" {key}")
        }
        println!();
    }
    Ok(())
}

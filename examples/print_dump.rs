use lammps_util_rust::{Dump, DumpParsingError};
use std::path::Path;

fn main() -> Result<(), DumpParsingError> {
    let dump = Dump::new(Path::new("examples/dump.simple"), &Vec::new())?;
    let keys = dump.get_keys();
    for tstep in dump.timesteps.keys() {
        println!("\ntimestep: {tstep}");
        print!("keys:");
        for key in keys.iter() {
            print!(" {key}")
        }
        println!();
        let mut table = Vec::new();
        for key in keys.iter() {
            table.push(dump.get_property(*tstep, key));
        }
        for i in 0..dump.timesteps[tstep].atom_count {
            for (pos, _) in keys.iter().enumerate() {
                print!("{}\t", table[pos][i]);
            }
            println!();
        }
    }
    Ok(())
}

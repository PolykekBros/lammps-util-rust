use lammps_util::DumpFile;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dump = DumpFile::read(Path::new("examples/dump.simple"), &Vec::new())?;
    for snapshot in dump.get_snapshots() {
        println!("\ntimestep: {}", snapshot.get_timestep());
        let keys = snapshot.get_keys();
        print!("keys:");
        for key in keys.iter() {
            print!(" {key}")
        }
        println!();
        for i in 0..snapshot.get_atoms_count() {
            for key in keys.iter() {
                print!("{}\t", snapshot.get_property(key)[i]);
            }
            println!();
        }
    }
    Ok(())
}

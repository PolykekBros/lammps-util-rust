use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use lammps_util_rust::DumpFile;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// input dump file
    dump_input_file: PathBuf,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let dump = DumpFile::read(&cli.dump_input_file, &[])?;
    let zero_lvl = dump.get_snapshots()[0].get_zero_lvl();
    println!("zero_lvl: {zero_lvl}");
    Ok(())
}

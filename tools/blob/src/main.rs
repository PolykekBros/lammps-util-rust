use std::{
    iter,
    path::{Path, PathBuf},
};

use anyhow::Result;
use clap::Parser;
use lammps_util_rust::{DumpFile, IteratorAvg};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    dump_path: PathBuf,
}

fn parse_dump(p: &Path) -> Result<()> {
    let dump = DumpFile::read(p, &[3000])?;
    let snapshot = dump.get_snapshots()[0];
    let ax = snapshot.get_property("x");
    let ay = snapshot.get_property("y");
    let ek = snapshot.get_property("c_atom_ke");
    let ek_avg = ek.iter().copied().avg().unwrap();
    let ax_avg = ax.iter().copied().avg().unwrap();
    let ay_avg = ay.iter().copied().avg().unwrap();
    println!("# <Ek>: {ek_avg}, <x>: {ax_avg}, <y>: {ay_avg}");
    println!("# | R | Ek |");
    iter::zip(ax, ay)
        .map(|(ax, ay)| {
            let x = ax - ax_avg;
            let y = ay - ay_avg;
            (x.powi(2) + y.powi(2)).sqrt()
        })
        .zip(ek)
        .map(|(r, ek)| (r, ek - ek_avg))
        .for_each(|(r, ek)| println!("{r}\t{ek}"));
    Ok(())
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    parse_dump(&cli.dump_path)?;
    Ok(())
}

use anyhow::Result;
use clap::Parser;
use lammps_util_rust::{crater_snapshot, DumpFile, DumpSnapshot};
use std::path::PathBuf;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// DUMP INPUT
    dump_input_file: PathBuf,

    /// DUMP FINAL
    dump_final_file: PathBuf,

    /// OUTPUT DIR
    output_dir: PathBuf,

    /// Max depth (A)
    #[arg(short, long, default_value_t = 50.0)]
    max_depth: f64,

    /// Cutoff (A)
    #[arg(short, long, default_value_t = 3.0)]
    cutoff: f64,
}

fn get_crater_info(snapshot: &DumpSnapshot, zero_lvl: f64) -> String {
    let z = snapshot.get_property("z");
    let mut crater_count = 0;
    let mut surface_count = 0;
    let mut z_avg = 0.0;
    let mut z_min = f64::INFINITY;
    for &z in z {
        if z > -2.4 * 0.707 + zero_lvl {
            surface_count += 1;
        }
        crater_count += 1;
        z_min = z_min.min(z - zero_lvl);
        z_avg += z - zero_lvl;
    }
    z_avg /= crater_count as f64;
    let volume = crater_count as f64 * 20.1;
    let surface = surface_count as f64 * 7.3712;
    format!("{crater_count} {volume} {surface} {z_avg} {z_min}")
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let dump_input = DumpFile::read(&cli.dump_input_file, &[])?;
    let snapshot_input = dump_input.get_snapshots()[0];
    let zero_lvl = snapshot_input.get_zero_lvl();

    let dump_final = DumpFile::read(&cli.dump_final_file, &[])?;
    let snapshot_final = dump_final.get_snapshots()[0];

    let snapshot_crater = crater_snapshot(snapshot_input, snapshot_final, 3.0);
    println!("crater atoms: {}", snapshot_crater.atoms_count);
    let info = get_crater_info(&snapshot_crater, zero_lvl);

    let dump_crater = DumpFile::new(vec![snapshot_crater]);
    dump_crater.save(&cli.output_dir.join("dump.crater"))?;
    println!("{info}");
    Ok(())
}

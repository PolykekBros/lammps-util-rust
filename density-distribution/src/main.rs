use clap::Parser;
use lammps_util_rust::{DumpFile, DumpSnapshot};
use plotters::prelude::LogScalable;
use std::{ops::Deref, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[arg(value_name = "DUMP_FILE")]
    dump_file: PathBuf,

    #[arg(short, long, value_name = "TIMESTEP")]
    timestep: Option<u64>,

    #[arg(short, long, value_name = "delta")]
    delta: Option<f64>,
}

fn get_distribution(dump: &DumpSnapshot, delta: f64) -> (Vec<f64>, Vec<Vec<f64>>) {
    let atom_id = dump.get_property("id");
    let atom_x = dump.get_property("x");
    let atom_y = dump.get_property("y");
    let atom_z = dump.get_property("z");
    let mut ids = Vec::new();
    for &id in atom_id {
        let id_usize = id as usize;
        if !ids.contains(&id_usize) {
            ids.push(id_usize.to_owned());
        }
    }
    let x_min = atom_x.iter().fold(f64::INFINITY, |acc, x| acc.min(*x));
    let x_max = atom_x.iter().fold(f64::NEG_INFINITY, |acc, x| acc.max(*x));
    let y_min = atom_y.iter().fold(f64::INFINITY, |acc, y| acc.min(*y));
    let y_max = atom_y.iter().fold(f64::NEG_INFINITY, |acc, y| acc.max(*y));
    let z_min = atom_z.iter().fold(f64::INFINITY, |acc, z| acc.min(*z));
    let z_max = atom_z.iter().fold(f64::NEG_INFINITY, |acc, z| acc.max(*z));
    let volume = (z_max - z_min) * (y_max - y_min) * (x_max - x_min);
    let count = (z_max - z_min).ceil() as usize;
    let mut plot_x = Vec::with_capacity(count);
    let mut plot_y = Vec::with_capacity(ids.len());
    for _ in &ids {
        plot_y.push(Vec::<f64>::with_capacity(count));
    }
    for i in 0..count {
        let i_f64 = i as f64;
        plot_x.push(i_f64 * delta + delta / 2.0f64);
        for &id in &ids {
            let atom_id_z = std::iter::zip(atom_id, atom_z);
            plot_y[id - 1].push(
                atom_id_z
                    .filter(|(&fid, &z)| {
                        (z >= i_f64 * delta && z < (i_f64 + 1.0f64) * delta) && id == (fid as usize)
                    })
                    .count() as f64
                    / volume,
            );
        }
    }
    (plot_x, plot_y)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    let dump_path = cli.dump_file;
    let timesteps = match cli.timestep {
        Some(timestep) => vec![timestep],
        _ => Vec::new(),
    };
    let dump = DumpFile::new(dump_path.as_path(), &timesteps)?;
    let (plot_x, plot_y) =
        get_distribution(dump.get_snapshots()[0], cli.delta.unwrap_or(5.43 * 2.0));
    Ok(())
}

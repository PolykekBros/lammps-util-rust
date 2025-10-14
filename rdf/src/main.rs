use lammps_util_rust::{DumpFile, DumpSnapshot, XYZ};
use anyhow::Result;
use clap::Parser;
use std::path::{Path, PathBuf};
use log::info;


#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    dump_file: PathBuf,

    #[arg(short, long)]
    timestep: Option<u64>,

    #[arg(short, long)]
    cutoff: f64,

    #[arg(short, long)]
    binsize: f64,
    
    #[arg(short, long)]
    output_dir: PathBuf,
}

fn get_bins(cutoff: f64, binsize: f64) -> Vec<f64> {
    let bin_num = (cutoff / binsize).floor() as usize;
    (1..=bin_num).map(| bin | (bin as f64 - 0.5)*binsize).collect()
}

fn get_rdf(cutoff: f64, binsize: f64, dump: &DumpSnapshot) -> Vec<(f64, usize)> {
    let bins = get_bins(cutoff, binsize);
    let mut rdf = bins.into_iter().map(| bin | (bin, 0)).collect::<Vec<(f64, usize)>>();
    let coords = dump.get_coordinates();
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords.clone());
    for atom in coords {
        rdf.iter_mut().for_each(| (r, n) | {
            if *r > 0.0 {
                let neigh = kdtree.within_radius(&atom, *r + 0.5*binsize).len();
                *n += neigh - kdtree.within_radius(&atom, *r - 0.5*binsize).len();
            }
        });
    };
    rdf
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let dump_path = cli.dump_file;
    let timesteps = match cli.timestep {
        Some(timestep) => vec![timestep],
        _ => vec![],
    };
    let dump = DumpFile::read(dump_path.as_path(), &timesteps)?;
    let snapshot = dump.get_snapshots()[0];
    Ok(())
}

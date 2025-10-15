use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use lammps_util_rust::{DumpFile, DumpSnapshot};
use rayon::prelude::*;
use std::{f64, iter, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    dump_file: PathBuf,

    #[arg(short, long)]
    timestep: Option<u64>,

    #[arg(short, long)]
    cutoff: f64,

    #[arg(short, long)]
    n_bins: usize,
}

fn get_bins(cutoff: f64, n: usize) -> Vec<(f64, f64)> {
    let delta = cutoff / n as f64;
    (0..n)
        .map(|i| i as f64)
        .map(|i| (i * delta, (i + 1.0) * delta))
        .collect()
}

fn normalize(
    rdf: impl IntoIterator<Item = ((f64, f64), f64)>,
    dump: &DumpSnapshot,
) -> Vec<(f64, f64)> {
    let vol = dump.sym_box.volume();
    let num = dump.atoms_count as f64;
    let rho = num / vol;
    rdf.into_iter()
        .map(|((lo, hi), mut n)| {
            let vshell = 4.0 / 3.0 * f64::consts::PI * (hi.powi(3) - lo.powi(3));
            let nnorm = rho * vshell;
            n /= nnorm * num;
            ((hi + lo) / 2.0, n)
        })
        .collect()
}

fn get_rdf(cutoff: f64, n: usize, dump: &DumpSnapshot) -> Vec<(f64, f64)> {
    let bins = get_bins(cutoff, n);
    let coords = dump.get_coordinates();
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords.clone());
    // TODO: acoount for periodic boundaries
    let rdf = coords
        .par_iter()
        .map(|atom| {
            bins.iter()
                .copied()
                .map(|(lo, hi)| {
                    let n_lo = kdtree.within_radius(atom, lo).len();
                    let n_hi = kdtree.within_radius(atom, hi).len();
                    let n = (n_hi - n_lo) as f64;
                    ((lo, hi), n)
                })
                .collect::<Vec<_>>()
        })
        .reduce(
            || {
                bins.iter()
                    .map(|&(lo, hi)| ((lo, hi), 0.0))
                    .collect::<Vec<_>>()
            },
            |mut a, b| {
                iter::zip(&mut a, b).for_each(|(a, b)| a.1 += b.1);
                a
            },
        );
    normalize(rdf, dump)
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
    let rdf = get_rdf(cli.cutoff, cli.n_bins, snapshot);
    let table = rdf
        .into_iter()
        .map(|(r, n)| [r, n].into_iter().map(|x| format!("{x:10.4}")).join("\t"))
        .join("\n");
    println!("# radius number\n{table}");
    Ok(())
}

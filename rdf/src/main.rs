use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use lammps_util_rust::{DumpFile, DumpSnapshot, XYZ};
use rayon::prelude::*;
use std::{iter, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    dump_file: PathBuf,

    #[arg(short, long)]
    timestep: Option<u64>,

    #[arg(short, long)]
    cutoff: f32,

    #[arg(short, long)]
    n_bins: usize,
}

fn get_bins(cutoff: f32, n: usize) -> Vec<(f32, f32)> {
    let delta = cutoff / n as f32;
    (0..n)
        .map(|i| i as f32)
        .map(|i| (i * delta, (i + 1.0) * delta))
        .collect()
}

fn normalize(
    rdf: impl IntoIterator<Item = ((f32, f32), usize)>,
    dump: &DumpSnapshot,
) -> Vec<(f32, f32)> {
    let vol = dump.sym_box.volume();
    let num = dump.atoms_count as f32;
    let rho = num / vol;
    rdf.into_iter()
        .map(|((lo, hi), n)| {
            let vshell = 4.0 / 3.0 * std::f32::consts::PI * (hi.powi(3) - lo.powi(3));
            let nnorm = rho * vshell;
            let n = n as f32 / (nnorm * num);
            (hi.midpoint(lo), n)
        })
        .collect()
}

fn get_rdf(cutoff: f32, n: usize, dump: &DumpSnapshot) -> Vec<(f32, f32)> {
    let bins = get_bins(cutoff, n);
    let mut coords = dump.get_coordinates();
    // XYZ::get_supercell_coords(&mut coords, &dump.sym_box, cutoff);
    println!("before tree");
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords);
    let rdf = kdtree
        .items()
        .par_iter()
        .filter(|atom| atom.index < dump.atoms_count)
        .map(|atom| {
            let d_sq = kdtree
                .within_radius(atom, cutoff)
                .into_iter()
                .map(|neigh| atom.distance_squared(neigh.coords))
                .collect::<Vec<_>>();
            bins.iter()
                .map(|&(lo, hi)| {
                    let lo_sq = lo.powi(2);
                    let hi_sq = hi.powi(2);
                    let n = d_sq
                        .iter()
                        .filter(|&&d_sq| d_sq >= lo_sq && d_sq < hi_sq && d_sq != 0.0)
                        .count();
                    ((lo, hi), n)
                })
                .collect::<Vec<_>>()
        })
        .reduce(
            || {
                bins.iter()
                    .map(|&(lo, hi)| ((lo, hi), 0))
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
    let timesteps = cli.timestep.map(|t| vec![t]).unwrap_or_default();
    println!("before read dump");
    let dump = DumpFile::read(dump_path.as_path(), &timesteps)?;
    println!("read dump");
    let snapshot = dump.get_snapshots()[0];
    let rdf = get_rdf(cli.cutoff, cli.n_bins, snapshot);
    let table = rdf
        .into_iter()
        .map(|vals| {
            <[_; 2]>::from(vals)
                .into_iter()
                .map(|x| format!("{x:10.4}"))
                .join("\t")
        })
        .join("\n");
    println!("# radius number\n{table}");
    Ok(())
}

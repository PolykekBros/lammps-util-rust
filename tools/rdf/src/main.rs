#![allow(clippy::cast_precision_loss)]
use clap::Parser;
use lammps_files::{Dump, Snapshot};
use lammps_util::{XYZ, xyz::xyz_vec_from_snapshot};
use rayon::prelude::*;

use std::{fmt, iter, path::PathBuf};

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
    rdf: impl IntoIterator<Item = ((f64, f64), usize)>,
    snapshot: &Snapshot,
) -> Vec<(f64, f64)> {
    let vol = snapshot.get_symbox().bbox.volume();
    let num = snapshot.get_atoms_count() as f64;
    let rho = num / vol;
    rdf.into_iter()
        .map(|((lo, hi), n)| {
            let vshell = 4.0 / 3.0 * std::f64::consts::PI * (hi.powi(3) - lo.powi(3));
            let nnorm = rho * vshell;
            let n = n as f64 / (nnorm * num);
            (hi.midpoint(lo), n)
        })
        .collect()
}

fn get_rdf(cutoff: f64, n: usize, snapshot: &Snapshot) -> Vec<(f64, f64)> {
    let bins = get_bins(cutoff, n);
    let mut coords = xyz_vec_from_snapshot(snapshot);
    XYZ::get_supercell_coords(&mut coords, snapshot.get_symbox(), cutoff);
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords);
    let rdf = kdtree
        .items()
        .par_iter()
        .filter(|atom| !atom.is_ghost)
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
    normalize(rdf, snapshot)
}

fn print_table<T: fmt::Display>(table: &[T], header: &[&str], rows: usize, cols: usize) {
    if cols == 0 || rows == 0 {
        return;
    }
    assert_eq!(header.len(), cols);
    assert_eq!(table.len(), cols * rows);
    let get_idx = |row_idx: usize, col_idx: usize| row_idx * cols + col_idx;
    let mut widths = header.iter().map(|h| h.chars().count()).collect::<Vec<_>>();
    for row_idx in 0..rows {
        let start = get_idx(row_idx, 0);
        let end = get_idx(row_idx + 1, 0);
        for (w, t) in iter::zip(&mut widths, &table[start..end]) {
            *w = format!("{t:.06}").len().max(*w);
        }
    }
    print!("#");
    for (h, width) in iter::zip(header, &widths) {
        print!("{h:>width$}", width = width + 1);
    }
    println!();
    for row_idx in 0..rows {
        print!(" ");
        for (col_idx, width) in iter::zip(0..cols, &widths) {
            print!(
                "{:>width$.06}",
                table[get_idx(row_idx, col_idx)],
                width = width + 1
            );
        }
        println!();
    }
}

fn main() {
    env_logger::init();
    let cli = Cli::parse();
    let dump_path = cli.dump_file;
    let timesteps = cli.timestep.map(|t| vec![t]).unwrap_or_default();
    let dump = Dump::open(dump_path, &timesteps).unwrap();
    let snapshot = &dump.get_snapshots()[0];
    let rdf = get_rdf(cli.cutoff, cli.n_bins, snapshot);
    let rows = rdf.len();
    let cols = 2;
    let header = ["r", "g(r)"];
    let table = rdf
        .into_iter()
        .flat_map(|(x, y)| [x, y])
        .collect::<Vec<_>>();
    print_table(&table, &header, rows, cols);
}

#![allow(clippy::cast_precision_loss)]
use clap::Parser;
use kd_tree::KdTree3;
use lammps_files::{Dump, Snapshot};
use lammps_util::{XYZ, xyz::xyz_vec_from_snapshot};
// use rayon::prelude::*;

use std::{fmt, iter, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    dump_file: PathBuf,

    /// types to compute rdf for - '1,2' -> 1-1 + 1-2 + 2-2 (all if none)
    #[arg(short, long, value_delimiter = ',')]
    types: Vec<usize>,

    #[arg(short, long)]
    timestep: Option<u64>,

    #[arg(short, long)]
    cutoff: f64,

    #[arg(short, long)]
    n_bins: usize,
}

fn get_bins(cutoff: f64, n: usize) -> Vec<f64> {
    let delta = cutoff / n as f64;
    (0..n)
        .map(|i| i as f64)
        .flat_map(|i| [i * delta, (i + 1.0) * delta])
        .collect()
}

fn make_table<'a>(
    bins: impl IntoIterator<Item = &'a [f64; 2]>,
    rows: usize,
    cols: usize,
) -> Vec<f64> {
    let mut table = vec![0.0; rows * cols];
    for (idx, &[lo, hi]) in bins.into_iter().enumerate() {
        table[idx * cols] = hi.midpoint(lo);
    }
    table
}

fn normalize_g<'a>(
    bins: impl IntoIterator<Item = &'a [f64; 2]>,
    g: impl IntoIterator<Item = usize>,
    snapshot: &Snapshot,
) -> Vec<f64> {
    let vol = snapshot.get_symbox().bbox.volume();
    let num = snapshot.get_atoms_count() as f64;
    let rho = num / vol;
    iter::zip(bins, g)
        .map(|([lo, hi], n)| {
            let vshell = 4.0 / 3.0 * std::f64::consts::PI * (hi.powi(3) - lo.powi(3));
            let nnorm = rho * vshell;
            n as f64 / (nnorm * num)
        })
        .collect()
}

fn get_rdf_new(snapshot: &Snapshot, n: usize, cutoff: f64, types: &[usize]) -> Vec<f64> {
    let bins_vec = get_bins(cutoff, n);
    let (bins, rem) = bins_vec.as_chunks::<2>();
    assert_eq!(rem.len(), 0);
    assert_eq!(bins.len(), n);
    assert_eq!(bins.len(), n * 2);
    let mut coords = xyz_vec_from_snapshot(snapshot);
    XYZ::get_supercell_coords(&mut coords, snapshot.get_symbox(), cutoff);
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords);
    match types.len() {
        0 => {
            let mut table = make_table(bins, bins.len(), 2);
            let g = calculate_rdf(snapshot, &kdtree, bins, cutoff, |_| true, |_| true);
            assert_eq!(g.len(), bins.len());
            for (idx, g) in g.into_iter().enumerate() {
                table[idx * 2 + 1] = g;
            }
            table
        }
        _ => {
            let cols = (0..types.len())
                .map(|i| (i..types.len()).count())
                .sum::<usize>();
            let mut table = make_table(bins, bins.len(), cols);
            let atypes = snapshot.get_property("type");
            for (col, (ti, tj)) in types
                .iter()
                .enumerate()
                .flat_map(|(i, ti)| types[i..].iter().map(move |tj| (ti, tj)))
                .enumerate()
            {
                let g = calculate_rdf(
                    snapshot,
                    &kdtree,
                    bins,
                    cutoff,
                    |atom| (atypes[atom.index] as usize).eq(ti),
                    |atom| (atypes[atom.index] as usize).eq(tj),
                );
                assert_eq!(g.len(), bins.len());
                for (row, g) in g.into_iter().enumerate() {
                    table[row * cols + col + 1] = g;
                }
            }
            todo!()
        }
    }
}

fn calculate_rdf<F1, F2>(
    snapshot: &Snapshot,
    kdtree: &KdTree3<XYZ>,
    bins: &[[f64; 2]],
    cutoff: f64,
    atom_1_filter: F1,
    atom_2_filter: F2,
) -> Vec<f64>
where
    F1: Fn(&&XYZ) -> bool,
    F2: Fn(&&XYZ) -> bool,
{
    let g = kdtree
        .items()
        .iter()
        .filter(|atom| !atom.is_ghost)
        .filter(atom_1_filter)
        .map(|atom| calculate_rdf_hist(kdtree, bins, cutoff, atom, &atom_2_filter))
        .fold(Vec::<usize>::new(), |mut g, part_g| {
            iter::zip(&mut g, part_g).for_each(|(g, p_g)| *g += p_g);
            g
        });
    normalize_g(bins, g, snapshot)
}

fn calculate_rdf_hist<F>(
    kdtree: &KdTree3<XYZ>,
    bins: &[[f64; 2]],
    cutoff: f64,
    atom: &XYZ,
    atom_2_filter: F,
) -> Vec<usize>
where
    F: Fn(&&XYZ) -> bool,
{
    let d_sq = kdtree
        .within_radius(atom, cutoff)
        .into_iter()
        .filter(atom_2_filter)
        .map(|neigh| atom.distance_squared(neigh.coords))
        .collect::<Vec<_>>();
    bins.iter()
        .map(|&[lo, hi]| {
            let lo_sq = lo.powi(2);
            let hi_sq = hi.powi(2);
            d_sq.iter()
                .filter(|&&d_sq| d_sq >= lo_sq && d_sq < hi_sq && d_sq != 0.0)
                .count()
        })
        .collect::<Vec<_>>()
}

fn print_table<T, H>(table: &[T], header: &[H], rows: usize, cols: usize)
where
    T: fmt::Display,
    H: fmt::Display,
{
    if cols == 0 || rows == 0 {
        return;
    }
    assert_eq!(header.len(), cols);
    assert_eq!(table.len(), cols * rows);
    let get_idx = |row_idx: usize, col_idx: usize| row_idx * cols + col_idx;
    let mut widths = header
        .iter()
        .map(|h| format!("{h}").len())
        .collect::<Vec<_>>();
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

fn make_rdf_table_header(types: &[usize]) -> Vec<String> {
    match types.len() {
        0 => vec!["r".to_string(), "g(r)".to_string()],
        _ => {
            let mut v = Vec::new();
            v.push("r".to_string());
            for s in types
                .iter()
                .enumerate()
                .flat_map(|(i, ti)| types[i..].iter().map(move |tj| (ti, tj)))
                .map(|(ti, tj)| format!("g({ti}-{tj})(r)"))
            {
                v.push(s);
            }

            v
        }
    }
}

fn main() {
    env_logger::init();
    let cli = Cli::parse();
    let dump_path = cli.dump_file;
    let timesteps = cli.timestep.map(|t| vec![t]).unwrap_or_default();
    let dump = Dump::open(dump_path, &timesteps).unwrap();
    let snapshot = &dump.get_snapshots()[0];
    let types = cli.types.as_slice();
    let table = get_rdf_new(snapshot, cli.n_bins, cli.cutoff, types);
    let header = make_rdf_table_header(types);
    print_table(
        &table,
        header.as_slice(),
        table.len() / header.len(),
        header.len(),
    );
}

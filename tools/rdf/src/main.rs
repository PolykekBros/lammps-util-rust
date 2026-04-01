#![deny(
    clippy::all,
    clippy::correctness,
    clippy::suspicious,
    clippy::complexity,
    clippy::perf,
    clippy::style,
    clippy::pedantic
)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_possible_truncation)]
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
    #[arg(short = 'T', long, value_delimiter = ',')]
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

fn make_table(bins: &[[f64; 2]], cols: usize) -> Vec<f64> {
    let mut table = vec![0.0; bins.len() * cols];
    for (idx, &[lo, hi]) in bins.iter().enumerate() {
        table[idx * cols] = hi.midpoint(lo);
    }
    table
}

fn normalize_g<'a>(
    g: &mut [f64],
    bins: impl IntoIterator<Item = &'a [f64; 2]>,
    snapshot: &Snapshot,
) {
    let vol = snapshot.get_symbox().bbox.volume();
    let num = snapshot.get_atoms_count() as f64;
    let rho = num / vol;
    iter::zip(g, bins).for_each(|(n, [lo, hi])| {
        let vshell = 4.0 / 3.0 * std::f64::consts::PI * (hi.powi(3) - lo.powi(3));
        let nnorm = rho * vshell;
        *n /= nnorm * num;
    });
}

fn get_rdf_new(
    snapshot: &Snapshot,
    n: usize,
    cutoff: f64,
    type_pairs: &[(usize, usize)],
) -> Vec<f64> {
    let bins_vec = get_bins(cutoff, n);
    assert_eq!(bins_vec.len(), n * 2);
    let (bins, rem) = bins_vec.as_chunks::<2>();
    assert_eq!(rem.len(), 0);
    assert_eq!(bins.len(), n);
    let mut coords = xyz_vec_from_snapshot(snapshot);
    XYZ::get_supercell_coords(&mut coords, snapshot.get_symbox(), cutoff);
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords);
    if type_pairs.is_empty() {
        let cols = 1 + 1;
        let mut table = make_table(bins, cols);
        let g = calculate_rdf(snapshot, &kdtree, bins, cutoff, |_| true, |_| true);
        assert_eq!(g.len(), bins.len());
        for (idx, g) in g.into_iter().enumerate() {
            table[idx * 2 + 1] = g;
        }
        table
    } else {
        let cols = type_pairs.len() + 1;
        let mut table = make_table(bins, cols);
        let atypes = snapshot.get_property("type");
        for (col, (ti, tj)) in type_pairs.iter().enumerate() {
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
        table
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
    let mut g = kdtree
        .items()
        .iter()
        .filter(|atom| !atom.is_ghost)
        .filter(atom_1_filter)
        .map(|atom| calculate_rdf_hist(kdtree, bins, cutoff, atom, &atom_2_filter))
        .fold(vec![0.0; bins.len()], |mut g, part_g| {
            iter::zip(&mut g, part_g).for_each(|(g, p_g)| *g += p_g);
            g
        });
    let ssum = g.iter().copied().sum::<f64>() as usize;
    println!("sum {ssum}");
    normalize_g(&mut g, bins, snapshot);
    g
}

fn calculate_rdf_hist<F>(
    kdtree: &KdTree3<XYZ>,
    bins: &[[f64; 2]],
    cutoff: f64,
    atom: &XYZ,
    atom_2_filter: F,
) -> Vec<f64>
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
                .count() as f64
        })
        .collect()
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

fn make_rdf_table_header(type_pairs: &[(usize, usize)]) -> Vec<String> {
    let mut v = vec!["r".to_string()];
    if type_pairs.is_empty() {
        v.push("g(r)".to_string());
    } else {
        v.extend(type_pairs.iter().map(|(ti, tj)| format!("g({ti}-{tj})(r)")));
    }
    v
}

fn get_type_pairs(mut types: Vec<usize>) -> Vec<(usize, usize)> {
    if types.is_empty() {
        return Vec::new();
    }
    types.sort_unstable();
    let mut uniq = vec![types[0]];
    let mut last = types[0];
    for &t in &types[1..] {
        if t == last {
            continue;
        }
        uniq.push(t);
        last = t;
    }
    let types = &uniq;
    let n = uniq.len();
    (0..n)
        .flat_map(|i| (i..n).map(move |j| (types[i], types[j])))
        .collect()
}

fn main() {
    env_logger::init();
    let Cli {
        dump_file,
        types,
        timestep,
        cutoff,
        n_bins,
    } = Cli::parse();
    let type_pairs = get_type_pairs(types);
    let timesteps = timestep.map(|t| vec![t]).unwrap_or_default();
    let dump = Dump::open(dump_file, &timesteps).unwrap();
    let snapshot = &dump.get_snapshots()[0];
    let table = get_rdf_new(snapshot, n_bins, cutoff, &type_pairs);
    let header = make_rdf_table_header(&type_pairs);
    print_table(
        &table,
        header.as_slice(),
        table.len() / header.len(),
        header.len(),
    );
}

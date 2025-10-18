use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use lammps_util_rust::{DumpFile, DumpSnapshot, SymBox, XYZ};
// use rayon::prelude::*;
use std::{array, f64, iter, ops::Deref, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    dump_file: PathBuf,

    #[arg(short, long)]
    timestep: Option<u64>,

    #[arg(short, long)]
    n_bins: usize,

    #[arg(short, long)]
    type_i: usize,

    #[arg(short, long)]
    type_j: usize,

    #[arg(short, long)]
    type_k: usize,

    #[arg(short, long)]
    cutoff_i: f64,

    #[arg(short, long)]
    cutoff_j: f64,
}

fn get_bins(n: usize) -> Vec<(f64, f64)> {
    let delta = 180.0 / n as f64;
    (0..n)
        .map(|i| i as f64)
        .map(|i| (i * delta, (i + 1.0) * delta))
        .collect()
}

fn get_cos((atom_i, atom_j, atom_k): (&XYZ, &XYZ, &XYZ)) -> f64 {
    let (atom_i, atom_j, atom_k) = (*atom_i.clone(), *atom_j.clone(), *atom_k.clone());
    let r1 = atom_i - atom_k;
    let r2 = atom_j - atom_k;
    0.0
}

fn get_adf(
    type_i: usize,
    type_j: usize,
    type_k: usize,
    cutoff_i: f32,
    cutoff_j: f32,
    n: usize,
    dump: &DumpSnapshot,
) -> Vec<(f64, f64)> {
    let bins = get_bins(n);
    let mut coords = dump.get_coordinates();
    let d_types = dump.get_property("type");
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords);
    kdtree
        .items()
        .iter()
        .filter(|atom| d_types[atom.index()] as usize == type_k)
        .flat_map(|atom_k| {
            let i_neigh = kdtree.within_radius(atom_k, cutoff_i);
            let mut j_neigh = kdtree.within_radius(atom_k, cutoff_j);
            j_neigh.retain(|atom_j| {
                d_types[atom_j.index()] as usize == type_j && atom_j.index() != atom_k.index()
            });
            i_neigh
                .into_iter()
                .filter(|atom_i| {
                    d_types[atom_i.index()] as usize == type_i && atom_i.index() != atom_k.index()
                })
                .flat_map(move |atom_i| j_neigh.iter().map(move |atom_j| (atom_i, *atom_j, atom_k)))
        })
        .map(get_cos);
    Vec::new()
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
    let adf = get_adf(
        cli.type_i,
        cli.type_j,
        cli.type_k,
        cli.cutoff_i as f32,
        cli.cutoff_j as f32,
        cli.n_bins,
        snapshot,
    );
    Ok(())
}

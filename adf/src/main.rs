use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use lammps_util_rust::{DumpFile, DumpSnapshot, XYZ};
use std::{f64, iter, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    dump_file: PathBuf,

    #[arg(short, long)]
    timestep: Option<u64>,

    #[arg(short, long)]
    n_bins: usize,

    #[arg(long)]
    type_i: usize,

    #[arg(long)]
    type_j: usize,

    #[arg(long)]
    type_k: usize,

    #[arg(long)]
    cutoff_i: f64,

    #[arg(long)]
    cutoff_j: f64,
}

#[allow(clippy::cast_precision_loss)]
fn get_bins(n: usize) -> Vec<(f64, f64)> {
    let delta = 180.0 / n as f64;
    (0..n)
        .map(|i| i as f64)
        .map(|i| (i * delta, (i + 1.0) * delta))
        .collect()
}

fn get_cos(atom_i: XYZ, atom_j: XYZ, atom_k: XYZ) -> f64 {
    let r1 = *atom_i - *atom_k;
    let r2 = *atom_j - *atom_k;
    (r1.dot(r2) / (r1.length() * r2.length())).clamp(-1.0, 1.0)
}

#[allow(clippy::cast_precision_loss)]
fn normalize(adf: Vec<((f64, f64), usize)>, n: usize) -> Vec<(f64, f64, usize)> {
    let int: f64 = adf
        .iter()
        .map(|(bonds, val)| {
            let theta = bonds.0.midpoint(bonds.1);
            let d_theta = bonds.1 - bonds.0;
            *val as f64 * theta.to_radians().sin() * d_theta
        })
        .sum();
    println!("{int}");
    adf.into_iter()
        .map(|((lo, hi), val)| (lo.midpoint(hi), val as f64 / int, val / n))
        .collect()
}

#[allow(clippy::cast_possible_truncation)]
#[allow(clippy::cast_sign_loss)]
fn get_adf(
    type_i: usize,
    type_j: usize,
    type_k: usize,
    cutoff_i: f64,
    cutoff_j: f64,
    n: usize,
    dump: &DumpSnapshot,
) -> Vec<(f64, f64, usize)> {
    let bins = get_bins(n);
    let coords = dump.get_coordinates();
    let d_types = dump.get_property("type");
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords);
    let adf = kdtree
        .items()
        .iter()
        .filter(|atom| d_types[atom.index] as usize == type_k)
        .map(|atom_k| {
            let i_neigh = kdtree.within_radius(atom_k, cutoff_i);
            let j_neigh = {
                let mut j_neigh = kdtree.within_radius(atom_k, cutoff_j);
                j_neigh.retain(|atom_j| {
                    d_types[atom_j.index] as usize == type_j && atom_j.index != atom_k.index
                });
                j_neigh
            };
            let j_neigh_ref = &j_neigh;
            let angles = i_neigh
                .iter()
                .filter(|atom_i| {
                    d_types[atom_i.index] as usize == type_i && atom_i.index != atom_k.index
                })
                .flat_map(move |atom_i| {
                    j_neigh_ref
                        .iter()
                        .filter(|atom_j| atom_j.index != atom_i.index)
                        .copied()
                        .map(move |atom_j| get_cos(**atom_i, *atom_j, *atom_k).acos().to_degrees())
                })
                .collect::<Vec<_>>();
            bins.iter()
                .map(|&(lo, hi)| {
                    let mut n = angles
                        .iter()
                        .filter(|angle| angle >= &&lo && angle < &&hi)
                        .count();
                    if type_i == type_j {
                        n /= 2;
                    }
                    ((lo, hi), n)
                })
                .collect::<Vec<_>>()
        })
        .fold(
            {
                bins.iter()
                    .map(|&(lo, hi)| ((lo, hi), 0))
                    .collect::<Vec<_>>()
            },
            |mut a, b| {
                iter::zip(&mut a, b).for_each(|(a, b)| a.1 += b.1);
                a
            },
        );
    normalize(adf, dump.atoms_count)
}

#[allow(clippy::cast_precision_loss)]
fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let dump_path = cli.dump_file;
    let timesteps = cli
        .timestep
        .map_or_else(Vec::new, |timestep| vec![timestep]);
    let dump = DumpFile::read(dump_path.as_path(), &timesteps)?;
    let snapshot = dump.get_snapshots()[0];
    let adf = get_adf(
        cli.type_i,
        cli.type_j,
        cli.type_k,
        cli.cutoff_i,
        cli.cutoff_j,
        cli.n_bins,
        snapshot,
    );
    let mut n_tot = 0.0;
    let table = adf
        .into_iter()
        .map(|(t, g, n)| {
            n_tot += n as f64;
            [t, g, n_tot]
                .into_iter()
                .map(|x| format!("{x:10.4}"))
                .join("\t")
        })
        .join("\n");
    println!("# theta g\n{table}");
    Ok(())
}

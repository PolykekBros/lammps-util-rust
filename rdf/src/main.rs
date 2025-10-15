use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use lammps_util_rust::{DumpFile, DumpSnapshot, SymBox, XYZ};
use rayon::prelude::*;
use std::{array, f64, iter, path::PathBuf};

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
    dump: &DumpSnapshot,
) -> Vec<(f64, f64)> {
    let vol = dump.sym_box.volume();
    let num = dump.atoms_count as f64;
    let rho = num / vol;
    rdf.into_iter()
        .map(|((lo, hi), n)| {
            let vshell = 4.0 / 3.0 * f64::consts::PI * (hi.powi(3) - lo.powi(3));
            let nnorm = rho * vshell;
            let n = n as f64 / (nnorm * num);
            ((hi + lo) / 2.0, n)
        })
        .collect()
}

fn get_supercell_coords(coords: &[XYZ], sym_box: &SymBox, cutoff: f64) -> Vec<XYZ> {
    let mut supercell_coords = coords.to_vec();
    let hi = [sym_box.xhi, sym_box.yhi, sym_box.zhi];
    let lo = [sym_box.xlo, sym_box.ylo, sym_box.zlo];
    let shift: [f64; 3] = array::from_fn(|i| hi[i] - lo[i]);
    (-1..=1)
        .flat_map(|px| (-1..=1).map(move |py| (px, py)))
        .flat_map(|(px, py)| (-1..=1).map(move |pz| [px, py, pz]))
        .filter(|periods| periods.iter().any(|&period| period != 0))
        .for_each(|periods| {
            let shift: [f64; 3] = array::from_fn(|i| shift[i] * periods[i] as f64);
            supercell_coords.extend(
                coords
                    .iter()
                    .filter(|atom| {
                        let flags: [bool; 3] = array::from_fn(|i| match periods[i] {
                            1 => atom[i] < lo[i] + cutoff,
                            -1 => atom[i] > hi[i] - cutoff,
                            _ => true,
                        });
                        flags.into_iter().all(|flag| flag)
                    })
                    .map(|atom| XYZ::from(array::from_fn(|i| atom[i] + shift[i]), 0)),
            );
        });
    supercell_coords
}

fn get_rdf(cutoff: f64, n: usize, dump: &DumpSnapshot) -> Vec<(f64, f64)> {
    let bins = get_bins(cutoff, n);
    let coords = dump.get_coordinates();
    let supercell_coords = get_supercell_coords(&coords, &dump.sym_box, cutoff);
    let kdtree = kd_tree::KdTree::build_by_ordered_float(supercell_coords);
    let rdf = coords
        .par_iter()
        .map(|atom| {
            let d_sq = kdtree
                .within_radius(atom, cutoff)
                .into_iter()
                .map(|neigh| atom.distance_squared(neigh))
                .collect::<Vec<_>>();
            bins.iter()
                .map(|&(lo, hi)| {
                    let lo_sq = lo.powi(2);
                    let hi_sq = hi.powi(2);
                    let n = d_sq
                        .iter()
                        .filter(|&&d_sq| d_sq >= lo_sq && d_sq < hi_sq)
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

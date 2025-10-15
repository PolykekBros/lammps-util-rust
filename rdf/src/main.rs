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

fn get_supercell_coords(coords: &mut Vec<XYZ>, sym_box: &SymBox, cutoff: f64) -> usize {
    let hi = [sym_box.xhi, sym_box.yhi, sym_box.zhi];
    let lo = [sym_box.xlo, sym_box.ylo, sym_box.zlo];
    let shift: [f64; 3] = array::from_fn(|i| hi[i] - lo[i]);
    let periods = (-1..=1)
        .flat_map(|px| (-1..=1).map(move |py| (px, py)))
        .flat_map(|(px, py)| (-1..=1).map(move |pz| [px, py, pz]))
        .filter(|periods| periods.iter().any(|&period| period != 0))
        .collect::<Vec<_>>();
    let shifts = periods
        .iter()
        .map(|periods| array::from_fn::<_, 3, _>(|i| shift[i] * periods[i] as f64))
        .collect::<Vec<_>>();
    let max_id = coords
        .iter()
        .map(|xyz| xyz.index())
        .max()
        .unwrap_or_default();
    let mut id = max_id;
    for atom_idx in 0..coords.len() {
        for (period, shift) in iter::zip(&periods, &shifts) {
            if (0..3).all(|i| match period[i] {
                1 => coords[atom_idx][i] < lo[i] + cutoff,
                -1 => coords[atom_idx][i] > hi[i] - cutoff,
                _ => true,
            }) {
                id += 1;
                coords.push(XYZ::from(
                    array::from_fn(|i| coords[atom_idx][i] + shift[i]),
                    id,
                ));
            }
        }
    }
    max_id
}

fn get_rdf(cutoff: f64, n: usize, dump: &DumpSnapshot) -> Vec<(f64, f64)> {
    let bins = get_bins(cutoff, n);
    let mut coords = dump.get_coordinates();
    let max_id = get_supercell_coords(&mut coords, &dump.sym_box, cutoff);
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords);
    let rdf = kdtree
        .items()
        .par_iter()
        .filter(|atom| atom.index() <= max_id)
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

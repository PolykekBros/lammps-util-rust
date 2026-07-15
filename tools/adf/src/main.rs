use anyhow::Result;
use clap::{Parser, Subcommand};
use itertools::Itertools;
use lammps_files::Dump;
use lammps_util::{xyz::xyz_vec_from_snapshot, XYZ};
use log::debug;
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

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Calculate ADF per slice along z
    Slice(SliceArgs),
}

#[derive(Parser)]
struct SliceArgs {
    /// Lowest z coordinate (default: minimum z)
    #[arg(long)]
    zlo: Option<f64>,

    /// Highest z coordinate (default: maximum z)
    #[arg(long)]
    zhi: Option<f64>,

    /// Slice width (default: 10)
    #[arg(long, default_value_t = 10.0)]
    d: f64,
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
    debug!("{int}");
    adf.into_iter()
        .map(|((lo, hi), val)| (lo.midpoint(hi), val as f64 / int, val / n))
        .collect()
}

#[allow(clippy::cast_possible_truncation)]
#[allow(clippy::cast_sign_loss)]
fn get_adf_for_slice(
    type_i: usize,
    type_j: usize,
    type_k: usize,
    cutoff_i: f64,
    cutoff_j: f64,
    n: usize,
    d_types: &[usize],
    tree: &kd_tree::KdTree3<XYZ>,
    z_lo: Option<f64>,
    z_hi: Option<f64>,
) -> Vec<(f64, f64, usize)> {
    let bins = get_bins(n);
    let n_slice_atoms: usize = tree
        .items()
        .iter()
        .filter(|atom| {
            d_types[atom.index] == type_k
                && z_lo.map_or(true, |lo| atom.coords[2] >= lo)
                && z_hi.map_or(true, |hi| atom.coords[2] < hi)
        })
        .count();
    if n_slice_atoms == 0 {
        return bins
            .iter()
            .map(|&(lo, hi)| (lo.midpoint(hi), 0.0, 0))
            .collect();
    }
    let adf = tree
        .items()
        .iter()
        .filter(|atom| {
            d_types[atom.index] == type_k
                && z_lo.map_or(true, |lo| atom.coords[2] >= lo)
                && z_hi.map_or(true, |hi| atom.coords[2] < hi)
        })
        .map(|atom_k| {
            let i_neigh = tree.within_radius(atom_k, cutoff_i);
            let j_neigh = {
                let mut j_neigh = tree.within_radius(atom_k, cutoff_j);
                j_neigh.retain(|atom_j| {
                    d_types[atom_j.index] == type_j && atom_j.index != atom_k.index
                });
                j_neigh
            };
            let j_neigh_ref = &j_neigh;
            let angles = i_neigh
                .iter()
                .filter(|atom_i| d_types[atom_i.index] == type_i && atom_i.index != atom_k.index)
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
    normalize(adf, n_slice_atoms)
}

fn format_adf_table(adf: Vec<(f64, f64, usize)>) -> String {
    let mut n_tot = 0.0;
    adf.into_iter()
        .map(|(t, g, n)| {
            n_tot += n as f64;
            [t, g, n_tot]
                .into_iter()
                .map(|x| format!("{x:10.4}"))
                .join("\t")
        })
        .join("\n")
}

fn get_kdtree<P>(
    dump_path: P,
    timesteps: &[u64],
    cutoff_i: f64,
    cutoff_j: f64,
) -> Result<(kd_tree::KdTree3<XYZ>, Vec<usize>)>
where
    P: AsRef<std::path::Path>,
{
    let dump = Dump::open(dump_path, &timesteps)?;
    let snapshot = &dump.get_snapshots()[0];
    let d_types = snapshot
        .get_property("type")
        .iter()
        .copied()
        .map(|t| t as usize)
        .collect::<Vec<_>>();
    let mut coords = xyz_vec_from_snapshot(snapshot);
    XYZ::get_supercell_coords(&mut coords, snapshot.get_symbox(), cutoff_i.max(cutoff_j));
    Ok((kd_tree::KdTree::build_by_ordered_float(coords), d_types))
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let dump_path = cli.dump_file;

    let mut header = format!(
        "# adf --n-bins {} --type-i {} --type-j {} --type-k {} --cutoff-i {} --cutoff-j {}",
        cli.n_bins, cli.type_i, cli.type_j, cli.type_k, cli.cutoff_i, cli.cutoff_j,
    );

    if let Some(timestep) = cli.timestep {
        header.push_str(&format!(" --timestep {}", timestep));
    }

    match &cli.command {
        None => {}
        Some(Commands::Slice(slice_args)) => {
            header.push_str(" slice");
            if let Some(zlo) = slice_args.zlo {
                header.push_str(&format!(" --zlo {}", zlo));
            }
            if let Some(zhi) = slice_args.zhi {
                header.push_str(&format!(" --zhi {}", zhi));
            }
            header.push_str(&format!(" --d {}", slice_args.d));
        }
    }

    header.push_str(&format!(" {}", dump_path.display()));
    println!("{header}");
    println!("# theta g");
    let timesteps = cli
        .timestep
        .map_or_else(Vec::new, |timestep| vec![timestep]);
    let (kdtree, d_type) = get_kdtree(dump_path, &timesteps, cli.cutoff_i, cli.cutoff_j)?;
    match cli.command {
        None => {
            let adf = get_adf_for_slice(
                cli.type_i,
                cli.type_j,
                cli.type_k,
                cli.cutoff_i,
                cli.cutoff_j,
                cli.n_bins,
                &d_type,
                &kdtree,
                None,
                None,
            );
            let table = format_adf_table(adf);
            println!("{table}");
        }
        Some(Commands::Slice(slice_args)) => {
            let zlo = slice_args.zlo.unwrap_or_else(|| {
                kdtree
                    .items()
                    .iter()
                    .map(|a| a.coords[2])
                    .fold(f64::INFINITY, f64::min)
            });
            let zhi = slice_args.zhi.unwrap_or_else(|| {
                kdtree
                    .items()
                    .iter()
                    .map(|a| a.coords[2])
                    .fold(f64::NEG_INFINITY, f64::max)
            });
            let d = slice_args.d;

            let mut slice_idx = 0;
            let mut current_zlo = zlo;
            while current_zlo < zhi {
                let slice_zhi = current_zlo + d;
                let adf = get_adf_for_slice(
                    cli.type_i,
                    cli.type_j,
                    cli.type_k,
                    cli.cutoff_i,
                    cli.cutoff_j,
                    cli.n_bins,
                    &d_type,
                    &kdtree,
                    Some(current_zlo),
                    Some(slice_zhi),
                );
                let table = format_adf_table(adf);
                println!(
                    "# slice: {} {} {}\n{table}",
                    slice_idx, current_zlo, slice_zhi
                );
                current_zlo += d;
                slice_idx += 1;
            }
        }
    }
    Ok(())
}

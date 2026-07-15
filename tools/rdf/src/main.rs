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
use clap::Subcommand;
use kd_tree::KdTree3;
use lammps_files::{snapshot::SymBox, Dump, Snapshot};
use lammps_util::{xyz::xyz_vec_from_snapshot, XYZ};
// use rayon::prelude::*;
use std::{
    fmt::{self, Write},
    iter,
    path::PathBuf,
};

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

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Calculate RDF per slice along z
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

#[derive(Clone, Copy)]
struct RdfQuery {
    type_i: Option<usize>,
    type_j: Option<usize>,
    z_lo: Option<f64>,
    z_hi: Option<f64>,
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
    symbox: &SymBox,
    ni: f64,
    nj: f64,
) {
    let vol = symbox.bbox.volume();
    let rho = nj / vol;
    iter::zip(g, bins).for_each(|(n, [lo, hi])| {
        let vshell = 4.0 / 3.0 * std::f64::consts::PI * (hi.powi(3) - lo.powi(3));
        let nnorm = rho * vshell;
        *n /= nnorm * ni;
    });
}

fn get_kdtree(snapshot: &Snapshot, cutoff: f64) -> (KdTree3<XYZ>, Vec<usize>, &SymBox) {
    let d_types = snapshot
        .get_property("type")
        .iter()
        .copied()
        .map(|t| t as usize)
        .collect();
    let mut coords = xyz_vec_from_snapshot(snapshot);
    let symbox = snapshot.get_symbox();
    XYZ::get_supercell_coords(&mut coords, symbox, cutoff);
    (
        kd_tree::KdTree::build_by_ordered_float(coords),
        d_types,
        symbox,
    )
}

fn get_rdf_new(
    kdtree: &KdTree3<XYZ>,
    d_types: &[usize],
    n: usize,
    cutoff: f64,
    type_pairs: &[(usize, usize)],
    query: RdfQuery,
    symbox: &SymBox,
) -> Vec<f64> {
    let bins_vec = get_bins(cutoff, n);
    assert_eq!(bins_vec.len(), n * 2);
    let (bins, rem) = bins_vec.as_chunks::<2>();
    assert_eq!(rem.len(), 0);
    assert_eq!(bins.len(), n);
    if type_pairs.is_empty() {
        let cols = 1 + 1;
        let mut table = make_table(bins, cols);
        let g = calculate_rdf(
            kdtree,
            bins,
            cutoff,
            d_types,
            RdfQuery { type_i: None, type_j: None, ..query },
            symbox,
        );
        assert_eq!(g.len(), bins.len());
        for (idx, g) in g.into_iter().enumerate() {
            table[idx * 2 + 1] = g;
        }
        table
    } else {
        let cols = type_pairs.len() + 1;
        let mut table = make_table(bins, cols);
        for (col, (ti, tj)) in type_pairs.iter().enumerate() {
            let g = calculate_rdf(
                kdtree,
                bins,
                cutoff,
                d_types,
                RdfQuery { type_i: Some(*ti), type_j: Some(*tj), ..query },
                symbox,
            );
            assert_eq!(g.len(), bins.len());
            for (row, g) in g.into_iter().enumerate() {
                table[row * cols + col + 1] = g;
            }
        }
        table
    }
}

fn calculate_rdf(
    kdtree: &KdTree3<XYZ>,
    bins: &[[f64; 2]],
    cutoff: f64,
    d_types: &[usize],
    query: RdfQuery,
    symbox: &SymBox,
) -> Vec<f64> {
    let type_i = query.type_i;
    let type_j = query.type_j;
    let z_lo = query.z_lo;
    let z_hi = query.z_hi;

    let mut g = kdtree
        .items()
        .iter()
        .filter(|atom| !atom.is_ghost)
        .filter(|atom| type_i.is_none_or(|ti| d_types[atom.index] == ti))
        .filter(|atom| {
            z_lo.is_none_or(|lo| atom.coords[2] >= lo) && z_hi.is_none_or(|hi| atom.coords[2] < hi)
        })
        .map(|atom| calculate_rdf_hist(kdtree, bins, cutoff, atom, d_types, type_j))
        .fold(vec![0.0; bins.len()], |mut g, part_g| {
            iter::zip(&mut g, part_g).for_each(|(g, p_g)| *g += p_g);
            g
        });
    let ni = kdtree
        .items()
        .iter()
        .filter(|atom| !atom.is_ghost)
        .filter(|atom| type_i.is_none_or(|ti| d_types[atom.index] == ti))
        .filter(|atom| {
            z_lo.is_none_or(|lo| atom.coords[2] >= lo) && z_hi.is_none_or(|hi| atom.coords[2] < hi)
        })
        .count() as f64;
    let nj = kdtree
        .items()
        .iter()
        .filter(|atom| !atom.is_ghost)
        .filter(|atom| type_j.is_none_or(|tj| d_types[atom.index] == tj))
        .count() as f64;
    normalize_g(&mut g, bins, symbox, ni, nj);
    g
}

fn calculate_rdf_hist(
    kdtree: &KdTree3<XYZ>,
    bins: &[[f64; 2]],
    cutoff: f64,
    atom: &XYZ,
    d_types: &[usize],
    type_j: Option<usize>,
) -> Vec<f64> {
    let d_sq = kdtree
        .within_radius(atom, cutoff)
        .iter()
        .filter(|neigh| type_j.is_none_or(|tj| d_types[neigh.index] == tj))
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
    let cli = Cli::parse();
    let timesteps = cli.timestep.map(|t| vec![t]).unwrap_or_default();
    let dump = Dump::open(&cli.dump_file, &timesteps).unwrap();
    let snapshot = &dump.get_snapshots()[0];

    // Reproducibility header
    let mut header = format!("# rdf --n-bins {} --cutoff {}", cli.n_bins, cli.cutoff,);
    if !cli.types.is_empty() {
        let types_str = cli
            .types
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",");
        write!(&mut header, " --types {types_str}").unwrap();
    }
    if let Some(timestep) = cli.timestep {
        write!(&mut header, " --timestep {timestep}").unwrap();
    }
    match &cli.command {
        None => {}
        Some(Commands::Slice(slice_args)) => {
            header.push_str(" slice");
            if let Some(zlo) = slice_args.zlo {
                write!(&mut header, " --zlo {zlo}").unwrap();
            }
            if let Some(zhi) = slice_args.zhi {
                write!(&mut header, " --zhi {zhi}").unwrap();
            }
            write!(&mut header, " --d {}", slice_args.d).unwrap();
        }
    }
    write!(&mut header, " {}", cli.dump_file.display()).unwrap();
    println!("{header}");

    let type_pairs = get_type_pairs(cli.types);
    let table_header = make_rdf_table_header(&type_pairs);

    let (kdtree, d_types, symbox) = get_kdtree(snapshot, cli.cutoff);
    match cli.command {
        None => {
            let table = get_rdf_new(
                &kdtree,
                &d_types,
                cli.n_bins,
                cli.cutoff,
                &type_pairs,
                RdfQuery { type_i: None, type_j: None, z_lo: None, z_hi: None },
                symbox,
            );
            print_table(
                &table,
                &table_header,
                table.len() / table_header.len(),
                table_header.len(),
            );
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
                let table = get_rdf_new(
                    &kdtree,
                    &d_types,
                    cli.n_bins,
                    cli.cutoff,
                    &type_pairs,
                    RdfQuery {
                        type_i: None,
                        type_j: None,
                        z_lo: Some(current_zlo),
                        z_hi: Some(slice_zhi),
                    },
                    symbox,
                );
                println!("# slice: {slice_idx} {current_zlo} {slice_zhi}");
                print_table(
                    &table,
                    &table_header,
                    table.len() / table_header.len(),
                    table_header.len(),
                );
                current_zlo += d;
                slice_idx += 1;
            }
        }
    }
}

use anyhow::Result;
use clap::Parser;
use itertools::Itertools;
use lammps_util_rust::{DumpFile, DumpSnapshot};
use std::{f32::consts::PI, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    dump_file: PathBuf,

    #[arg(short, long)]
    timestep: Option<u64>,

    #[arg(short, long)]
    cutoff: f64,

    #[arg(short, long)]
    binsize: f64,

    #[arg(short, long)]
    output_dir: PathBuf,
}

fn get_bins(cutoff: f64, binsize: f64) -> Vec<f64> {
    let bin_num = (cutoff / binsize).floor() as usize;
    (1..=bin_num)
        .map(|bin| (bin as f64 - 0.5) * binsize)
        .collect()
}

fn normalize(rdf: &mut [(f64, f64)], binsize: f64, dump: &DumpSnapshot) {
    let vol = dump.sym_box.volume();
    let num = dump.atoms_count as f64;
    println!("{}", num);
    let rho = num / vol;
    rdf.iter_mut()
        .for_each(|(r, n)| *n /= 4.0 * PI as f64 * (r.powi(2)) * binsize * rho);
}

fn get_rdf(cutoff: f64, binsize: f64, dump: &DumpSnapshot) -> Vec<(f64, f64)> {
    let bins = get_bins(cutoff, binsize);
    let mut rdf = bins
        .into_iter()
        .map(|bin| (bin, 0.0))
        .collect::<Vec<(f64, f64)>>();
    let coords = dump.get_coordinates();
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords.clone());
    for atom in coords {
        rdf[1..].iter_mut().for_each(|(r, n)| {
            let neigh = kdtree.within_radius(&atom, *r + 0.5 * binsize).len();
            let n_cur = neigh - kdtree.within_radius(&atom, *r - 0.5 * binsize).len();
            *n += n_cur as f64;
        });
    }
    normalize(&mut rdf, binsize, dump);
    rdf
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
    let rdf = get_rdf(cli.cutoff, cli.binsize, snapshot);
    let table = rdf
        .into_iter()
        .map(|(r, n)| {
            [r, n as f64]
                .into_iter()
                .map(|x| format!("{x:10.4}"))
                .join("\t")
        })
        .join("\n");
    println!("# radius number\n{table}");
    Ok(())
}

#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]
use anyhow::Result;
use clap::Parser;
use lammps_util_rust::{DumpFile, DumpSnapshot, RunDir, get_clusters, process_results_dir};
use std::{collections::HashMap, iter, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    results_dir: PathBuf,

    /// Atom types "<type 1>,<type 2>,...,<type N>", ex. "Si,C,O"
    #[arg(short, long)]
    particles: String,
}

struct Atom {
    _id: usize,
    atype: usize,
    _coords: [f64; 3],
    velocity: [f64; 3],
    mass: f64,
}

impl Atom {
    const fn new(id: usize, atype: usize, coords: [f64; 3], velocity: [f64; 3], mass: f64) -> Self {
        Self {
            _id: id,
            atype,
            _coords: coords,
            velocity,
            mass,
        }
    }
}

#[derive(Default)]
struct Cluster {
    mass: f64,
    counts: Vec<usize>,
    momentum: [f64; 3],
    ek: f64,
    angle: f64,
}

impl Cluster {
    fn new(atoms: &[Atom], types_map: &HashMap<usize, String>) -> Self {
        let mut counts = vec![0; types_map.len()];
        let mut mass = 0.0;
        let mut momentum = [0.0; 3];
        for atom in atoms {
            mass += atom.mass;
            iter::zip(&mut momentum, atom.velocity).for_each(|(m, v)| *m += v * atom.mass);
            counts[atom.atype - 1] += 1;
        }
        let ek =
            2.0 * 5.1875 * 1e-5 * momentum.iter().map(|m| m.powi(2)).sum::<f64>() / (2.0 * mass);
        let angle = (momentum[2] / momentum[0].hypot(momentum[1])).atan();
        Self {
            mass,
            counts,
            momentum,
            ek,
            angle,
        }
    }
}

fn parse_types(s: &str) -> (HashMap<usize, String>, Vec<String>) {
    let type_names = s
        .split(',')
        .map(|s| s.trim().to_string())
        .collect::<Vec<_>>();
    let types_map = type_names.iter().map(String::from).enumerate().collect();
    (types_map, type_names)
}

fn analyze_clusters(dump: &DumpSnapshot, types_map: &HashMap<usize, String>) -> Vec<Cluster> {
    let id = dump.get_property("id");
    let atype = dump.get_property("type");
    let x = dump.get_property("x");
    let y = dump.get_property("y");
    let z = dump.get_property("z");
    let vx = dump.get_property("vx");
    let vy = dump.get_property("vy");
    let vz = dump.get_property("vz");
    let mass = dump.get_property("mass");
    get_clusters(dump)
        .into_values()
        .map(|atoms_idx| {
            let atoms = atoms_idx
                .into_iter()
                .map(|atom_idx| {
                    Atom::new(
                        id[atom_idx] as usize,
                        atype[atom_idx] as usize,
                        [x[atom_idx], y[atom_idx], z[atom_idx]],
                        [vx[atom_idx], vy[atom_idx], vz[atom_idx]],
                        mass[atom_idx],
                    )
                })
                .collect::<Vec<_>>();
            Cluster::new(&atoms, types_map)
        })
        .collect()
}

fn do_single_dir(dir: &RunDir, types_map: &HashMap<usize, String>) -> Result<Vec<Cluster>> {
    let dump = DumpFile::read(&dir.path.join("dump.sputter"), &[])?;
    let clusters = analyze_clusters(dump.get_snapshots()[0], types_map);
    Ok(clusters)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let (types_map, type_names) = parse_types(&cli.particles);
    let clusters = process_results_dir(&cli.results_dir, |run_dir| {
        let result = do_single_dir(run_dir, &types_map);
        eprintln!(
            "{} {:?}",
            run_dir.path.to_string_lossy(),
            result.as_ref().map(Vec::len)
        );
        result
    })?;
    println!("# № {} ∑", type_names.join(" "));
    for (run_dir, clusters) in clusters {
        let counts = clusters.into_iter().map(|cluster| cluster.counts).fold(
            vec![0; types_map.len()],
            |mut agg, counts| {
                iter::zip(&mut agg, counts).for_each(|(a, b)| *a += b);
                agg
            },
        );
        let sum = counts.iter().sum::<usize>();
        print!("{}\t", run_dir.num);
        for count in counts {
            print!("{count}\t");
        }
        println!("{sum}");
    }
    Ok(())
}

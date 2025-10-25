use anyhow::Result;
use clap::Parser;
use lammps_util_rust::{DumpFile, DumpSnapshot, RunDir, process_results_dir};
use std::{collections::HashMap, iter, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    results_dir: PathBuf,

    /// Atom types "<type 1>,<type 2>,...,<type N>", ex. "Si,C,O"
    #[arg(short, long)]
    particles: String,

    /// Number of threads to run in parallel
    #[arg(short, long, default_value_t = 2)]
    threads: usize,
}

struct Atom {
    _id: usize,
    atype: usize,
    _coords: [f64; 3],
    velocity: [f64; 3],
    mass: f64,
}

impl Atom {
    fn new(id: usize, atype: usize, coords: [f64; 3], velocity: [f64; 3], mass: f64) -> Self {
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
    counts: HashMap<usize, usize>,
    momentum: [f64; 3],
    ek: f64,
    angle: f64,
}

impl Cluster {
    fn new(atoms: &[Atom], types_map: &HashMap<usize, String>) -> Self {
        let mut cluster = Cluster {
            counts: types_map.keys().map(|i| (*i, 0)).collect(),
            ..Default::default()
        };
        for atom in atoms {
            cluster.mass += atom.mass;
            iter::zip(&mut cluster.momentum, atom.velocity).for_each(|(m, v)| *m += v * atom.mass);
            cluster.counts.entry(atom.atype).and_modify(|n| *n += 1);
        }
        cluster.ek = 2.0 * 5.1875 * 1e-5 * cluster.momentum.iter().map(|m| m.powi(2)).sum::<f64>()
            / (2.0 * cluster.mass);
        cluster.angle = (cluster.momentum[2]
            / (cluster.momentum[0].powi(2) + cluster.momentum[1].powi(2)).sqrt())
        .atan();
        cluster
    }
}

fn parse_types(s: &str) -> (HashMap<usize, String>, Vec<String>) {
    let type_names = s
        .split(":")
        .map(|s| s.trim().to_string())
        .collect::<Vec<_>>();
    let types_map = type_names.iter().map(String::from).enumerate().collect();
    (types_map, type_names)
}

fn get_clusters(dump: &DumpSnapshot, types_map: &HashMap<usize, String>) -> Vec<Cluster> {
    let id = dump.get_property("id");
    let atype = dump.get_property("type");
    let x = dump.get_property("x");
    let y = dump.get_property("y");
    let z = dump.get_property("z");
    let vx = dump.get_property("vx");
    let vy = dump.get_property("vy");
    let vz = dump.get_property("vz");
    let mass = dump.get_property("mass");
    let cluster = dump.get_property("cluster");
    let mut clusters = HashMap::new();
    for i in 0..dump.atoms_count {
        clusters
            .entry(cluster[i] as usize)
            .and_modify(|v: &mut Vec<_>| {
                v.push(Atom::new(
                    id[i] as usize,
                    atype[i] as usize,
                    [x[i], y[i], z[i]],
                    [vx[i], vy[i], vz[i]],
                    mass[i],
                ));
            })
            .or_insert(Vec::new());
    }
    clusters
        .values()
        .map(|atoms| Cluster::new(atoms, types_map))
        .collect()
}

fn do_single_dir(dir: &RunDir, types_map: &HashMap<usize, String>) -> Result<Vec<Cluster>> {
    let dump = DumpFile::read(&dir.path.join("dump.sputter"), &[])?;
    let clusters = get_clusters(dump.get_snapshots()[0], types_map);
    Ok(clusters)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let (types_map, type_names) = parse_types(&cli.particles);
    let clusters = process_results_dir(&cli.results_dir, cli.threads, |run_dir| {
        do_single_dir(run_dir, &types_map)
    })?;
    println!("# № {} ∑", type_names.join(" "));
    for (run_dir, clusters) in clusters.into_iter() {
        for cluster in clusters.into_iter() {
            let sum = cluster.counts.values().copied().sum::<usize>();
            let mut counts = cluster.counts.iter().collect::<Vec<_>>();
            counts.sort_by(|a, b| a.0.cmp(b.0));
            let counts_s = counts
                .into_iter()
                .map(|(_, c)| c.to_string())
                .collect::<Vec<_>>()
                .join("\n");
            println!("{} {} {sum}", run_dir.num, counts_s);
        }
    }
    Ok(())
}

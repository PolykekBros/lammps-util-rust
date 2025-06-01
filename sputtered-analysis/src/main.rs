use anyhow::Result;
use clap::Parser;
use lammps_util_rust::{DumpFile, DumpSnapshot, RunDir, get_runs_dirs};
use rayon::{ThreadPoolBuilder, prelude::*};
use std::{
    collections::HashMap,
    path::{Path, PathBuf},
};

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
    fn new(atoms: &[Atom]) -> Self {
        let mut cluster = Cluster::default();
        for atom in atoms {
            cluster.mass += atom.mass;
            for (i, vel) in atom.velocity.iter().enumerate() {
                cluster.momentum[i] += vel * atom.mass;
            }
            cluster
                .counts
                .entry(atom.atype)
                .and_modify(|n| *n += 1)
                .or_default();
        }
        cluster.ek = 2.0
            * 5.1875
            * 1e-5
            * (cluster.momentum[0].powi(2)
                + cluster.momentum[1].powi(2)
                + cluster.momentum[2].powi(2))
            / (2.0 * cluster.mass);
        cluster.angle = (cluster.momentum[2]
            / (cluster.momentum[0].powi(2) + cluster.momentum[1].powi(2)).sqrt())
        .atan();
        cluster
    }
}

fn parse_types(s: &str) -> Vec<String> {
    s.split(":").map(|s| s.trim().to_string()).collect()
}

fn get_clusters(dump: &DumpSnapshot) -> Vec<Cluster> {
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
    clusters.values().map(|atoms| Cluster::new(atoms)).collect()
}

fn do_single_dir(dir: RunDir) -> Result<(usize, Vec<Cluster>)> {
    let dump = DumpFile::read(&dir.path.join("dump.sputter"), &[])?;
    let clusters = get_clusters(dump.get_snapshots()[0]);
    Ok((dir.num, clusters))
}

fn do_results_dir(results_dir: &Path, threads: usize) -> Result<Vec<(usize, Vec<Cluster>)>> {
    let tp = ThreadPoolBuilder::new().num_threads(threads).build()?;
    let run_dirs = get_runs_dirs(results_dir)?;
    let mut results = tp.install(|| {
        run_dirs
            .into_par_iter()
            .map(do_single_dir)
            .collect::<Result<Vec<_>>>()
    })?;
    results.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(results)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let type_names = parse_types(&cli.particles);
    let clusters = do_results_dir(&cli.results_dir, cli.threads)?;
    println!("# № {} ∑", type_names.join(" "));
    for (sim_num, clusters) in clusters.into_iter() {
        for cluster in clusters.into_iter() {
            let sum = cluster.counts.values().copied().sum::<usize>();
            let mut counts = cluster.counts.iter().collect::<Vec<_>>();
            counts.sort_by(|a, b| a.0.cmp(b.0));
            let counts_s = counts
                .into_iter()
                .map(|(_, c)| c.to_string())
                .collect::<Vec<_>>()
                .join("\n");
            println!("{sim_num} {} {sum}", counts_s);
        }
    }
    Ok(())
}

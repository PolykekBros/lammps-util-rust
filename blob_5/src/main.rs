use anyhow::{Result, bail};
use clap::Parser;
use lammps_util_rust::{DumpFile, RunDir, get_avg_with_std, get_runs_dirs};
use rayon::{ThreadPoolBuilder, prelude::*};
use std::path::{Path, PathBuf};

const N: usize = 5;

#[derive(Clone)]
struct Particle {
    pos: [f64; 3],
    _vel: [f64; 3],
    ek: f64,
    id: usize,
}

impl Particle {
    fn new(pos: [f64; 3], vel: [f64; 3], ek: f64, id: usize) -> Self {
        Self {
            pos,
            _vel: vel,
            ek,
            id,
        }
    }
}

struct Timestep {
    particles: Vec<Particle>,
    _step: usize,
}

impl Timestep {
    fn new(particles: Vec<Particle>, step: usize) -> Self {
        Self {
            particles,
            _step: step,
        }
    }

    fn get_ids_sorted_by_z(&self) -> Vec<usize> {
        let mut particles = self.particles.clone();
        particles.sort_by(|a, b| a.pos[2].total_cmp(&b.pos[2]));
        particles.into_iter().map(|p| p.id).collect()
    }

    fn get_energies_by_ids(&self, ids: &[usize]) -> Vec<f64> {
        self.particles
            .iter()
            .filter(|p| ids.contains(&p.id))
            .map(|p| p.ek)
            .collect()
    }

    fn get_bottom_ids(&self) -> [usize; N] {
        let ids_sorted = self.get_ids_sorted_by_z();
        let mut ids = [0; N];
        ids.clone_from_slice(&ids_sorted[..N]);
        ids
    }

    fn get_top_ids(&self) -> [usize; N] {
        let ids_sorted = self.get_ids_sorted_by_z();
        let mut ids = [0; N];
        ids.clone_from_slice(&ids_sorted[ids_sorted.len() - N..]);
        ids
    }
}

struct Run {
    timesteps: Vec<Timestep>,
    _num: usize,
}

impl Run {
    fn new(timesteps: Vec<Timestep>, num: usize) -> Self {
        Self {
            timesteps,
            _num: num,
        }
    }

    fn get_avg_energies_by_ids(&self, ids: &[usize]) -> Vec<f64> {
        self.timesteps
            .iter()
            .map(|r| r.get_energies_by_ids(ids).iter().sum::<f64>() / N as f64)
            .collect()
    }

    fn get_avg_bottom_energies(&self) -> Vec<f64> {
        self.timesteps
            .first()
            .map(|timestep| {
                let ids = timestep.get_bottom_ids();
                self.get_avg_energies_by_ids(&ids)
            })
            .unwrap_or_default()
    }

    fn get_avg_top_energies(&self) -> Vec<f64> {
        self.timesteps
            .first()
            .map(|timestep| {
                let ids = timestep.get_top_ids();
                self.get_avg_energies_by_ids(&ids)
            })
            .unwrap_or_default()
    }
}

fn process_run_dir(run_dir: RunDir) -> Result<Run> {
    let dump = DumpFile::read(&run_dir.path.join("dump.during"), &[])?;
    let timesteps = dump
        .get_snapshots()
        .iter()
        .map(|s| {
            let x = s.get_property("x");
            let y = s.get_property("y");
            let z = s.get_property("z");
            let vx = s.get_property("vx");
            let vy = s.get_property("vy");
            let vz = s.get_property("vz");
            let ek = s.get_property("c_atom_ke");
            let id = s.get_property("id");
            let particles = (0..s.atoms_count)
                .map(|i| {
                    let pos = [x[i], y[i], z[i]];
                    let vel = [vx[i], vy[i], vz[i]];
                    Particle::new(pos, vel, ek[i], id[i] as usize)
                })
                .collect();
            Timestep::new(particles, s.step as usize)
        })
        .collect();
    Ok(Run::new(timesteps, run_dir.num))
}

fn data_transpose(data: &[Vec<f64>]) -> Result<Vec<Vec<f64>>> {
    let m = data.first().map(|inner| inner.len()).unwrap_or_default();
    for inner in data.iter() {
        let inner_m = inner.len();
        if inner_m != m {
            bail!("Inner dimensions do not much: {m} != {inner_m}");
        }
    }
    Ok((0..m)
        .map(|i| data.iter().map(|inner| inner[i]).collect())
        .collect())
}

fn parse_data(data: &[Vec<f64>]) -> Result<Vec<(f64, f64)>> {
    Ok(data_transpose(data)?
        .into_iter()
        .filter_map(|inner| get_avg_with_std(&inner))
        .collect())
}

fn get_bottom_data(runs: &[Run]) -> Result<Vec<(f64, f64)>> {
    let data = runs
        .iter()
        .map(|r| r.get_avg_bottom_energies())
        .collect::<Vec<_>>();
    parse_data(&data)
}

fn get_top_data(runs: &[Run]) -> Result<Vec<(f64, f64)>> {
    let data = runs
        .iter()
        .map(|r| r.get_avg_top_energies())
        .collect::<Vec<_>>();
    parse_data(&data)
}

fn get_data(results_dir: &Path, threads: usize) -> Result<Vec<Run>> {
    let tp = ThreadPoolBuilder::new().num_threads(threads).build()?;
    let dirs = get_runs_dirs(results_dir)?;
    tp.install(|| {
        dirs.into_par_iter()
            .map(process_run_dir)
            .collect::<Result<Vec<_>>>()
    })
}

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Results dir
    #[arg(short, long)]
    results_dir: PathBuf,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 2)]
    threads: usize,
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let data = get_data(&cli.results_dir, cli.threads)?;
    let data_bottom = get_bottom_data(&data)?;
    let data_top = get_top_data(&data)?;
    println!("# N ek_avg_bottom std ek_avg_top std");
    for i in 0..data_bottom.len() {
        let bottom = data_bottom[i];
        let top = data_top[i];
        println!(
            "{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            i as f32 * 10.0,
            bottom.0,
            bottom.1,
            top.0,
            top.1
        );
    }
    Ok(())
}

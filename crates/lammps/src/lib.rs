#![allow(clippy::cast_precision_loss)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::missing_errors_doc)]
use std::marker::PhantomData;
use std::path::{Path, PathBuf};
use anyhow::Result;

pub trait DirTask: Clone {
    type Output;
    fn run<P: AsRef<Path>>(&self, dir: P) -> Result<Self::Output>;
}

pub trait Task {
    type Output;
    fn run(&self) -> Result<Self::Output>;
}

pub struct MainWrapper<T> {
    cli: PhantomData<T>,
}

impl<T> Default for MainWrapper<T> {
    fn default() -> Self {
        Self { cli: PhantomData }
    }
}

impl<T: clap::Parser> MainWrapper<T> {
    pub fn run<F: FnMut(T) -> Box<dyn Task<Output = ()>>>(self, f: F) -> Result<()> {
        env_logger::init();
        let cli = T::parse();
        let mut f = f;
        f(cli).run()
    }
}

pub struct SingleTask<T> {
    pub run_dir: PathBuf,
    pub task: T,
}

impl<T> SingleTask<T> {
    pub fn new(run_dir: PathBuf, task: T) -> Self {
        Self { run_dir, task }
    }
}

impl<T> Task for SingleTask<T>
where
    T: DirTask,
{
    type Output = T::Output;
    fn run(&self) -> Result<Self::Output> {
        self.task.run(&self.run_dir)
    }
}
mod clusterizer;
mod math;
pub mod xyz;

use lammps_files::Snapshot;
use log::debug;
use rayon::prelude::*;
use std::{
    fs::read_dir,
    io,
};

pub use clusterizer::{clusterize_snapshot, get_clusters};
pub use geomutil_util;
pub use math::{range, IteratorAvg};
pub use xyz::XYZ;

pub use lammps_files::snapshot::copy_snapshot_with_indices;
pub use lammps_files::Snapshot as DumpSnapshot;

#[derive(Clone)]
pub struct DumpFile(pub lammps_files::Dump);

impl DumpFile {
    pub fn new(snapshots: Vec<lammps_files::Snapshot>) -> Self {
        Self(lammps_files::Dump::new(snapshots))
    }

    pub fn read<P: AsRef<Path>>(path: P, timesteps: &[u64]) -> Result<Self> {
        let dump = lammps_files::Dump::open(path, timesteps)?;
        Ok(Self(dump))
    }
}

impl std::ops::Deref for DumpFile {
    type Target = lammps_files::Dump;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for DumpFile {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

use crate::xyz::xyz_vec_from_snapshot;

pub struct RunDir {
    pub path: PathBuf,
    pub num: usize,
}

impl RunDir {
    const fn new(path: PathBuf, num: usize) -> Self {
        Self { path, num }
    }
}

pub fn get_runs_dirs(results_dir: &Path) -> io::Result<impl Iterator<Item = io::Result<RunDir>>> {
    Ok(read_dir(results_dir)?.filter_map(|e| {
        e.map(|e| {
            let p = e.path();
            let i = p
                .file_name()?
                .to_str()?
                .strip_prefix("run_")?
                .parse::<usize>()
                .ok()?;
            Some(RunDir::new(p, i))
        })
        .transpose()
    }))
}

pub fn process_results_dir<T, F>(dir: &Path, processor: F) -> Result<Vec<(RunDir, T)>>
where
    T: Send,
    F: Fn(&RunDir) -> Result<T> + Send + Sync,
{
    let run_dirs_iter = get_runs_dirs(dir)?.collect::<Vec<_>>();
    let mut results = run_dirs_iter
        .into_par_iter()
        .map(|run_dir| {
            let run_dir = run_dir?;
            processor(&run_dir).map(|r| (run_dir, r))
        })
        .collect::<Result<Vec<_>>>()?;
    results.sort_by_key(|a| a.0.num);
    Ok(results)
}

fn crater_candidates_snapshot(
    initial_snapshot: &Snapshot,
    final_snapshot: &Snapshot,
    candidate_cutoff: f64,
    cluster_cutoff: f64,
) -> Snapshot {
    let initial_coords = xyz_vec_from_snapshot(initial_snapshot);
    let final_coords = xyz_vec_from_snapshot(final_snapshot);
    let kdtree = kd_tree::KdTree::build_by_ordered_float(final_coords);
    let mut indices = Vec::new();
    for atom in initial_coords {
        if kdtree.within_radius(&atom, candidate_cutoff).is_empty() {
            indices.push(atom.index);
        }
    }
    let candidates_snapshot = copy_snapshot_with_indices(initial_snapshot, indices.into_iter());
    debug!(
        "crater candidates atom count: {}",
        candidates_snapshot.get_atoms_count()
    );
    clusterize_snapshot(&candidates_snapshot, cluster_cutoff)
}

#[must_use]
pub fn crater_snapshot(
    initial_snapshot: &Snapshot,
    final_snapshot: &Snapshot,
    candidate_cutoff: f64,
    cluster_cutoff: f64,
) -> Snapshot {
    let candidates_snapshot = crater_candidates_snapshot(
        initial_snapshot,
        final_snapshot,
        candidate_cutoff,
        cluster_cutoff,
    );
    let clusters = get_clusters(&candidates_snapshot);
    let atoms = clusters
        .into_values()
        .max_by(|a, b| a.len().cmp(&b.len()))
        .unwrap();
    copy_snapshot_with_indices(&candidates_snapshot, atoms.into_iter())
}

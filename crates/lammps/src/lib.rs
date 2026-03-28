#![allow(clippy::cast_precision_loss)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::missing_errors_doc)]
mod clusterizer;
mod math;
pub mod xyz;

use anyhow::Result;
use lammps_files::{snapshot::copy_snapshot_with_indices, Snapshot};
use log::debug;
use rayon::prelude::*;
use std::{
    fs::read_dir,
    io,
    path::{Path, PathBuf},
};

pub use clusterizer::{clusterize_snapshot, get_clusters};
pub use geomutil_util;
pub use math::{range, IteratorAvg};
pub use xyz::XYZ;

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

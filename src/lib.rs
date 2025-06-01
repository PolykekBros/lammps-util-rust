mod clusterizer;
mod dump_file;
mod dump_snapshot;
mod xyz;

use anyhow::Result;
use log::debug;
use regex::Regex;
use std::{
    fs::read_dir,
    path::{Path, PathBuf},
};

pub use clusterizer::{clusterize_snapshot, get_cluster_counts, get_max_cluster_id};
pub use dump_file::DumpFile;
pub use dump_snapshot::{
    copy_snapshot, copy_snapshot_with_indices, copy_snapshot_with_indices_with_keys,
    copy_snapshot_with_keys, DumpSnapshot, SymBox,
};
pub use xyz::{check_cutoff, XYZ};

pub struct RunDir {
    pub path: PathBuf,
    pub num: usize,
}

impl RunDir {
    fn new(path: PathBuf, num: usize) -> Self {
        Self { path, num }
    }
}

pub fn get_runs_dirs(results_dir: &Path) -> Result<Vec<RunDir>> {
    let re = Regex::new(r"^run_(\d+)$")?;
    let paths = read_dir(results_dir)?.collect::<Result<Vec<_>, _>>()?;
    let dirs = paths
        .into_iter()
        .map(|e| e.path())
        .filter_map(|p| {
            p.file_name()
                .and_then(|name| name.to_str())
                .and_then(|s| re.captures(s))
                .and_then(|caps| caps.get(1))
                .and_then(|m| m.as_str().parse::<usize>().ok())
                .map(|n| RunDir::new(p, n))
        })
        .collect();
    Ok(dirs)
}

pub fn range_f64(begin: f64, end: f64, count: usize) -> Vec<f64> {
    assert!(begin <= end);
    let step = (end - begin) / 1.max(count - 1) as f64;
    (0..count).map(|n| begin + n as f64 * step).collect()
}

fn crater_candidates_snapshot(
    initial_snapshot: &DumpSnapshot,
    final_snapshot: &DumpSnapshot,
    candidate_cutoff: f64,
    cluster_cutoff: f64,
) -> DumpSnapshot {
    let initial_coords = initial_snapshot.get_coordinates();
    let final_coords = final_snapshot.get_coordinates();
    let kdtree = kd_tree::KdTree::build_by_ordered_float(final_coords.clone());
    let mut indices = Vec::new();
    for atom in initial_coords {
        if kdtree.within_radius(&atom, candidate_cutoff).is_empty() {
            indices.push(atom.index());
        }
    }
    let candidates_snapshot = copy_snapshot_with_indices(initial_snapshot, indices.into_iter());
    debug!(
        "crater candidates atom count: {}",
        candidates_snapshot.atoms_count
    );
    clusterize_snapshot(&candidates_snapshot, cluster_cutoff)
}

pub fn crater_snapshot(
    initial_snapshot: &DumpSnapshot,
    final_snapshot: &DumpSnapshot,
    candidate_cutoff: f64,
    cluster_cutoff: f64,
) -> DumpSnapshot {
    let candidates_snapshot = &crater_candidates_snapshot(
        initial_snapshot,
        final_snapshot,
        candidate_cutoff,
        cluster_cutoff,
    );
    let max_cluster = get_max_cluster_id(candidates_snapshot);
    let cluster = candidates_snapshot.get_property("cluster");
    let indices = cluster
        .iter()
        .enumerate()
        .filter(|(_, &cluster)| cluster as usize == max_cluster)
        .map(|(i, _)| i);
    copy_snapshot_with_indices(candidates_snapshot, indices)
}

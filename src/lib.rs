mod clusterizer;
mod dump_file;
mod dump_snapshot;
mod math;
mod xyz;

use log::debug;
use std::{
    fs::read_dir,
    io,
    path::{Path, PathBuf},
};

pub use clusterizer::{clusterize_snapshot, get_cluster_counts, get_max_cluster_id};
pub use dump_file::DumpFile;
pub use dump_snapshot::{
    copy_snapshot, copy_snapshot_with_indices, copy_snapshot_with_indices_with_keys,
    copy_snapshot_with_keys, DumpSnapshot, SymBox,
};
pub use math::{get_avg, get_avg_with_std, get_std, range_f32, range_f64};
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

pub fn get_runs_dirs(results_dir: &Path) -> io::Result<impl Iterator<Item = io::Result<RunDir>>> {
    Ok(read_dir(results_dir)?.flat_map(|e| {
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

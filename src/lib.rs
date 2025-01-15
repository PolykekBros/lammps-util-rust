mod clusterizer;
mod dump_file;
mod dump_snapshot;
mod xyz;

use std::iter;

pub use clusterizer::{clusterize_snapshot, get_max_cluster};
pub use dump_file::DumpFile;
pub use dump_snapshot::{DumpSnapshot, SymBox};
pub use xyz::{check_cutoff, XYZ};

pub fn range_f64(begin: f64, end: f64, count: usize) -> Vec<f64> {
    assert!(begin <= end);
    let step = (end - begin) / 1.max(count - 1) as f64;
    (0..count).map(|n| begin + n as f64 * step).collect()
}

pub fn copy_snapshot(input_snapshot: &DumpSnapshot) -> DumpSnapshot {
    copy_snapshot_with_indices_with_keys(
        input_snapshot,
        iter::empty(),
        0..input_snapshot.atoms_count,
    )
}

pub fn copy_snapshot_with_indices(
    input_snapshot: &DumpSnapshot,
    indices: impl Iterator<Item = usize>,
) -> DumpSnapshot {
    copy_snapshot_with_indices_with_keys(input_snapshot, iter::empty(), indices)
}

pub fn copy_snapshot_with_keys<'a>(
    input_snapshot: &DumpSnapshot,
    additional_keys: impl Iterator<Item = &'a str>,
) -> DumpSnapshot {
    copy_snapshot_with_indices_with_keys(
        input_snapshot,
        additional_keys,
        0..input_snapshot.atoms_count,
    )
}

pub fn copy_snapshot_with_indices_with_keys<'a>(
    input_snapshot: &DumpSnapshot,
    additional_keys: impl Iterator<Item = &'a str>,
    indices: impl Iterator<Item = usize>,
) -> DumpSnapshot {
    let mut keys = input_snapshot.get_keys_map().clone();
    for key in additional_keys {
        keys.insert(key.to_string(), keys.len());
    }
    let indices = indices.collect::<Vec<_>>();
    let mut snapshot = DumpSnapshot::new(
        keys,
        input_snapshot.step,
        indices.len(),
        input_snapshot.sym_box.clone(),
    );
    for (new_i, i) in indices.into_iter().enumerate() {
        for (j, _) in input_snapshot.get_keys().iter().enumerate() {
            snapshot.set_atom_value(j, new_i, input_snapshot.get_atom_value(j, i));
        }
    }
    snapshot
}

fn crater_candidates_snapshot(
    initial_snapshot: &DumpSnapshot,
    final_snapshot: &DumpSnapshot,
    cutoff: f64,
) -> DumpSnapshot {
    let initial_coords = initial_snapshot.get_coordinates();
    let final_coords = final_snapshot.get_coordinates();
    let kdtree = kd_tree::KdTree::build_by_ordered_float(final_coords.clone());
    let mut indices = Vec::new();
    for atom in initial_coords {
        if kdtree.within_radius(&atom, cutoff).is_empty() {
            indices.push(atom.index());
        }
    }
    let candidates_snapshot = copy_snapshot_with_indices(initial_snapshot, indices.into_iter());
    clusterize_snapshot(&candidates_snapshot, 3.0)
}

pub fn crater_snapshot(
    initial_snapshot: &DumpSnapshot,
    final_snapshot: &DumpSnapshot,
    cutoff: f64,
) -> DumpSnapshot {
    let candidates_snapshot = &crater_candidates_snapshot(initial_snapshot, final_snapshot, cutoff);
    let max_cluster = get_max_cluster(candidates_snapshot);
    let cluster = candidates_snapshot.get_property("cluster");
    let indices = cluster
        .iter()
        .enumerate()
        .filter(|(_, &cluster)| cluster as usize == max_cluster)
        .map(|(i, _)| i);
    copy_snapshot_with_indices(candidates_snapshot, indices)
}

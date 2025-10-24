use log::debug;

use crate::{copy_snapshot_with_keys, DumpSnapshot, XYZ};
use std::{
    collections::{HashMap, HashSet},
    iter,
};

#[must_use]
pub fn clusterize_snapshot(snapshot: &DumpSnapshot, cutoff: f32) -> DumpSnapshot {
    assert!(cutoff >= 0.0);
    let mut snapshot = copy_snapshot_with_keys(snapshot, iter::once("cluster"));
    let coords = snapshot.get_coordinates();
    let clusters = clusterize_coords(&coords, cutoff);
    let cluster_j = snapshot.get_property_index("cluster");
    let id_j = snapshot.get_property_index("id");
    for atom_i in 0..snapshot.atoms_count {
        let (&cluster_i, _) = clusters
            .iter()
            .find(|(_, cluster)| cluster.contains(&atom_i))
            .unwrap();
        let cluster_id = snapshot.get_atom_value(id_j, cluster_i);
        snapshot.set_atom_value(cluster_j, atom_i, cluster_id);
    }
    snapshot
}

fn clusterize_coords(coords: &[XYZ], cutoff: f32) -> HashMap<usize, HashSet<usize>> {
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords.to_vec());
    let mut visited = vec![false; coords.len()];
    let mut map = HashMap::new();
    for atom in coords {
        if visited[atom.index] {
            continue;
        }
        visited[atom.index] = true;
        let cluster = map.entry(atom.index).or_insert_with(HashSet::new);
        let mut stack = vec![atom];
        while let Some(atom) = stack.pop() {
            cluster.insert(atom.index);
            for neigh in kdtree.within_radius(atom, cutoff) {
                if !visited[neigh.index] {
                    visited[neigh.index] = true;
                    stack.push(neigh);
                }
            }
        }
    }
    debug!("clusters {map:?}");
    map
}

#[must_use]
pub fn get_cluster_counts(snapshot: &DumpSnapshot) -> HashMap<usize, usize> {
    let clusters = snapshot.get_property("cluster");
    let mut cluster_cnt = HashMap::new();
    for cluster in clusters {
        #[allow(clippy::cast_possible_truncation)]
        #[allow(clippy::cast_sign_loss)]
        let cnt = cluster_cnt.entry(*cluster as usize).or_insert(0);
        *cnt += 1;
    }
    debug!("clusters count {cluster_cnt:?}");
    cluster_cnt
}

#[must_use]
pub fn get_max_cluster_id(snapshot: &DumpSnapshot) -> usize {
    let cluster_cnt = get_cluster_counts(snapshot);
    let (max_cluster, _) = cluster_cnt
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1))
        .expect("Cluster snapshot is empty");
    max_cluster
}

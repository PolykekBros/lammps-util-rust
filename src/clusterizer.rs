use log::debug;

use crate::{copy_snapshot_with_keys, DumpSnapshot, SymBox, XYZ};
use std::{
    collections::{HashMap, HashSet},
    iter,
};

#[must_use]
pub fn clusterize_snapshot(snapshot: &DumpSnapshot, cutoff: f64) -> DumpSnapshot {
    let mut snapshot = copy_snapshot_with_keys(snapshot, iter::once("cluster"));
    let coords = snapshot.get_coordinates();
    let clusters = clusterize_coords(coords, &snapshot.sym_box, cutoff);
    let cluster_property_idx = snapshot.get_property_index("cluster");
    let id_property_idx = snapshot.get_property_index("id");
    for (atom_idx, cluster_idx) in clusters.into_iter().flat_map(|(cluster_idx, atoms_idx)| {
        atoms_idx
            .into_iter()
            .map(move |atom_idx| (atom_idx, cluster_idx))
    }) {
        let cluster_id = snapshot.get_atom_value(id_property_idx, cluster_idx);
        snapshot.set_atom_value(cluster_property_idx, atom_idx, cluster_id);
    }
    snapshot
}

fn clusterize_coords(
    mut coords: Vec<XYZ>,
    sym_box: &SymBox,
    cutoff: f64,
) -> HashMap<usize, HashSet<usize>> {
    let mut visited = vec![false; coords.len()];
    XYZ::get_supercell_coords(&mut coords, sym_box, cutoff);
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords);
    let mut map = HashMap::new();
    for atom_idx in kdtree
        .items()
        .iter()
        .filter(|atom| !atom.is_ghost)
        .map(|atom| atom.index)
    {
        if visited[atom_idx] {
            continue;
        }
        visited[atom_idx] = true;
        let cluster = map.entry(atom_idx).or_insert_with(HashSet::new);
        let mut stack = vec![atom_idx];
        while let Some(atom_idx) = stack.pop() {
            cluster.insert(atom_idx);
            for neigh in kdtree.within_radius(&kdtree.items()[atom_idx], cutoff) {
                if !visited[neigh.index] {
                    visited[neigh.index] = true;
                    stack.push(neigh.index);
                }
            }
        }
    }
    debug!("clusters {map:?}");
    map
}

#[allow(clippy::cast_possible_truncation)]
#[allow(clippy::cast_sign_loss)]
#[must_use]
pub fn get_clusters(snapshot: &DumpSnapshot) -> HashMap<usize, HashSet<usize>> {
    let clusters = snapshot.get_property("cluster");
    let ids = snapshot.get_property("id");
    let cluster_idx_map = iter::zip(ids, clusters)
        .enumerate()
        .filter(|(_, (id, cluster))| id.eq(cluster))
        .map(|(atom_idx, (id, _))| (*id as usize, atom_idx))
        .collect::<HashMap<_, _>>();
    let mut clusters_map: HashMap<usize, HashSet<usize>> = HashMap::new();
    clusters
        .iter()
        .copied()
        .enumerate()
        .for_each(|(atom_idx, cluster)| {
            let cluster_idx = cluster_idx_map[&(cluster as usize)];
            let atoms = clusters_map.entry(cluster_idx).or_default();
            atoms.insert(atom_idx);
        });
    clusters_map
}

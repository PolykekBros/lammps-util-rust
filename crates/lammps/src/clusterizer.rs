use lammps_files::{
    snapshot::{copy_snapshot_with_keys, SymBox},
    Snapshot,
};
use log::debug;

use crate::{xyz::xyz_vec_from_snapshot, XYZ};
use std::{
    collections::{HashMap, HashSet},
    iter,
};

#[must_use]
pub fn clusterize_snapshot(snapshot: &Snapshot, cutoff: f64) -> Snapshot {
    let mut snapshot = copy_snapshot_with_keys(snapshot, iter::once("cluster"));
    let coords = xyz_vec_from_snapshot(&snapshot);
    let clusters = clusterize_coords(coords, snapshot.get_symbox(), cutoff);
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
    let mut indices = vec![0; coords.len()];
    XYZ::get_supercell_coords(&mut coords, sym_box, cutoff);
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords);
    kdtree
        .items()
        .iter()
        .enumerate()
        .filter(|(_, atom)| !atom.is_ghost)
        .for_each(|(i, atom)| {
            indices[atom.index] = i;
        });
    let mut map = HashMap::new();
    for (atom_idx, atom) in indices.iter().copied().map(|i| (i, kdtree.items()[i])) {
        if visited[atom.index] {
            continue;
        }
        visited[atom.index] = true;
        let cluster = map.entry(atom.index).or_insert_with(HashSet::new);
        let mut stack = vec![atom_idx];
        while let Some(atom_idx) = stack.pop() {
            let atom = kdtree.items()[atom_idx];
            cluster.insert(atom.index);
            for neigh in kdtree.within_radius(&atom, cutoff) {
                if !visited[neigh.index] {
                    visited[neigh.index] = true;
                    stack.push(indices[neigh.index]);
                }
            }
        }
    }
    debug!(
        "clusters {:?}",
        map.iter()
            .map(|(cluster_id, clusters)| (*cluster_id, clusters.len()))
            .collect::<HashMap<_, _>>()
    );
    map
}

#[allow(clippy::cast_possible_truncation)]
#[allow(clippy::cast_sign_loss)]
#[must_use]
pub fn get_clusters(snapshot: &Snapshot) -> HashMap<usize, HashSet<usize>> {
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

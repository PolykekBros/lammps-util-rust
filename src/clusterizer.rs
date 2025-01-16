use log::debug;

use crate::{copy_snapshot_with_keys, DumpSnapshot, XYZ};
use std::collections::{HashMap, HashSet};

pub fn clusterize_snapshot(snapshot: &DumpSnapshot, cutoff: f64) -> DumpSnapshot {
    assert!(cutoff >= 0.0);
    let mut snapshot = copy_snapshot_with_keys(snapshot, ["cluster"].into_iter());
    let coords = snapshot.get_coordinates();
    let clusters = get_clusters(&coords, cutoff);
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

fn get_clusters(coords: &Vec<XYZ>, cutoff: f64) -> HashMap<usize, HashSet<usize>> {
    let kdtree = kd_tree::KdTree::build_by_ordered_float(coords.clone());
    let mut visited = vec![false; coords.len()];
    let mut map = HashMap::new();
    for atom in coords {
        if visited[atom.index()] {
            continue;
        }
        visited[atom.index()] = true;
        let cluster = map.entry(atom.index()).or_insert(HashSet::new());
        let mut stack = vec![atom];
        while let Some(atom) = stack.pop() {
            cluster.insert(atom.index());
            for neigh in kdtree.within_radius(atom, cutoff) {
                if !visited[neigh.index()] {
                    visited[neigh.index()] = true;
                    stack.push(neigh);
                }
            }
        }
    }
    debug!("clusters {map:?}");
    map
}

pub fn get_max_cluster(snapshot: &DumpSnapshot) -> usize {
    let cluster = snapshot.get_property("cluster");
    let mut cluster_cnt = HashMap::new();
    for cluster in cluster {
        let cnt = cluster_cnt.entry(*cluster as usize).or_insert(0);
        *cnt += 1;
    }
    debug!("clusters count {cluster_cnt:?}");
    let (max_cluster, _) = cluster_cnt
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1))
        .expect("Cluster snapshot is empty");
    max_cluster
}

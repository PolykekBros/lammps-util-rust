use super::{DumpFile, DumpSnapshot, XYZ};
use std::collections::{HashMap, HashSet};

pub struct Clusterizer {
    snapshot: DumpSnapshot,
    x_j: usize,
    y_j: usize,
    z_j: usize,
    cutoff: f64,
}

impl Clusterizer {
    pub fn new(snapshot_input: &DumpSnapshot, cutoff: f64) -> Self {
        assert!(cutoff >= 0.0);
        let mut keys = snapshot_input.get_keys_map().clone();
        let cluster_j = keys.len();
        keys.entry("cluster".to_string()).insert_entry(cluster_j);
        let x_j = keys["x"];
        let y_j = keys["y"];
        let z_j = keys["z"];
        let snapshot = (0..snapshot_input.atoms_count).fold(
            DumpSnapshot::new(
                keys,
                snapshot_input.step,
                snapshot_input.atoms_count,
                snapshot_input.sym_box.clone(),
            ),
            |mut snapshot, i| {
                snapshot_input
                    .get_keys()
                    .iter()
                    .enumerate()
                    .for_each(|(j, _)| {
                        snapshot.set_atom_value(j, i, snapshot_input.get_atom_value(j, i));
                    });
                snapshot
            },
        );
        Self {
            snapshot,
            x_j,
            y_j,
            z_j,
            cutoff,
        }
    }

    fn get_atom_xyz(&self, atom_i: usize) -> XYZ {
        XYZ::from([
            self.snapshot.get_atom_value(self.x_j, atom_i),
            self.snapshot.get_atom_value(self.y_j, atom_i),
            self.snapshot.get_atom_value(self.z_j, atom_i),
        ])
    }

    pub fn clusterize(mut self) -> DumpFile {
        let atoms: Vec<XYZ> = (0..self.snapshot.atoms_count)
            .map(|i| self.get_atom_xyz(i))
            .collect();
        let atoms_map =
            atoms
                .iter()
                .copied()
                .enumerate()
                .fold(HashMap::new(), |mut atoms, (i, xyz)| {
                    atoms.entry(xyz).insert_entry(i);
                    atoms
                });
        let kdtree = kd_tree::KdTree::build_by_ordered_float(atoms.clone());
        let mut clusters = HashMap::new();
        let mut visited = vec![false; kdtree.len()];
        atoms.iter().enumerate().for_each(|(i, atom)| {
            if !visited[i] {
                let mut current_cluster = HashSet::new();
                let mut stack = vec![(i, atom)];
                visited[i] = true;
                while let Some((current_i, atom)) = stack.pop() {
                    current_cluster.insert(current_i);
                    kdtree
                        .within_radius(atom, self.cutoff)
                        .iter()
                        .for_each(|xyz| {
                            let neigh_i = atoms_map[xyz];
                            println!("{i} {neigh_i} {xyz:?} {current_cluster:?}");
                            if !visited[neigh_i] {
                                visited[neigh_i] = true;
                                stack.push((neigh_i, xyz));
                            }
                        });
                }
                clusters
                    .entry(*current_cluster.iter().take(1).next().unwrap())
                    .insert_entry(current_cluster);
            }
        });
        let cluster_j = self.snapshot.get_property_index("cluster");
        let id_j = self.snapshot.get_property_index("id");
        (0..self.snapshot.atoms_count).for_each(|atom_i| {
            let (cluster_i, _) = clusters
                .iter()
                .find(|(_, cluster)| cluster.contains(&atom_i))
                .unwrap();
            self.snapshot.set_atom_value(
                cluster_j,
                atom_i,
                self.snapshot.get_atom_value(id_j, *cluster_i),
            );
        });
        DumpFile::new(vec![self.snapshot])
    }
}

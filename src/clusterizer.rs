use super::{check_cutoff, DumpFile, DumpSnapshot};
use std::collections::{HashMap, HashSet};

pub struct Clusterizer {
    snapshot: DumpSnapshot,
    x_j: usize,
    y_j: usize,
    z_j: usize,
}

impl Clusterizer {
    pub fn new(snapshot_input: &DumpSnapshot) -> Self {
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
        }
    }

    fn get_atom_xyz(&self, atom_i: usize) -> (f64, f64, f64) {
        (
            self.snapshot.get_atom_value(self.x_j, atom_i),
            self.snapshot.get_atom_value(self.y_j, atom_i),
            self.snapshot.get_atom_value(self.z_j, atom_i),
        )
    }

    fn find_cluster_i(&self, cluster: &HashSet<usize>, atoms: &HashSet<usize>) -> Option<usize> {
        cluster
            .iter()
            .copied()
            .filter_map(|cluster_i| {
                let cluster_xyz = self.get_atom_xyz(cluster_i);
                atoms.iter().copied().find(|atom_i| {
                    let atom_xyz = self.get_atom_xyz(*atom_i);
                    check_cutoff(atom_xyz, cluster_xyz, 3.0)
                })
            })
            .next()
    }

    pub fn clusterize(mut self) -> DumpFile {
        let mut indices: HashSet<usize> = (0..self.snapshot.atoms_count).collect();
        let clusters: HashMap<usize, HashSet<usize>> =
            (0..self.snapshot.atoms_count).fold(HashMap::new(), |mut clusters, _| {
                let (cluster_i, atom_i) = clusters
                    .iter()
                    .filter_map(|(cluster_i, cluster)| {
                        self.find_cluster_i(cluster, &indices)
                            .map(|atom_i| (*cluster_i, atom_i))
                    })
                    .next()
                    .unwrap_or({
                        let atom_i = *indices.iter().next().unwrap();
                        (atom_i, atom_i)
                    });
                clusters
                    .entry(cluster_i)
                    .and_modify(|cluster| {
                        cluster.insert(atom_i);
                    })
                    .or_insert(HashSet::from([atom_i]));
                indices.remove(&atom_i);
                clusters
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

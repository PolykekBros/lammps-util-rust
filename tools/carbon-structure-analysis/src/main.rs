use anyhow::Result;
use clap::Parser;
use lammps_files::Dump;
use lammps_util::{xyz::xyz_vec_from_snapshot, MainWrapper, Task, XYZ};
use log::info;
use std::{
    collections::{HashMap, VecDeque},
    path::PathBuf,
};

const OUTER_CUTOFF: f64 = 1.75;
const INNER_CUTOFF: f64 = 1.5;

/// Analyze properties of carbon structres in a .dump file
#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Dump file
    dump_file: PathBuf,

    /// Carbon atom type Id
    #[arg(short, long)]
    carbon_id: usize,
}

#[derive(Clone, Debug)]
struct CarbonAtomsExtractor {
    dump_path: PathBuf,
    carbon_type_id: usize,
}

impl CarbonAtomsExtractor {
    fn new(dump_path: PathBuf, carbon_type_id: usize) -> Self {
        Self {
            dump_path,
            carbon_type_id,
        }
    }

    fn carbon_atoms(self) -> Result<Vec<XYZ>> {
        let dump = Dump::open(&self.dump_path, &[])?;
        let snapshot = &dump.get_snapshots()[0];
        let types = snapshot.get_property("type");
        let mut coords = xyz_vec_from_snapshot(snapshot)
            .into_iter()
            .filter(|xyz| (types[xyz.index] as usize).eq(&self.carbon_type_id))
            .collect();
        XYZ::get_supercell_coords(&mut coords, snapshot.get_symbox(), OUTER_CUTOFF);
        Ok(coords.into())
    }
}

struct RingsFinder {
    adjacency_list: Vec<Vec<usize>>,
    atoms: Vec<XYZ>,
    visited: Vec<bool>,
}

impl RingsFinder {
    pub fn new(atoms: Vec<XYZ>) -> Self {
        let tree = kd_tree::KdTree::build_by_ordered_float(atoms.clone());
        let atoms = tree.items().to_vec();
        let index_map = atoms
            .iter()
            .enumerate()
            .filter(|(_, atom)| !atom.is_ghost)
            .map(|(idx, atom)| (atom.index, idx))
            .collect::<HashMap<_, _>>();
        let adjacency_list = atoms
            .iter()
            .filter(|atom| !atom.is_ghost)
            .map(|atom| {
                let inner = tree.within_radius(atom, INNER_CUTOFF);
                let outer = tree.within_radius(atom, OUTER_CUTOFF);
                outer
                    .into_iter()
                    .filter(|atom| !inner.contains(atom))
                    .map(|atom| index_map[&atom.index])
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let visited = vec![false; atoms.len()];
        info!("Adjacency list ({}) built", adjacency_list.len());
        Self {
            adjacency_list,
            atoms,
            visited,
        }
    }

    pub fn find(mut self) -> HashMap<usize, usize> {
        let mut rings_cnt = HashMap::new();
        for (start, neighbours) in self.adjacency_list.iter().enumerate() {
            self.visited[start] = true;
            let mut local_visisted = self.visited.clone();
            for cnt in neighbours
                .iter()
                .copied()
                .filter(|end| !self.visited[*end])
                .filter_map(|end| self.bfs(start, end, &mut local_visisted))
            {
                *rings_cnt.entry(cnt).or_insert(0) += 1;
            }
        }
        rings_cnt
    }

    fn bfs(&self, start: usize, end: usize, local_visited: &mut Vec<bool>) -> Option<usize> {
        assert_ne!(start, end);
        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.extend(
            self.adjacency_list[start]
                .iter()
                .filter(|&&idx| !self.visited[idx] && idx.ne(&end)),
        );
        let mut count = 1;
        while let Some(current) = queue.pop_front() {
            count += 1;
            if current == end {
                assert_ne!(count, 2);
                return Some(count);
            }
            local_visited[current] = true;
            for neighbor in self.adjacency_list[current]
                .iter()
                .filter(|&&idx| !local_visited[idx])
            {
                queue.push_front(*neighbor);
            }
        }
        None
    }
}

struct RingsFinderTask {
    carbon_atoms_extractor: CarbonAtomsExtractor,
}

impl Task for RingsFinderTask {
    type Output = ();
    fn run(&self) -> Result<Self::Output> {
        let carbon_atoms = self.carbon_atoms_extractor.clone().carbon_atoms()?;
        let rings = RingsFinder::new(carbon_atoms).find();
        println!("Rings: {rings:?}");
        Ok(())
    }
}

impl RingsFinderTask {
    fn new(carbon_atoms_extractor: CarbonAtomsExtractor) -> Self {
        Self {
            carbon_atoms_extractor,
        }
    }
}

fn main() -> Result<()> {
    MainWrapper::<Cli>::default().run(|cli| {
        Box::new(RingsFinderTask::new(CarbonAtomsExtractor::new(
            cli.dump_file,
            cli.carbon_id,
        )))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use lammps_util::geomutil_util::Point3;

    #[test]
    fn test_triangle_graph() {
        let node0 = XYZ::new(Point3::from([0.0, 0.0, 0.0]), 0, false);
        let node1 = XYZ::new(Point3::from([1.6, 0.0, 0.0]), 1, false);
        let node2 = XYZ::new(Point3::from([0.8, 1.3856, 0.0]), 2, false);
        let nodes = vec![node0, node1, node2];
        let finder = RingsFinder::new(nodes);
        assert_eq!(finder.adjacency_list.len(), 3);
        assert_eq!(finder.atoms.len(), 3);
        assert_eq!(finder.visited.len(), 3);
        for neighbours in finder.adjacency_list.iter() {
            assert_eq!(neighbours.len(), 2);
        }
        let rings = finder.find();
        let expected_rings: HashMap<usize, usize> = HashMap::from([(3, 1)]);
        assert_eq!(rings, expected_rings);
    }

    #[test]
    fn test_two_triangles_graph() {
        // TODO: fix this test
        let node0 = XYZ::new(Point3::from([0.0, 0.0, 0.0]), 0, false);
        let node1 = XYZ::new(Point3::from([1.6, 0.0, 0.0]), 1, false);
        let node2 = XYZ::new(Point3::from([0.8, 1.3856, 0.0]), 2, false);
        let node3 = XYZ::new(Point3::from([0.8, -1.3856, 0.0]), 3, false);
        let nodes = vec![node0, node1, node2, node3];
        let finder = RingsFinder::new(nodes);
        let mut neigh_cnt = HashMap::from([(3, 2), (2, 2)]);
        for neighs in finder.adjacency_list.iter() {
            let cnt = neighs.len();
            assert!(neigh_cnt.contains_key(&cnt));
            neigh_cnt.entry(cnt).and_modify(|cnt| *cnt -= 1);
        }
        let neigh_cnt_expected = HashMap::from([(3, 0), (2, 0)]);
        assert_eq!(neigh_cnt_expected, neigh_cnt);
        let rings = finder.find();
        let expected_rings: HashMap<usize, usize> = HashMap::from([(3, 2)]);
        assert_eq!(rings, expected_rings);
    }
}

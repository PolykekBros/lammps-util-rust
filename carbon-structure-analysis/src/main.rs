use anyhow::{Context, Result};
use clap::Parser;
use lammps_util_rust::{DumpFile, DumpSnapshot, XYZ};
use log::{debug, info};
use std::collections::{HashMap, HashSet, VecDeque};
use std::path::{Path, PathBuf};

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

fn load_snapshot(path: &Path) -> Result<DumpSnapshot> {
    let dump = DumpFile::read(path, &[]).context(format!(
        "Failed to read .dump file: {}",
        path.to_string_lossy()
    ))?;
    let snapshot = dump.get_snapshots()[0];
    Ok(snapshot.to_owned())
}

fn get_carbon_atoms(snapshot: &DumpSnapshot, type_id: usize) -> Vec<XYZ> {
    let coords = snapshot.get_coordinates();
    let types = snapshot.get_property("type");
    coords
        .iter()
        .filter(|xyz| (types[xyz.index()] as usize).eq(&type_id))
        .copied()
        .collect()
}

#[derive(Debug, Clone)]
struct Neighbours {
    atom: XYZ,
    neighbours: Vec<XYZ>,
}

struct RingsFinder {
    adjecency_list: HashMap<XYZ, Vec<XYZ>>,
}

impl RingsFinder {
    pub fn new(atoms: Vec<XYZ>) -> Self {
        let mut adjecency_list = HashMap::new();
        let tree = kd_tree::KdTree::build_by_ordered_float(atoms.clone());
        for atom in atoms {
            let inner = tree.within_radius(&atom, 1.5);
            let outer = tree.within_radius(&atom, 1.7);
            let neighbours = outer
                .iter()
                .filter(|atom| !inner.contains(atom))
                .map(|atom| **atom)
                .collect::<Vec<_>>();
            adjecency_list.entry(atom).insert_entry(neighbours);
        }
        Self { adjecency_list }
    }

    pub fn find(&self) -> HashMap<usize, usize> {
        let mut rings = Vec::new();
        for (&atom_begin, neighbours) in self.adjecency_list.iter() {
            for &atom_end in neighbours {
                if let Some(ring) = self.bfs(atom_begin, atom_end) {
                    if !rings.contains(&ring) {
                        rings.push(ring);
                    }
                }
            }
        }
        let mut rings_cnt = HashMap::new();
        for ring in rings {
            rings_cnt
                .entry(ring.len())
                .and_modify(|c| *c += 1)
                .or_insert(1);
        }
        rings_cnt
    }

    fn bfs(&self, atom_begin: XYZ, atom_end: XYZ) -> Option<HashSet<XYZ>> {
        let mut queue = VecDeque::from([(atom_begin, HashSet::from([atom_begin]))]);
        let mut visited = HashSet::new();
        while let Some((atom, ring)) = queue.pop_front() {
            if visited.contains(&atom) {
                continue;
            }
            visited.insert(atom);
            for &atom_next in self.adjecency_list[&atom].iter() {
                let mut ring = ring.clone();
                if atom == atom_begin && atom_next == atom_end {
                    continue;
                }
                ring.insert(atom_next);
                if atom_next == atom_end {
                    return Some(ring);
                }
                queue.push_back((atom_next, ring));
            }
        }
        None
    }
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let snapshot = load_snapshot(&cli.dump_file)?;
    let atoms = get_carbon_atoms(&snapshot, cli.carbon_id);
    info!("Loaded {} carbon atoms", atoms.len());
    let rings = RingsFinder::new(atoms).find();
    println!("Rings: {rings:?}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_triangle_graph() {
        let node0 = XYZ::from([0.0, 0.0, 0.0], 0);
        let node1 = XYZ::from([1.6, 0.0, 0.0], 0);
        let node2 = XYZ::from([0.8, 1.3856, 0.0], 0);

        let nodes = vec![node0, node1, node2];
        let finder = RingsFinder::new(nodes);

        for (_, neighbours) in finder.adjecency_list.iter() {
            assert_eq!(neighbours.len(), 2);
        }

        let rings = finder.find();
        let expected_rings: HashMap<usize, usize> = HashMap::from([(3, 1)]);

        assert_eq!(rings, expected_rings);
    }

    #[test]
    fn test_two_triangles_graph() {
        let node0 = XYZ::from([0.0, 0.0, 0.0], 0);
        let node1 = XYZ::from([1.6, 0.0, 0.0], 0);
        let node2 = XYZ::from([0.8, 1.3856, 0.0], 0);
        let node3 = XYZ::from([0.8, -1.3856, 0.0], 0);

        let nodes = vec![node0, node1, node2, node3];
        let finder = RingsFinder::new(nodes);

        assert_eq!(finder.adjecency_list[&node0].len(), 3);
        assert_eq!(finder.adjecency_list[&node1].len(), 3);
        assert_eq!(finder.adjecency_list[&node2].len(), 2);
        assert_eq!(finder.adjecency_list[&node3].len(), 2);

        let rings = finder.find();
        let expected_rings: HashMap<usize, usize> = HashMap::from([(3, 2)]);

        assert_eq!(rings, expected_rings);
    }
}

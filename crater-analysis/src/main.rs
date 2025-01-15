use anyhow::Result;
use clap::Parser;
use lammps_util_rust::{clusterize_snapshot, copy_snapshot_with_indices, DumpFile, DumpSnapshot};
use std::{collections::HashMap, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// DUMP INPUT
    #[arg()]
    dump_input_file: PathBuf,

    /// DUMP FINAL
    #[arg()]
    dump_final_file: PathBuf,

    /// OUTPUT DIR
    #[arg()]
    output_dir: PathBuf,

    /// Max depth (A)
    #[arg(short, long, default_value_t = 50.0)]
    max_depth: f64,

    /// Cutoff (A)
    #[arg(short, long, default_value_t = 3.0)]
    cutoff: f64,

    /// Stripe width (A)
    #[arg(short, long, default_value_t = 5.43 / 2.0)]
    width: f64,
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
    copy_snapshot_with_indices(initial_snapshot, indices.into_iter())
}

fn get_crater_info(snapshot: &DumpSnapshot, zero_lvl: f64) -> String {
    let cluster = snapshot.get_property("cluster");
    let cluster_cnt = cluster
        .iter()
        .copied()
        .fold(HashMap::new(), |mut indices, cluster| {
            indices
                .entry(cluster as u64)
                .and_modify(|cnt| *cnt += 1)
                .or_insert(1);
            indices
        });
    let (cluster_max, _) = cluster_cnt.iter().fold(
        (0, u64::MIN),
        |(cluster_max, max_cnt), (cluster, cnt)| match *cnt > max_cnt {
            true => (*cluster, *cnt),
            false => (cluster_max, max_cnt),
        },
    );
    let z = snapshot.get_property("z");
    let mut crater_count = 0;
    let mut surface_count = 0;
    let mut z_avg = 0.0;
    let mut z_min = f64::INFINITY;
    for i in 0..snapshot.atoms_count {
        let z = z[i];
        let cluster = cluster[i];
        if cluster as u64 == cluster_max {
            if z > -2.4 * 0.707 + zero_lvl {
                surface_count += 1;
            }
            crater_count += 1;
            z_min = z_min.min(z - zero_lvl);
            z_avg += z - zero_lvl;
        }
    }
    z_avg /= crater_count as f64;
    let volume = crater_count as f64 * 20.1;
    let surface = surface_count as f64 * 7.3712;
    format!("{crater_count} {volume} {surface} {z_avg} {z_min}")
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let dump_input = DumpFile::read(&cli.dump_input_file, &[])?;
    let snapshot_input = dump_input.get_snapshots()[0];
    let zero_lvl = snapshot_input.get_zero_lvl();

    let dump_final = DumpFile::read(&cli.dump_final_file, &[])?;
    let snapshot_final = dump_final.get_snapshots()[0];

    let snapshot_candidates =
        crater_candidates_snapshot(snapshot_input, snapshot_final, cli.cutoff);
    let snapshot_candidates_clusterized = clusterize_snapshot(&snapshot_candidates, 3.0);
    let info = get_crater_info(&snapshot_candidates_clusterized, zero_lvl);

    let dump_candidates_clusterized = DumpFile::new(vec![snapshot_candidates_clusterized]);
    dump_candidates_clusterized.save(&cli.output_dir.join("dump.crater_candidates"))?;
    println!("{info}");
    Ok(())
}

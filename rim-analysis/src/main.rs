use anyhow::{Context, Result};
use clap::Parser;
use lammps_util_rust::{
    clusterize_snapshot, copy_snapshot_with_indices, get_cluster_counts, DumpFile, DumpSnapshot,
};
use log::info;
use nalgebra::{Vector1, Vector2};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

const RIM_THRESHOLD: usize = 50;
const CLUSTER_TRAJECTORY_LINE: usize = 38;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// DUMP INPUT
    dump_input_file: PathBuf,

    /// DUMP FINAL
    dump_final_file: PathBuf,

    /// OUTPUT DIR
    output_dir: PathBuf,

    /// Cutoff (A)
    #[arg(short, long, default_value_t = 3.0)]
    cutoff: f64,
}

fn get_above_zero(snap: &DumpSnapshot, zero_lvl: f64) -> DumpSnapshot {
    let indices = snap
        .get_property("z")
        .iter()
        .enumerate()
        .filter(|(_, &z)| z > zero_lvl)
        .map(|(i, _)| i);
    copy_snapshot_with_indices(snap, indices)
}

fn filter_clusters(snap: &DumpSnapshot) -> HashSet<usize> {
    let cnts = get_cluster_counts(&snap);
    cnts.iter()
        .filter(|(_, &cnt)| cnt >= RIM_THRESHOLD)
        .map(|(&id, _)| id)
        .collect()
}

fn rim_snapshot(
    initial_snapshot: &DumpSnapshot,
    final_snapshot: &DumpSnapshot,
    cutoff: f64,
) -> DumpSnapshot {
    let zero_lvl = initial_snapshot.get_zero_lvl();
    let above_zero_lvl = get_above_zero(final_snapshot, zero_lvl);
    let clusters = clusterize_snapshot(&above_zero_lvl, cutoff);
    let clusters_selected = filter_clusters(&clusters);
    let indices = clusters
        .get_property("cluster")
        .iter()
        .copied()
        .enumerate()
        .filter(|(_, cluster)| clusters_selected.contains(&(*cluster as usize)))
        .map(|(i, _)| i);
    copy_snapshot_with_indices(&clusters, indices)
}

fn get_center_pos(cluster_xyz_path: &Path) -> Result<(f64, f64)> {
    let file = File::open(cluster_xyz_path).with_context(|| {
        format!(
            "Failed to read cluster trajectory from {}",
            cluster_xyz_path.display()
        )
    })?;
    let reader = BufReader::new(file);
    let line =
        reader
            .lines()
            .nth(CLUSTER_TRAJECTORY_LINE - 1)
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Failed to read line {CLUSTER_TRAJECTORY_LINE} from cluster trajectory",
                )
            })??;
    let mut tokens = line.split_whitespace();
    tokens.next();
    let x = tokens
        .next()
        .ok_or_else(|| anyhow::anyhow!("Missing X coordinate in line {CLUSTER_TRAJECTORY_LINE}"))?
        .parse::<f64>()
        .with_context(|| {
            format!("Failed to parse X coordinate in line {CLUSTER_TRAJECTORY_LINE}",)
        })?;
    let y = tokens
        .next()
        .ok_or_else(|| anyhow::anyhow!("Missing Y coordinate in line {CLUSTER_TRAJECTORY_LINE}"))?
        .parse::<f64>()
        .with_context(|| {
            format!("Failed to parse Y coordinate in line {CLUSTER_TRAJECTORY_LINE}",)
        })?;
    Ok((x, y))
}

fn get_angle_distribution(snap: &DumpSnapshot, center_x: f64, center_y: f64) -> Vec<usize> {
    let center = vector![center_x, center_y];
    let coords = snap.get_coordinates().into_iter().map(|xyz| Vector2)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let dump_input = DumpFile::read(&cli.dump_input_file, &[])?;
    let snapshot_input = dump_input.get_snapshots()[0];

    let dump_final = DumpFile::read(&cli.dump_final_file, &[])?;
    let snapshot_final = dump_final.get_snapshots()[0];

    let (center_x, center_y) = get_center_pos(&cli.output_dir.join("cluster_xyz.txt"))
        .context("Failed to get center position")?;
    info!("cluster center: ({center_x}, {center_y})");
    let rim = rim_snapshot(snapshot_input, snapshot_final, cli.cutoff);

    info!("rim count: {}", rim.atoms_count);

    let dump_rim = DumpFile::new(vec![rim]);
    dump_rim.save(&cli.output_dir.join("dump.rim"))?;
    Ok(())
}

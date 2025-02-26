use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use lammps_util_rust::{
    clusterize_snapshot, copy_snapshot_with_indices, get_cluster_counts, DumpFile, DumpSnapshot,
};
use log::info;
use nalgebra::{vector, Vector2};
use std::{
    array::from_fn,
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    str::FromStr,
};

const RIM_THRESHOLD: usize = 10;
const CLUSTER_TRAJECTORY_LINE: usize = 38;
const ANGLE_ROTATION: usize = 10;
const DATA_LEN: usize = 360 / ANGLE_ROTATION;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Cutoff (A)
    #[arg(short, long, default_value_t = 3.0)]
    cutoff: f64,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Analyse a single run dir
    Single {
        /// Run directory
        dir: PathBuf,
    },

    /// Analyse a whole results dir
    All {
        /// Results directory
        dir: PathBuf,
    },
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

fn get_center_pos(cluster_xyz_path: &Path) -> Result<Vector2<f64>> {
    let reader = File::open(cluster_xyz_path)
        .map(BufReader::new)
        .with_context(|| {
            format!(
                "Failed to read cluster trajectory from {}",
                cluster_xyz_path.display()
            )
        })?;
    let Some(Ok(line)) = reader.lines().nth(CLUSTER_TRAJECTORY_LINE - 1) else {
        bail!("Failed to read line {CLUSTER_TRAJECTORY_LINE} from cluster trajectory");
    };
    let mut iter = line.split_whitespace().skip(1).map(f64::from_str);
    let Some((Ok(x), Ok(y))) = iter.next().zip(iter.next()) else {
        bail!("Failed to parse XY coordinates in the line {CLUSTER_TRAJECTORY_LINE}");
    };
    Ok(vector![x, y])
}

fn angle_between_vectors(a: Vector2<f64>, b: Vector2<f64>) -> f64 {
    let angle = a.y.atan2(a.x) - b.y.atan2(b.x);
    angle.to_degrees() + 180.0
}

fn get_angle_distribution(snap: &DumpSnapshot, center: Vector2<f64>) -> [Vec<f64>; DATA_LEN] {
    let start = vector![0.0, -1.0];
    let mut radii = from_fn(|_| Vec::new());
    for coord in snap
        .get_coordinates()
        .into_iter()
        .map(|xyz| vector![xyz.x, xyz.y] - center)
    {
        let angle = angle_between_vectors(start, coord);
        let index = (angle / ANGLE_ROTATION as f64) as usize;
        let radius = coord.magnitude();
        radii[index].push(radius);
    }
    radii
}

fn get_mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn get_std(values: &[f64]) -> f64 {
    let mean = get_mean(values);
    let std = values.iter().map(|r| (r - mean) * (r - mean)).sum::<f64>() / values.len() as f64;
    std.sqrt()
}

fn get_data(radii: [Vec<f64>; DATA_LEN]) -> [[f64; 5]; DATA_LEN] {
    let mut data = from_fn(|_| [0.0; 5]);
    for (i, values) in radii.iter().enumerate() {
        let mean = get_mean(&values);
        let cnt = values.len() as f64;
        data[i] = [
            (i * ANGLE_ROTATION + ANGLE_ROTATION / 2) as f64,
            cnt,
            mean,
            cnt * mean,
            get_std(&values),
        ]
    }
    data
}

fn parse_run_dir(dir: &Path, cutoff: f64) -> Result<[Vec<f64>; DATA_LEN]> {
    let dump_input = DumpFile::read(&dir.join("dump.initial"), &[])?;
    let snap_input = dump_input.get_snapshots()[0];
    let dump_final = DumpFile::read(&dir.join("dump.final_no_cluster"), &[])?;
    let snap_final = dump_final.get_snapshots()[0];
    let center =
        get_center_pos(&dir.join("cluster_xyz.txt")).context("Failed to get center position")?;
    info!("cluster center: {center}");
    let rim = rim_snapshot(snap_input, snap_final, cutoff);
    info!("rim count: {}", rim.atoms_count);
    let radii = get_angle_distribution(&rim, center);
    let dump_rim = DumpFile::new(vec![rim]);
    dump_rim.save(&dir.join("dump.rim"))?;
    Ok(radii)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let radii = match cli.command {
        Commands::Single { dir } => parse_run_dir(&dir, cli.cutoff)?,
        Commands::All { dir } => parse_run_dir(&dir, cli.cutoff)?,
    };
    let data = get_data(radii);
    info!("data: {data:?}");
    Ok(())
}

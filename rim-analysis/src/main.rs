use anyhow::{anyhow, bail, Context, Result};
use clap::Parser;
use lammps_util_rust::{
    clusterize_snapshot, copy_snapshot_with_indices, get_cluster_counts, DumpFile, DumpSnapshot,
};
use log::info;
use nalgebra::{vector, Vector2};
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    str::FromStr,
};

const RIM_THRESHOLD: usize = 50;
const CLUSTER_TRAJECTORY_LINE: usize = 38;
const ANGLE_ROTATION: usize = 10;

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
    vector![x, y];
}

fn angle_between_vectors(a: Vector2<f64>, b: Vector2<f64>) -> f64 {
    let dot = a.dot(&b);
    let a_len = a.norm();
    let b_len = b.norm();
    let cos_theta = dot / (a_len * b_len);
    let cos_theta = cos_theta.clamp(-1.0, 1.0);
    let angle_radians = cos_theta.acos();
    let angle_degrees = angle_radians * 180.0 / PI;
    angle_degrees
}

fn get_angle_distribution(
    snap: &DumpSnapshot,
    center: Vector2<f64>,
) -> (Vec<Vec<f64>>, Vec<usize>) {
    let start = vector![0, -1];
    let mut radii = vec![Vec::new(); 360 / ANGLE_ROTATION];
    let mut count = vec![0; 360 / ANGLE_ROTATION];
    for coord in snap
        .get_coordinates()
        .into_iter()
        .map(|xyz| vector![xyz.x, xyz.y] - center)
    {
        let angle = start.angle(&coord).to_degrees();
        let angle = (angle + 360.0) % 360.0;
        let index = (angle / ANGLE_ROTATION as f64) as usize;
        let radius = coord.magnitude();
        count[index] += 1;
        radii[index].push(radius);
    }
    (radii, count)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let dump_input = DumpFile::read(&cli.dump_input_file, &[])?;
    let snapshot_input = dump_input.get_snapshots()[0];

    let dump_final = DumpFile::read(&cli.dump_final_file, &[])?;
    let snapshot_final = dump_final.get_snapshots()[0];

    let center = get_center_pos(&cli.output_dir.join("cluster_xyz.txt"))
        .context("Failed to get center position")?;
    info!("cluster center: {center}");
    let rim = rim_snapshot(snapshot_input, snapshot_final, cli.cutoff);

    info!("rim count: {}", rim.atoms_count);

    let dump_rim = DumpFile::new(vec![rim]);
    dump_rim.save(&cli.output_dir.join("dump.rim"))?;
    Ok(())
}

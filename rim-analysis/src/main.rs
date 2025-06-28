use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand};
use core::f32;
use geomutil_util::Point2;
use itertools::{izip, Itertools};
use lammps_util_rust::{
    clusterize_snapshot, copy_snapshot_with_indices, get_cluster_counts, process_results_dir,
    DumpFile, DumpSnapshot, IteratorAvg,
};
use log::info;
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader},
    iter::zip,
    path::{Path, PathBuf},
    str::FromStr,
};

const RIM_THRESHOLD: usize = 10;
const CLUSTER_TRAJECTORY_LINE: usize = 38;
const ANGLE_ROTATION: usize = 10;
const SECTORS_LEN: usize = 360 / ANGLE_ROTATION;
const MASS: [f32; 3] = [0.0, 28.08553, 12.011];
type Sectors = [Sector; SECTORS_LEN];

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
    /// Crater analysis for a single run dir
    Single(SingleCMD),

    /// Crater analysis for the whole results folder
    Multi(MultiCMD),
}

#[derive(Args)]
struct SingleCMD {
    run_dir: PathBuf,
}

#[derive(Args)]
struct MultiCMD {
    results_dir: PathBuf,

    /// Number of threads to run in parallel
    #[arg(short, long, default_value_t = 2)]
    threads: usize,
}

#[derive(Debug, Clone, Copy)]
struct Atom {
    coords: Point2,
    atype: usize,
}

impl Atom {
    fn new(coords: Point2, atype: usize) -> Self {
        Self { coords, atype }
    }
}

#[derive(Debug, Clone)]
struct RimValues {
    atoms: Vec<Atom>,
    center: Point2,
}

impl RimValues {
    fn new(atoms: Vec<Atom>, center: Point2) -> Self {
        let mut atoms = atoms;
        atoms.iter_mut().for_each(|a| a.coords -= center);
        Self { atoms, center }
    }

    fn get_sectors(&self) -> Sectors {
        std::array::from_fn(|i| {
            Sector::new(
                self.atoms
                    .iter()
                    .filter(|a| {
                        ((point_rotation(a.coords) / ANGLE_ROTATION as f32).floor() as usize)
                            .clamp(0, SECTORS_LEN)
                            == i
                    })
                    .copied()
                    .collect(),
            )
        })
    }
}

#[derive(Debug, Clone, Default)]
struct Sector {
    radius: Vec<f32>,
    mass: Vec<f32>,
    count: Vec<usize>,
}

impl Sector {
    fn new(atoms: Vec<Atom>) -> Self {
        Self {
            radius: atoms.iter().map(|a| a.coords.length()).collect(),
            mass: atoms.iter().map(|a| MASS[a.atype]).collect(),
            count: vec![atoms.len()],
        }
    }
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

fn get_rim_snapshot(
    initial_snapshot: &DumpSnapshot,
    final_snapshot: &DumpSnapshot,
    cutoff: f64,
) -> DumpSnapshot {
    let zero_lvl = initial_snapshot.get_zero_lvl();
    let above_zero_lvl = get_above_zero(final_snapshot, zero_lvl);
    let clusters = clusterize_snapshot(&above_zero_lvl, cutoff);
    let clusters_selected = get_cluster_counts(&clusters)
        .iter()
        .filter(|(_, &cnt)| cnt >= RIM_THRESHOLD)
        .map(|(&id, _)| id)
        .collect::<HashSet<_>>();
    let indices = clusters
        .get_property("cluster")
        .iter()
        .copied()
        .enumerate()
        .filter(|(_, cluster)| clusters_selected.contains(&(*cluster as usize)))
        .map(|(i, _)| i);
    copy_snapshot_with_indices(&clusters, indices)
}

fn get_rim_atoms(snap: &DumpSnapshot) -> Vec<Atom> {
    izip!(
        snap.get_property("x").iter().copied(),
        snap.get_property("y").iter().copied(),
        snap.get_property("type").iter().copied(),
    )
    .map(|(x, y, atype)| Atom::new(Point2::from([x as f32, y as f32]), atype as usize))
    .collect()
}

fn get_center_pos(path: &Path) -> Result<Point2> {
    let reader = File::open(path)
        .map(BufReader::new)
        .with_context(|| format!("Failed to read cluster trajectory from {}", path.display()))?;
    let Some(Ok(line)) = reader.lines().nth(CLUSTER_TRAJECTORY_LINE - 1) else {
        bail!("Failed to read line {CLUSTER_TRAJECTORY_LINE} from cluster trajectory");
    };
    let mut iter = line.split_whitespace().skip(1).map(f64::from_str);
    let Some((Ok(x), Ok(y))) = iter.next().zip(iter.next()) else {
        bail!("Failed to parse XY coordinates in the line {CLUSTER_TRAJECTORY_LINE}");
    };
    Ok(Point2::from([x as f32, y as f32]))
}

fn get_rim_values(dir: &Path, cutoff: f64) -> Result<RimValues> {
    let dump_input = DumpFile::read(&dir.join("dump.initial"), &[])?;
    let snap_input = dump_input.get_snapshots()[0];
    let dump_final = DumpFile::read(&dir.join("dump.final_no_cluster"), &[])?;
    let snap_final = dump_final.get_snapshots()[0];
    let snap_rim = get_rim_snapshot(snap_input, snap_final, cutoff);
    let atoms = get_rim_atoms(&snap_rim);
    let dump_rim = DumpFile::new(vec![snap_rim]);
    dump_rim.save(&dir.join("dump.rim"))?;
    let center =
        get_center_pos(&dir.join("cluster_xyz.txt")).context("Failed to get center position")?;
    info!("cluster center: {center:?}");
    Ok(RimValues::new(atoms, center))
}

fn rim_atom_rotation(v: Point2) -> f32 {
    let angle = (-v.y / v.length()).acos();
    if v.x > 0.0 {
        angle
    } else {
        2.0 * f32::consts::PI - angle
    }
}

fn point_rotation(v: Point2) -> f32 {
    rim_atom_rotation(v).to_degrees()
}

fn parse_run_dir(dir: &Path, cutoff: f64) -> Result<Sectors> {
    let rim_values = get_rim_values(dir, cutoff)?;
    info!("rim count: {}", rim_values.atoms.len());
    Ok(rim_values.get_sectors())
}

fn run_single(dir: &Path, cutoff: f64) -> Result<Sectors> {
    parse_run_dir(dir, cutoff)
}

fn run_multi(dir: &Path, threads: usize, cutoff: f64) -> Result<Sectors> {
    Ok(
        process_results_dir(dir, threads, |dir| parse_run_dir(&dir.path, cutoff))?
            .into_iter()
            .map(|(_, sectors)| sectors)
            .reduce(|mut acc, sectors| {
                zip(acc.iter_mut(), sectors).for_each(|(a, b)| {
                    a.mass.extend(b.mass);
                    a.radius.extend(b.radius);
                    a.count.extend(b.count);
                });
                acc
            })
            .unwrap(),
    )
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let values = match cli.command {
        Commands::Single(args) => run_single(&args.run_dir, cli.cutoff),
        Commands::Multi(args) => run_multi(&args.results_dir, args.threads, cli.cutoff),
    }?;
    let table = values
        .into_iter()
        .enumerate()
        .map(|(i, row)| {
            let (n_avg, n_std) = row
                .count
                .into_iter()
                .map(|x| x as f64)
                .avg_with_std()
                .unwrap();
            let (m_avg, m_std) = row.mass.into_iter().avg_with_std().unwrap();
            let (r_avg, r_std) = row.radius.into_iter().avg_with_std().unwrap();
            let vals = [n_avg, n_std, m_avg, m_std, r_avg, r_std]
                .into_iter()
                .map(|x| format!("{x:10.4}"))
                .join("\t");
            format!("{}\t{}", i * ANGLE_ROTATION, vals)
        })
        .join("\n");
    println!("# φ N σ(N) m σ(m) r σ(r)\n{table}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use assert_float_eq::assert_float_absolute_eq;

    use super::*; // Import the `add` function from the parent module

    #[test]
    fn rim_atom_rotation_test() {
        assert_float_absolute_eq!(point_rotation(Point2::from([1.0, -1.0])), 45.0, 1e-4);
        assert_float_absolute_eq!(point_rotation(Point2::from([1.0, 0.0])), 90.0, 1e-4);
        assert_float_absolute_eq!(point_rotation(Point2::from([1.0, 1.0])), 135.0, 1e-4);
        assert_float_absolute_eq!(point_rotation(Point2::from([0.0, 1.0])), 180.0, 1e-4);
        assert_float_absolute_eq!(point_rotation(Point2::from([-1.0, 1.0])), 225.0, 1e-4);
        assert_float_absolute_eq!(point_rotation(Point2::from([-1.0, 0.0])), 270.0, 1e-4);
        assert_float_absolute_eq!(point_rotation(Point2::from([-1.0, -1.0])), 315.0, 1e-4);
        assert_float_absolute_eq!(point_rotation(Point2::from([0.0, -1.0])), 360.0, 1e-4);
    }
}

use std::{iter, path::PathBuf};

use anyhow::Result;
use clap::Parser;
use lammps_util_rust::{crater_snapshot, geomutil_util::Point3, DumpFile, DumpSnapshot, XYZ};
use log::debug;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Dump initial
    dump_initial: PathBuf,

    /// Dump initial
    dump_final: PathBuf,
}

fn get_coords_shift(a: &[XYZ], b: &[XYZ]) -> (usize, Point3, Point3) {
    let deltas = iter::zip(a, b)
        .map(|(a, b)| **b - **a)
        .collect::<Vec<Point3>>();
    let count = deltas.len();
    let sum = deltas.iter().copied().reduce(|a, b| a + b).unwrap();
    let sum2 = deltas
        .into_iter()
        .map(|a| a * a)
        .reduce(|a, b| a + b)
        .unwrap();
    (count, sum, sum2)
}

fn get_ids(input_snapshot: &DumpSnapshot, final_snapshot: &DumpSnapshot) -> Vec<f64> {
    let crater_snapshot = crater_snapshot(input_snapshot, final_snapshot, 1.75, 3.0);
    let final_ids = final_snapshot.get_property("id");
    let crater_ids = crater_snapshot.get_property("id");
    let ids = crater_ids
        .iter()
        .filter(|id| final_ids.contains(id))
        .copied()
        .collect::<Vec<_>>();
    debug!("ids.len: {}", ids.len());
    ids
}

fn get_coords_filtered(snapshot: &DumpSnapshot, ids: &[f64]) -> Vec<XYZ> {
    iter::zip(snapshot.get_coordinates(), snapshot.get_property("id"))
        .filter(|(_, id)| ids.contains(id))
        .map(|(xyz, _)| xyz)
        .collect()
}

fn get_coords(
    input_snapshot: &DumpSnapshot,
    final_snapshot: &DumpSnapshot,
) -> (Vec<XYZ>, Vec<XYZ>) {
    let ids = get_ids(input_snapshot, final_snapshot);
    (
        get_coords_filtered(input_snapshot, &ids),
        get_coords_filtered(final_snapshot, &ids),
    )
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let dump_initial = DumpFile::read(&cli.dump_initial, &[])?;
    let initial_snapshot = dump_initial.get_snapshots()[0];

    let dump_final = DumpFile::read(&cli.dump_final, &[])?;
    let final_snapshot = dump_final.get_snapshots()[0];

    let (input_coords, final_coords) = get_coords(initial_snapshot, final_snapshot);
    let (cnt, sum, sum2) = get_coords_shift(&input_coords, &final_coords);
    println!("{cnt}");
    println!("{} {} {}", sum.x, sum.y, sum.z);
    println!("{} {} {}", sum2.x, sum2.y, sum2.z);
    Ok(())
}

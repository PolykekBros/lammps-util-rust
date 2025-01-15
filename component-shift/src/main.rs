use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use itertools::izip;
use lammps_util_rust::{DumpFile, DumpSnapshot, XYZ, clusterize_snapshot};
use nalgebra::{Point3, Vector2, Vector3};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Dump initial
    dump_initial: PathBuf,

    /// Dump initial
    dump_final: PathBuf,
}

fn calc_shift(
    a: &[XYZ],
    a_id: &[usize],
    b: &[XYZ],
    b_map: &HashMap<usize, usize>,
) -> (usize, Vector3<f64>, Vector3<f64>) {
    let deltas = izip!(a, a_id)
        .filter(|(_, a_id)| b_map.contains_key(*a_id))
        .map(|(a, a_id)| **a - *b[b_map[a_id]])
        .collect::<Vec<Vector3<f64>>>();
    let sum = deltas.iter().copied().sum::<Vector3<f64>>();
    let sum2 = deltas
        .iter()
        .copied()
        .map(|a| a.component_mul(&a))
        .sum::<Vector3<f64>>();
    (sum.len(), sum, sum2)
}

fn get_properties(input_snapshot: &DumpSnapshot, final_snapshot: &DumpSnapshot) -> (f64, f64, f64) {
    let input_coords = input_snapshot.get_coordinates();
    let input_ids = input_snapshot.get_property("id");
    let final_coords = final_snapshot.get_coordinates();
    let final_ids = final_snapshot.get_property("id");
    let cluster_snapshot
    let f_map = final_ids
        .iter()
        .copied()
        .enumerate()
        .fold(HashMap::new(), |mut map, (i, id)| {
            map.insert(id as usize, i);
            map
        });
    todo!()
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let dump_initial = DumpFile::read(&cli.dump_initial, &[])?;
    let dump_initial_snapshot = dump_initial.get_snapshots()[0];
    let dump_final = DumpFile::read(&cli.dump_final, &[])?;
    let dump_final_snapshot = dump_final.get_snapshots()[0];
    // let (rms_x, rms_y, rms_z) = calc_rms_all_components(dump_initial_snapshot, dump_final_snapshot);
    // println!("{rms_x:.6} {rms_y:.6} {rms_z:.6}");
    Ok(())
}

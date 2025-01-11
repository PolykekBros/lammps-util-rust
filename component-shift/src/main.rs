use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use itertools::izip;
use lammps_util_rust::{DumpFile, DumpSnapshot};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Dump initial
    #[arg()]
    dump_initial: PathBuf,

    /// Dump initial
    #[arg()]
    dump_final: PathBuf,
}

fn get_properties(snapshot: &DumpSnapshot) -> (&[f64], &[f64], &[f64], &[f64]) {
    (
        snapshot.get_property("x"),
        snapshot.get_property("y"),
        snapshot.get_property("z"),
        snapshot.get_property("id"),
    )
}

fn calc_rms_single_component(
    a: &[f64],
    a_id: &[f64],
    b: &[f64],
    b_map: &HashMap<usize, usize>,
) -> f64 {
    let a = a.iter().copied();
    let a_id = a_id.iter().copied();
    let (sum, count) = izip!(a, a_id)
        .filter_map(|(a, a_id)| {
            let a_id = a_id as usize;
            if b_map.contains_key(&a_id) {
                let b = b[b_map[&a_id]];
                Some(a - b)
            } else {
                None
            }
        })
        .fold((0.0, 0.0), |(sum, count), a| (sum + a, count + 1.0));
    (sum / count).sqrt()
}

fn calc_rms_all_components(
    initial_snapshot: &DumpSnapshot,
    final_snapshot: &DumpSnapshot,
) -> (f64, f64, f64) {
    let (i_x, i_y, i_z, i_id) = get_properties(initial_snapshot);
    let (f_x, f_y, f_z, f_id) = get_properties(final_snapshot);
    let f_map = f_id
        .iter()
        .copied()
        .enumerate()
        .fold(HashMap::new(), |mut map, (i, id)| {
            map.insert(id as usize, i);
            map
        });
    (
        calc_rms_single_component(i_x, i_id, f_x, &f_map),
        calc_rms_single_component(i_y, i_id, f_y, &f_map),
        calc_rms_single_component(i_z, i_id, f_z, &f_map),
    )
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let dump_initial = DumpFile::read(&cli.dump_initial, &[])?;
    let dump_initial_snapshot = dump_initial.get_snapshots()[0];
    let dump_final = DumpFile::read(&cli.dump_final, &[])?;
    let dump_final_snapshot = dump_final.get_snapshots()[0];
    let (rms_x, rms_y, rms_z) = calc_rms_all_components(dump_initial_snapshot, dump_final_snapshot);
    println!("{rms_x:.6} {rms_y:.6} {rms_z:.6}");
    Ok(())
}

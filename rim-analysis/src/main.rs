use anyhow::Result;
use clap::Parser;
use lammps_util_rust::{
    clusterize_snapshot, copy_snapshot_with_indices, get_max_cluster, DumpFile, DumpSnapshot,
};
use log::info;
use std::path::PathBuf;

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

fn rim_snapshot(
    initial_snapshot: &DumpSnapshot,
    final_snapshot: &DumpSnapshot,
    cutoff: f64,
) -> DumpSnapshot {
    let zero_lvl = initial_snapshot.get_zero_lvl();
    let above_zero_lvl = get_above_zero(final_snapshot, zero_lvl);
    let clusters = clusterize_snapshot(&above_zero_lvl, cutoff);
    let max_cluster = get_max_cluster(&clusters);
    let indices = clusters
        .get_property("cluster")
        .iter()
        .copied()
        .enumerate()
        .filter(|(_, cluster)| *cluster as usize == max_cluster)
        .map(|(i, _)| i);
    copy_snapshot_with_indices(&clusters, indices)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let dump_input = DumpFile::read(&cli.dump_input_file, &[])?;
    let snapshot_input = dump_input.get_snapshots()[0];

    let dump_final = DumpFile::read(&cli.dump_final_file, &[])?;
    let snapshot_final = dump_final.get_snapshots()[0];

    let rim = rim_snapshot(snapshot_input, snapshot_final, cli.cutoff);
    info!("rim count: {}", rim.atoms_count);

    let dump_rim = DumpFile::new(vec![rim]);
    dump_rim.save(&cli.output_dir.join("dump.rim"))?;
    Ok(())
}

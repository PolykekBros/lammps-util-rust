use anyhow::Result;
use clap::{Args, Parser, Subcommand};
use lammps_util_rust::{
    DumpFile, clusterize_snapshot, copy_snapshot_with_indices, get_cluster_counts,
};
use rayon::{ThreadPoolBuilder, prelude::*};
use regex::Regex;
use std::{
    fs::read_dir,
    path::{Path, PathBuf},
};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
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

fn do_run_dir(run_dir: &Path) -> Result<()> {
    let dump_final = DumpFile::read(&run_dir.join("dump.final"), &[])?;
    let snapshot_final = clusterize_snapshot(dump_final.get_snapshots()[0], 3.0);
    let counts = get_cluster_counts(&snapshot_final);
    let cluster_ids = counts
        .into_iter()
        .filter(|(_, count)| *count < 1000)
        .map(|(id, _)| id)
        .collect::<Vec<_>>();
    let clusters = snapshot_final.get_property("cluster");
    let mut sputter_indices = Vec::new();
    let mut no_sputter_indices = Vec::new();
    for (i, cluster) in clusters
        .iter()
        .map(|id| *id as usize)
        .enumerate()
        .take(snapshot_final.atoms_count)
    {
        if cluster_ids.contains(&cluster) {
            sputter_indices.push(i);
        } else {
            no_sputter_indices.push(i);
        }
    }
    let snapshot_sputter = copy_snapshot_with_indices(&snapshot_final, sputter_indices.into_iter());
    println!("sputtered: {}", snapshot_sputter.atoms_count);
    let dump_sputter = DumpFile::new(vec![snapshot_sputter]);
    dump_sputter.save(&run_dir.join("dump.sputter"))?;
    let snapshot_no_sputter =
        copy_snapshot_with_indices(&snapshot_final, no_sputter_indices.into_iter());
    let dump_no_sputter = DumpFile::new(vec![snapshot_no_sputter]);
    dump_no_sputter.save(&run_dir.join("dump.no_sputter"))?;
    Ok(())
}

fn do_results_dir(results_dir: &Path, threads: usize) -> Result<()> {
    let tp = ThreadPoolBuilder::new().num_threads(threads).build()?;
    let re = Regex::new(r"^run_(\d+)$")?;
    let mut entries = Vec::new();
    for entry in read_dir(results_dir)? {
        let run_dir = entry?.path();
        if re.is_match(&run_dir.file_name().unwrap_or_default().to_string_lossy()) {
            entries.push(run_dir);
        }
    }
    tp.install(|| {
        entries
            .into_par_iter()
            .map(|dir| do_run_dir(&dir))
            .collect::<Result<Vec<_>>>()
    })?;
    Ok(())
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match &cli.command {
        Commands::Single(args) => do_run_dir(&args.run_dir)?,
        Commands::Multi(args) => do_results_dir(&args.results_dir, args.threads)?,
    };
    Ok(())
}

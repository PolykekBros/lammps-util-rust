use anyhow::Result;
use clap::{Args, Parser, Subcommand};
use lammps_util_rust::{crater_snapshot, DumpFile, DumpSnapshot};
use log::debug;
use rayon::{prelude::*, ThreadPoolBuilder};
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

    /// Max depth (A)
    #[arg(short, long, default_value_t = 50.0)]
    max_depth: f64,

    /// Cutoff (A)
    #[arg(short, long, default_value_t = 1.75)]
    cutoff: f64,
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

fn get_crater_info(snapshot: &DumpSnapshot, zero_lvl: f64) -> String {
    let z = snapshot.get_property("z");
    let mut crater_count = 0;
    let mut surface_count = 0;
    let mut z_avg = 0.0;
    let mut z_min = f64::INFINITY;
    for &z in z {
        if z > -2.4 * 0.707 + zero_lvl {
            surface_count += 1;
        }
        crater_count += 1;
        z_min = z_min.min(z - zero_lvl);
        z_avg += z - zero_lvl;
    }
    z_avg /= crater_count as f64;
    let volume = crater_count as f64 * 20.1;
    let surface = surface_count as f64 * 7.3712;
    format!("{crater_count} {volume} {surface} {z_avg} {z_min}")
}

fn analyze_single_dir(dir: &Path, cutoff: f64, _depth: f64) -> Result<String> {
    let dump_input = DumpFile::read(&dir.join("dump.initial"), &[])?;
    let snapshot_input = dump_input.get_snapshots()[0];
    let zero_lvl = snapshot_input.get_zero_lvl();

    let dump_final = DumpFile::read(&dir.join("dump.final_no_cluster"), &[])?;
    let snapshot_final = dump_final.get_snapshots()[0];

    let snapshot_crater = crater_snapshot(snapshot_input, snapshot_final, cutoff, 3.0);
    debug!("crater atoms: {}", snapshot_crater.atoms_count);
    let info = get_crater_info(&snapshot_crater, zero_lvl);

    let dump_crater = DumpFile::new(vec![snapshot_crater]);
    dump_crater.save(&dir.join("dump.crater"))?;

    Ok(info)
}

fn analyze_results_dir(dir: &Path, threads: usize, cutoff: f64, _depth: f64) -> Result<String> {
    let tp = ThreadPoolBuilder::new().num_threads(threads).build()?;
    let re = Regex::new(r"^run_(\d+)$")?;
    let mut entries = Vec::new();
    for entry in read_dir(dir)? {
        let run_dir = entry?.path();
        if let Some(num) = re
            .captures(&run_dir.file_name().unwrap_or_default().to_string_lossy())
            .and_then(|caps| caps.get(1))
            .and_then(|m| m.as_str().parse::<usize>().ok())
        {
            entries.push((num, run_dir));
        };
    }
    let mut results = tp.install(|| {
        entries
            .into_par_iter()
            .map(|(num, dir)| analyze_single_dir(&dir, cutoff, _depth).map(|info| (num, info)))
            .collect::<Result<Vec<_>>>()
    })?;
    results.sort_by(|a, b| a.0.cmp(&b.0));
    let result = results
        .iter()
        .map(|(num, info)| format!("{num} {info}"))
        .collect::<Vec<_>>()
        .join("\n");
    Ok(result)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let info = match &cli.command {
        Commands::Single(args) => analyze_single_dir(&args.run_dir, cli.cutoff, cli.max_depth)?,
        Commands::Multi(args) => {
            analyze_results_dir(&args.results_dir, args.threads, cli.cutoff, cli.max_depth)?
        }
    };
    println!("{info}");
    Ok(())
}

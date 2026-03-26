use anyhow::Result;
use clap::{Args, Parser, Subcommand};
use lammps_util_rust::{crater_snapshot, process_results_dir, DumpFile, DumpSnapshot, IteratorAvg};
use log::debug;
use std::path::{Path, PathBuf};

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
    let surface_threshold = -2.4 * 0.707;
    let z = snapshot
        .get_property("z")
        .iter()
        .map(|z| z - zero_lvl)
        .collect::<Vec<_>>();
    let crater_count = z.len();
    let surface_count = z.iter().filter(|&&z| z > surface_threshold).count();
    let z_min = z.iter().copied().reduce(|acc, z| acc.min(z)).unwrap();
    let z_avg = z.iter().copied().avg().unwrap();
    let volume = crater_count as f64 * 20.1;
    let surface = surface_count as f64 * 7.3712;
    format!("{crater_count} {volume} {surface} {z_avg} {z_min}")
}

fn analyze_single_run(dir: &Path, cutoff: f64, _depth: f64) -> Result<String> {
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

fn analyze_results_dir(dir: &Path, threads: usize, cutoff: f64, depth: f64) -> Result<String> {
    let results = process_results_dir(dir, threads, |dir| {
        analyze_single_run(&dir.path, cutoff, depth)
    })?;
    let info = results
        .iter()
        .map(|(dir, info)| format!("{} {info}", dir.num))
        .collect::<Vec<_>>()
        .join("\n");
    Ok(info)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    let info = match &cli.command {
        Commands::Single(args) => analyze_single_run(&args.run_dir, cli.cutoff, cli.max_depth)?,
        Commands::Multi(args) => {
            analyze_results_dir(&args.results_dir, args.threads, cli.cutoff, cli.max_depth)?
        }
    };
    println!("{info}");
    Ok(())
}

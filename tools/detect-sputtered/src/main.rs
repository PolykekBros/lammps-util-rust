use anyhow::Result;
use clap::{Args, Parser, Subcommand};
use lammps_util::{
    DumpFile, clusterize_snapshot, copy_snapshot_with_indices, get_clusters, process_results_dir,
};
use std::path::{Path, PathBuf};

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
}

fn do_run_dir(run_dir: &Path) -> Result<usize> {
    let dump_final = DumpFile::read(&run_dir.join("dump.final"), &[])?;
    let snapshot_final = clusterize_snapshot(&dump_final.get_snapshots()[0], 3.0);
    let clusters = get_clusters(&snapshot_final);
    let sputter_indices = clusters
        .values()
        .filter(|atoms| atoms.len() < 1000)
        .flat_map(|atoms| atoms.iter().copied());
    let no_sputter_indices = clusters
        .values()
        .filter(|atoms| atoms.len() >= 1000)
        .flat_map(|atoms| atoms.iter().copied());
    let snapshot_sputter = copy_snapshot_with_indices(&snapshot_final, sputter_indices);
    let sputter_count = snapshot_sputter.get_atoms_count();
    DumpFile::new(vec![snapshot_sputter]).save(&run_dir.join("dump.sputter"))?;
    let snapshot_no_sputter =
        copy_snapshot_with_indices(&snapshot_final, no_sputter_indices.into_iter());
    DumpFile::new(vec![snapshot_no_sputter]).save(&run_dir.join("dump.no_sputter"))?;
    Ok(sputter_count)
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    match &cli.command {
        Commands::Single(args) => {
            do_run_dir(&args.run_dir)?;
        }
        Commands::Multi(args) => {
            process_results_dir(&args.results_dir, |run_dir| {
                let result = do_run_dir(&run_dir.path);
                println!("{} {result:?}", run_dir.path.to_string_lossy());
                result
            })?;
        }
    }
    Ok(())
}

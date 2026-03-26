use anyhow::Result;
use clap::{Parser, Subcommand};
use lammps_util_rust::{
    DumpFile, DumpSnapshot, clusterize_snapshot, copy_snapshot_with_indices, get_cluster_counts,
};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

/// A utility for processing LAMMPS simulation files.
#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Processes a LAMMPS data file based on a final dump file to remove sputtered atoms.
    Input(InputArgs),
    /// Processes a single LAMMPS dump file.
    Dump(DumpArgs),
}

/// Arguments for the 'input' subcommand.
#[derive(Parser)]
struct InputArgs {
    /// Lammps data file from which to remove sputtered atoms
    input_file: PathBuf,

    /// Final simulation dump file
    dump_final: PathBuf,

    /// Resulting data file
    output_file: PathBuf,
}

/// Arguments for the 'dump' subcommand.
#[derive(Parser)]
struct DumpArgs {
    /// Path to the LAMMPS dump file to be processed.
    dump_final: PathBuf,
}

fn get_indices_to_delete(snapshot: &DumpSnapshot) -> Vec<usize> {
    let snapshot = clusterize_snapshot(snapshot, 3.0);
    let clusters_to_delete = get_cluster_counts(&snapshot)
        .into_iter()
        .filter(|(_, count)| *count < 1000)
        .map(|(id, _)| id)
        .collect::<Vec<_>>();
    snapshot
        .get_property("cluster")
        .iter()
        .enumerate()
        .filter(|(_, c)| clusters_to_delete.contains(&(**c as usize)))
        .map(|(i, _)| i)
        .collect()
}

fn get_ids_to_delete(snapshot: &DumpSnapshot) -> Vec<usize> {
    let indices_to_delete = get_indices_to_delete(snapshot);
    snapshot
        .get_property("id")
        .iter()
        .enumerate()
        .filter(|(i, _)| indices_to_delete.contains(i))
        .map(|(_, id)| *id as usize)
        .collect()
}

fn delete_atoms(in_file: &Path, out_file: &Path, ids: &[usize]) -> Result<()> {
    let reader = BufReader::new(File::open(in_file)?);
    let mut writer = BufWriter::new(File::create(out_file)?);
    for (cnt, line) in reader.lines().enumerate() {
        let cnt = cnt + 1;
        let line = line?;
        let mut tokens = line.split(' ');
        let first_token = tokens.next().unwrap();
        if cnt == 3 {
            let atom_cnt = first_token.parse::<usize>()?;
            let new_atom_cnt = atom_cnt - ids.len();
            writeln!(writer, "{} atoms", new_atom_cnt)?;
        } else if cnt < 17 || tokens.next().is_none() {
            writeln!(writer, "{}", line)?;
        } else {
            let atom_id = first_token.parse::<usize>()?;
            if !ids.contains(&atom_id) {
                writeln!(writer, "{}", line)?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}

fn process_input(args: &InputArgs) -> Result<()> {
    let dump_final = DumpFile::read(&args.dump_final, &[])?;
    let ids_to_delete = get_ids_to_delete(dump_final.get_snapshots()[0]);
    println!("about to delete {} atoms", ids_to_delete.len());
    delete_atoms(&args.input_file, &args.output_file, &ids_to_delete)?;
    println!("deleted {} atoms", ids_to_delete.len());
    Ok(())
}

fn process_dump(args: &DumpArgs) -> Result<()> {
    let dump_final = DumpFile::read(&args.dump_final, &[])?;
    let snapshot = dump_final.get_snapshots()[0];
    let indices_to_delete = get_indices_to_delete(snapshot);
    println!("about to delete {} atoms", indices_to_delete.len());
    let indices_to_keep = (0..snapshot.atoms_count).filter(|i| !indices_to_delete.contains(i));
    let snapshot = copy_snapshot_with_indices(snapshot, indices_to_keep);
    let dump_no_sputter = DumpFile::new(vec![snapshot]);
    let path = args.dump_final.with_file_name(format!(
        "{}_no_sputter",
        args.dump_final
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
    ));
    dump_no_sputter.save(&path)?;
    Ok(())
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match &cli.command {
        Commands::Input(args) => process_input(args)?,
        Commands::Dump(args) => process_dump(args)?,
    }

    Ok(())
}

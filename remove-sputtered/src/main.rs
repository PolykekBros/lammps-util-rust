use anyhow::Result;
use clap::Parser;
use lammps_util_rust::{DumpFile, clusterize_snapshot, get_cluster_counts};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Lammps data file from which to remove sputtered atoms
    input_file: PathBuf,

    /// Final simulation dump file
    dump_final: PathBuf,

    /// Resulting data file
    output_file: PathBuf,
}

fn get_ids_to_delete(dump_path: &Path) -> Result<Vec<usize>> {
    let dump = DumpFile::read(dump_path, &[])?;
    let snapshot = clusterize_snapshot(dump.get_snapshots()[0], 3.0);
    let counts = get_cluster_counts(&snapshot);
    let ids_to_delete = counts
        .into_iter()
        .filter(|(_, count)| *count < 1000)
        .map(|(id, _)| id)
        .collect::<Vec<_>>();
    Ok(ids_to_delete)
}

fn delete_atoms(in_file: &Path, out_file: &Path, ids: &[usize]) -> Result<()> {
    let reader = BufReader::new(File::open(in_file)?);
    let mut writer = BufWriter::new(File::open(out_file)?);
    for (cnt, line) in reader.lines().enumerate() {
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

fn main() -> Result<()> {
    let cli = Cli::parse();
    let ids_to_delete = get_ids_to_delete(&cli.dump_final)?;
    println!("about to delete {} atoms", ids_to_delete.len());
    delete_atoms(&cli.input_file, &cli.output_file, &ids_to_delete)?;
    println!("deleted {} atoms", ids_to_delete.len());
    Ok(())
}

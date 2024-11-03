use clap::Parser;
use lammps_util_rust::{DumpFile, DumpSnapshot};
use plotters::prelude::*;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[arg(value_name = "DUMP_FILE")]
    dump_file: PathBuf,

    #[arg(short, long, value_name = "TIMESTEP")]
    timestep: Option<u64>,

    #[arg(short, long, value_name = "DELTA", default_value_t = 5.43 * 2.0)]
    delta: f64,

    #[arg(short, long, value_name = "OUTPUT_DIR")]
    output_dir: PathBuf,
}

fn get_distribution(dump: &DumpSnapshot, delta: f64) -> (Vec<f64>, Vec<Vec<f64>>) {
    println!("mii");
    let atom_type = dump.get_property("type");
    let atom_x = dump.get_property("x");
    let atom_y = dump.get_property("y");
    let atom_z = dump.get_property("z");
    println!("mii");
    let types = atom_type
        .iter()
        .copied()
        .map(|t| t as usize)
        .collect::<HashSet<_>>();
    println!("nipah");
    let x_min = atom_x.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = atom_x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let y_min = atom_y.iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = atom_y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let z_min = atom_z.iter().copied().fold(f64::INFINITY, f64::min);
    let z_max = atom_z.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    println!("{x_min}, {x_max} | {y_min}, {y_max} | {z_min}, {z_max}");
    let volume = delta * (y_max - y_min) * (x_max - x_min);
    let count = (z_max - z_min).ceil() as usize;
    let plot_x = (0..count)
        .map(|i| (i as f64) * delta + delta / 2.0f64)
        .collect();
    let plot_y = types
        .iter()
        .map(|id| (id, std::iter::zip(atom_type, atom_z)))
        .map(|(&id, a_t_z)| {
            (0..count)
                .map(|i| i as f64)
                .map(|i| (i * delta, (i + 1.0f64) * delta))
                .map(|(start, end)| {
                    a_t_z
                        .clone()
                        .filter(|(&a_id, &a_z)| {
                            (a_z >= start && a_z < end) && (a_id as usize) == id
                        })
                        .count() as f64
                        / volume
                })
                .collect()
        })
        .collect();
    (plot_x, plot_y)
}

fn plot_distribution(
    dir: &Path,
    plot_x: &Vec<f64>,
    plot_y: &Vec<Vec<f64>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let plot_path = dir.join("distribution.png");
    let root = BitMapBackend::new(&plot_path, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let y_max = plot_y.iter().flatten().fold(0.0f64, |acc, &y| acc.max(y));
    let mut chart = ChartBuilder::on(&root)
        .caption("Density Distribution", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(
            *plot_x.first().unwrap()..*plot_x.last().unwrap(),
            0f64..y_max,
        )?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(
        std::iter::zip(plot_x.to_owned(), plot_y[0].to_owned()),
        &RED,
    ))?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    let dump_path = cli.dump_file;
    let timesteps = match cli.timestep {
        Some(timestep) => vec![timestep],
        _ => Vec::new(),
    };
    println!("start");
    let dump = DumpFile::new(dump_path.as_path(), &timesteps)?;
    println!("read: {:?}", dump.get_snapshots()[0]);
    let (plot_x, plot_y) = get_distribution(dump.get_snapshots()[0], cli.delta);
    println!("calculated distribution dump");
    plot_distribution(&cli.output_dir, &plot_x, &plot_y)?;
    Ok(())
}

use clap::Parser;
use geomutil::{
    geomutil_triangulation::alpha_shape_2d,
    geomutil_util::{points_bounding_box, Point2D},
};
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

fn plot_slices(dump: &DumpSnapshot, delta: f64) {
    let coords = dump.get_coordinates();
    let z_min = 40.0;
    let z_max = coords
        .iter()
        .map(|c| c.z)
        .max_by(|a, b| a.total_cmp(b))
        .unwrap();
    let n = ((z_max - z_min) / delta).ceil() as usize;
    let mut slices = vec![Vec::new(); n];
    for c in coords {
        if c.z < z_min {
            continue;
        }
        let i = ((c.z - z_min) / delta).floor() as usize;
        slices[i].push(Point2D::new(c.x as f32, c.y as f32));
    }
    for (i, points) in slices.into_iter().enumerate() {
        let file_name = format!("triangles_{}.png", i);
        let root = BitMapBackend::new(&file_name, (800, 600)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let (low_boundary, up_boundary) = points_bounding_box(&points).unwrap();
        let mut chart = ChartBuilder::on(&root)
            .margin(10)
            .build_cartesian_2d(
                low_boundary.x - 1.0..up_boundary.x + 1.0,
                low_boundary.y - 1.0..up_boundary.y + 1.0,
            )
            .unwrap();
        chart.configure_mesh().draw().unwrap();

        let shapes = alpha_shape_2d(&points, 1.0).unwrap();
        for t in shapes.into_iter().map(|s| s.triangles).flatten() {
            let line_series = LineSeries::new(
                vec![t.a.xy(), t.b.xy(), t.c.xy(), t.a.xy()],
                BLACK.stroke_width(1),
            );
            chart.draw_series(line_series).unwrap();
        }
    }
}

fn get_distribution(dump: &DumpSnapshot, delta: f64) -> (Vec<f64>, Vec<Vec<f64>>) {
    let atom_type = dump.get_property("type");
    let atom_x = dump.get_property("x");
    let atom_y = dump.get_property("y");
    let atom_z = dump.get_property("z");
    let types = atom_type
        .iter()
        .copied()
        .map(|t| t as usize)
        .collect::<HashSet<_>>();
    let x_min = atom_x.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = atom_x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let y_min = atom_y.iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = atom_y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let z_min = atom_z.iter().copied().fold(f64::INFINITY, f64::min);
    let z_max = atom_z.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    println!("{x_min}, {x_max} | {y_min}, {y_max} | {z_min}, {z_max}");
    let volume = delta * (y_max - y_min) * (x_max - x_min);
    let count = ((z_max - z_min) / delta).ceil() as usize;
    let plot_x = (0..count)
        .map(|i| (i as f64) * delta + delta / 2.0f64)
        .collect();
    let plot_y = types
        .into_iter()
        .map(|t| {
            (0..count)
                .map(|i| i as f64)
                .map(|i| (i * delta, (i + 1.0f64) * delta))
                .map(|(start, end)| {
                    atom_type
                        .iter()
                        .zip(atom_z)
                        .filter(|(&a_t, &a_z)| (a_z >= start && a_z < end) && (a_t as usize) == t)
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
        .caption("Density Distribution", ("sans-serif", 30).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(
            *plot_x.first().unwrap()..*plot_x.last().unwrap(),
            0f64..y_max,
        )?;
    chart
        .configure_mesh()
        .x_max_light_lines(0)
        .y_max_light_lines(0)
        .draw()?;
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
    let dump = DumpFile::read(dump_path.as_path(), &timesteps)?;
    let snapshot = dump.get_snapshots()[0];
    plot_slices(&snapshot, cli.delta);
    // let (plot_x, plot_y) = get_distribution(dump.get_snapshots()[0], cli.delta);
    // println!("calculated distribution dump");
    // plot_distribution(&cli.output_dir, &plot_x, &plot_y)?;
    Ok(())
}

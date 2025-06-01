use clap::Parser;
use geomutil::{
    geomutil_triangulation::alpha_shape_2d,
    geomutil_util::{points_bounding_box, Point2D, Shape2D},
};
use lammps_util_rust::{DumpFile, DumpSnapshot};
use log::info;
use plotters::prelude::*;
use std::collections::{HashMap, HashSet};
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

    /// Save slice pictures
    #[arg(short, long)]
    pictures: bool,

    #[arg(short, long, value_name = "OUTPUT_DIR")]
    output_dir: PathBuf,
}

#[derive(Default, Clone)]
struct Particle {
    xy: Point2D,
    _id: usize,
    ptype: usize,
}

impl Particle {
    fn new(xy: Point2D, id: usize, ptype: usize) -> Self {
        Self { xy, _id: id, ptype }
    }
}

#[derive(Default, Clone)]
struct Slice {
    particles: Vec<Particle>,
    z_min: f64,
    z_max: f64,
    shapes: Vec<Shape2D>,
}

impl Slice {
    fn bounding_box(&self) -> (Point2D, Point2D) {
        let points = self.points();
        points_bounding_box(&points).unwrap()
    }

    fn points(&self) -> Vec<Point2D> {
        self.particles.iter().map(|p| p.xy).collect()
    }
}

fn get_slices(dump: &DumpSnapshot, delta: f64) -> Vec<Slice> {
    let coords = dump.get_coordinates();
    let types = dump.get_property("type");
    let z_min = coords
        .iter()
        .map(|c| c.z)
        .min_by(|a, b| a.total_cmp(b))
        .unwrap();
    let z_max = coords
        .iter()
        .map(|c| c.z)
        .max_by(|a, b| a.total_cmp(b))
        .unwrap();
    let n = ((z_max - z_min) / delta).ceil() as usize;
    let mut slices = vec![Slice::default(); n];
    for c in coords {
        if c.z < z_min {
            continue;
        }
        let i = ((c.z - z_min) / delta).floor() as usize;
        if slices[i].particles.is_empty() {
            slices[i].z_min = z_min + (i as f64) * delta;
            slices[i].z_max = slices[i].z_min + delta;
        }
        slices[i].particles.push(Particle::new(
            Point2D::new(c.x as f32, c.y as f32),
            c.index(),
            types[c.index()] as usize,
        ));
    }
    for slice in slices.iter_mut() {
        let points = &slice.points();
        info!("before shape, points.len: {}", points.len());
        slice.shapes = alpha_shape_2d(&slice.points(), 0.5).unwrap();
        info!("after shape");
    }
    slices
}

fn plot_slices(out_dir: &Path, slices: &[Slice]) {
    for (i, slice) in slices.iter().enumerate() {
        let path = out_dir.join(format!("slice_{}.png", i));
        let root = BitMapBackend::new(&path, (800, 800)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let (low_boundary, up_boundary) = slice.bounding_box();
        let mut chart = ChartBuilder::on(&root)
            .margin(10)
            .build_cartesian_2d(
                low_boundary.x - 1.0..up_boundary.x + 1.0,
                low_boundary.y - 1.0..up_boundary.y + 1.0,
            )
            .unwrap();
        chart.configure_mesh().draw().unwrap();
        for t in slice.shapes.iter().flat_map(|s| &s.triangles) {
            let line_series = LineSeries::new(
                vec![t.a.xy(), t.b.xy(), t.c.xy(), t.a.xy()],
                BLACK.stroke_width(1),
            );
            chart.draw_series(line_series).unwrap();
        }
    }
}

fn distribution_data(slices: &[Slice], delta: f64) -> (Vec<f64>, Vec<Vec<f64>>) {
    let x = slices.iter().map(|s| (s.z_max + s.z_min) / 2.0).collect();
    let mut types = HashSet::new();
    for s in slices.iter() {
        for p in s.particles.iter() {
            types.insert(p.ptype);
        }
    }
    let y = slices
        .iter()
        .map(|s| {
            let volume = s.shapes.iter().map(|s| s.area()).sum::<f32>() as f64 * delta;
            let mut counts = HashMap::new();
            info!("volume: {}", volume);
            for t in types.iter() {
                counts.insert(*t, 0.0);
            }
            for p in s.particles.iter() {
                counts.entry(p.ptype).and_modify(|c| *c += 1.0);
            }
            info!("counts: {counts:?}");
            for t in types.iter() {
                counts.entry(*t).and_modify(|c| *c /= volume);
            }
            counts
        })
        .collect::<Vec<HashMap<usize, f64>>>();
    let mut y_transposed = Vec::with_capacity(types.len());
    for t in types.iter() {
        y_transposed.push(y.iter().map(|v| v[t]).collect());
    }
    (x, y_transposed)
}

fn plot_distribution(
    dir: &Path,
    plot_x: &Vec<f64>,
    plot_y: &[Vec<f64>],
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
        .y_desc("n / volume")
        .x_desc("z coordinate (A)")
        .x_max_light_lines(0)
        .y_max_light_lines(0)
        .draw()?;
    chart
        .draw_series(LineSeries::new(
            std::iter::zip(plot_x.to_owned(), plot_y[0].to_owned()),
            &RED,
        ))?
        .label("Si");
    chart
        .draw_series(LineSeries::new(
            std::iter::zip(plot_x.to_owned(), plot_y[1].to_owned()),
            &BLUE,
        ))?
        .label("C");
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let cli = Cli::parse();
    let dump_path = cli.dump_file;
    let timesteps = match cli.timestep {
        Some(timestep) => vec![timestep],
        _ => vec![],
    };
    let dump = DumpFile::read(dump_path.as_path(), &timesteps)?;
    let snapshot = dump.get_snapshots()[0];
    let slices = get_slices(snapshot, cli.delta);
    if cli.pictures {
        plot_slices(&cli.output_dir, &slices);
    }
    let (plot_x, plot_y) = distribution_data(&slices, cli.delta);
    plot_distribution(&cli.output_dir, &plot_x, &plot_y)?;
    for (i, x) in plot_x.iter().enumerate() {
        print!("{}", x);
        print!("\t{}", plot_y[0][i]);
        println!("\t{}", plot_y[1][i]);
    }
    Ok(())
}

use anyhow::Result;
use clap::{Args, Parser, Subcommand};
use colorgrad::preset::viridis;
use geomutil_util::{Point2, Point3};
use lammps_util_rust::{process_results_dir, DumpFile, IteratorAvg, XYZ};
use plotters::{
    chart::ChartBuilder,
    prelude::{BitMapBackend, IntoDrawingArea},
    style::WHITE,
};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

mod heatmap;
use heatmap::{heatmap, Colorbar, Domain};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Width of the approximation square (A)
    #[arg(short, long, default_value_t = 5.43 / 2.0)]
    width: f64,

    /// Zero level of the crystal surface (A)
    #[arg(short, long)]
    zero_lvl: f64,
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

struct SurfaceValues {
    data: Vec<f32>,
    domain: Domain,
    square_width: f32,
    x_count: usize,
    y_count: usize,
}

impl SurfaceValues {
    pub fn new(domain: Domain, square_width: f32) -> Self {
        let x_count = (domain.width() / square_width).ceil() as usize;
        let y_count = (domain.height() / square_width).ceil() as usize;
        Self {
            data: vec![f32::NAN; x_count * y_count],
            square_width,
            domain,
            x_count,
            y_count,
        }
    }

    fn populate(&mut self, xyz: impl IntoIterator<Item = XYZ>) {
        let x_lo = self.domain.lo().x;
        let y_lo = self.domain.lo().y;
        let x_cnt_f32 = self.x_count as f32;
        let y_cnt_f32 = self.y_count as f32;
        let square_width = self.square_width;
        xyz.into_iter()
            .map(|p| {
                Point3::from([
                    (p.x as f32 - x_lo) / square_width,
                    (p.y as f32 - y_lo) / square_width,
                    p.z as f32,
                ])
            })
            .filter(|p| p.x >= 0.0 && p.x <= x_cnt_f32 && p.y >= 0.0 && p.y <= y_cnt_f32)
            .for_each(|p| {
                let x_i = (p.x as usize).min(self.x_count - 1);
                let y_i = (p.y as usize).min(self.y_count - 1);
                *self.at_mut(x_i, y_i) = self.at(x_i, y_i).max(p.z);
            });
    }

    fn interpolate(&mut self) {
        (0..self.x_count)
            .flat_map(|x_i| (0..self.y_count).map(move |y_i| (x_i, y_i)))
            .filter(|(x_i, y_i)| self.at(*x_i, *y_i).is_nan())
            .collect::<Vec<_>>()
            .into_iter()
            .for_each(|(x_i, y_i)| {
                let (sum, cnt) = (-1..=1)
                    .flat_map(|x_i_offset| (-1..=1).map(move |y_i_offset| (x_i_offset, y_i_offset)))
                    .map(|(x_i_offset, y_i_offset)| {
                        *self.at_i64(x_i as i64 + x_i_offset, y_i as i64 + y_i_offset)
                    })
                    .filter(|value| !value.is_nan())
                    .fold((0.0, 0), |(sum, cnt), value| (sum + value, cnt + 1));
                if cnt != 0 {
                    *self.at_mut(x_i, y_i) = sum / cnt as f32;
                }
            });
    }

    fn get_i(&self, x_i: usize, y_i: usize) -> usize {
        x_i * self.y_count + y_i
    }

    pub fn at(&self, x_i: usize, y_i: usize) -> &f32 {
        assert!(x_i < self.x_count);
        assert!(y_i < self.y_count);
        let i = self.get_i(x_i, y_i);
        &self.data[i]
    }

    pub fn at_mut(&mut self, x_i: usize, y_i: usize) -> &mut f32 {
        assert!(x_i < self.x_count);
        assert!(y_i < self.y_count);
        let i = self.get_i(x_i, y_i);
        &mut self.data[i]
    }

    pub fn at_i64(&self, x_i: i64, y_i: i64) -> &f32 {
        if x_i < 0 || y_i < 0 || x_i >= self.x_count as i64 || y_i >= self.y_count as i64 {
            &f32::NAN
        } else {
            self.at(x_i as usize, y_i as usize)
        }
    }
}

fn get_surface_values(
    xyz: impl IntoIterator<Item = XYZ>,
    domain: &Domain,
    square_width: f32,
    zero_lvl: f32,
) -> SurfaceValues {
    assert!(square_width > 0.0);
    let mut values = SurfaceValues::new(domain.clone(), square_width);
    values.populate(xyz);
    log::debug!(
        "intial NaNs: {}",
        values.data.iter().filter(|v| v.is_nan()).count()
    );
    values.interpolate();
    log::debug!(
        "after interpolation NaNs: {}",
        values.data.iter().filter(|v| v.is_nan()).count()
    );
    values.data.iter_mut().for_each(|z| *z -= zero_lvl);
    values
}

fn plot_surface_2d(output_dir: &Path, values: &SurfaceValues) -> Result<()> {
    let plot_width: u32 = 640;
    let plot_height: u32 = 512;
    let plot_color_width: u32 = 80;

    let plot_path = output_dir.join("surface_2d.png");
    let drawing_area =
        BitMapBackend::new(&plot_path, (plot_width, plot_height)).into_drawing_area();
    drawing_area.fill(&WHITE)?;

    let pixel_width = drawing_area.dim_in_pixel().0;
    let (left, right) = drawing_area.split_horizontally(pixel_width - plot_color_width);
    let left = left.shrink(
        ((plot_width - plot_color_width - plot_height) / 2, 0),
        (plot_height, plot_height),
    );
    let colorbar = Colorbar::new(-20.0, 10.0, viridis());
    colorbar.draw(ChartBuilder::on(&right));

    let mut chart_builder = ChartBuilder::on(&left);
    chart_builder.margin(10);
    heatmap(values, &colorbar, chart_builder)?;
    Ok(())
}

fn write_surface_coords(path: &Path, values: &SurfaceValues) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "{} {}", values.domain.lo().x, values.domain.lo().y)?;
    writeln!(writer, "{} {}", values.domain.hi().x, values.domain.hi().y)?;
    for x_i in 0..values.x_count {
        write!(writer, "{}", values.at(x_i, 0))?;
        for y_i in 1..values.y_count {
            write!(writer, "\t{}", values.at(x_i, y_i))?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

fn save_results(path: &Path, values: &SurfaceValues) -> Result<()> {
    plot_surface_2d(path, values)?;
    write_surface_coords(&path.join("surface_coords.txt"), values)?;
    Ok(())
}

fn analyze_single_run(path: &Path, square_width: f64, zero_lvl: f64) -> Result<SurfaceValues> {
    let dump_final = DumpFile::read(&path.join("dump.final_no_cluster"), &[])?;
    let snapshot = dump_final.get_snapshots()[0];
    let domain = Domain::new(
        Point2::from([snapshot.sym_box.xlo as f32, snapshot.sym_box.ylo as f32]),
        Point2::from([snapshot.sym_box.xhi as f32, snapshot.sym_box.yhi as f32]),
    );
    let coords = snapshot.get_coordinates();
    let values = get_surface_values(coords, &domain, square_width as f32, zero_lvl as f32);
    save_results(path, &values)?;
    Ok(values)
}

fn avg_surface_values(values: &[SurfaceValues]) -> Option<(SurfaceValues, SurfaceValues)> {
    if values.is_empty() {
        return None;
    }
    let data = (0..values[0].data.len())
        .map(|i| {
            values
                .iter()
                .map(|values| values.data[i])
                .avg_with_std()
                .unwrap()
        })
        .collect::<Vec<_>>();
    Some((
        SurfaceValues {
            data: data.iter().map(|(avg, _)| avg).copied().collect(),
            domain: values[0].domain.clone(),
            ..values[0]
        },
        SurfaceValues {
            data: data.iter().map(|(_, std)| std).copied().collect(),
            domain: values[0].domain.clone(),
            ..values[0]
        },
    ))
}

fn analyze_results_dir(dir: &Path, threads: usize, square_width: f64, zero_lvl: f64) -> Result<()> {
    let values = process_results_dir(dir, threads, move |dir| {
        analyze_single_run(&dir.path, square_width, zero_lvl)
    })?;
    let (avg, std) = avg_surface_values(
        &values
            .into_iter()
            .map(|(_, values)| values)
            .collect::<Vec<_>>(),
    )
    .unwrap();
    save_results(dir, &avg)?;
    write_surface_coords(&dir.join("surface_coords_errors.txt"), &std)?;
    Ok(())
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    match &cli.command {
        Commands::Single(args) => {
            analyze_single_run(&args.run_dir, cli.width, cli.zero_lvl)?;
        }
        Commands::Multi(args) => {
            analyze_results_dir(&args.results_dir, args.threads, cli.width, cli.zero_lvl)?
        }
    };
    Ok(())
}

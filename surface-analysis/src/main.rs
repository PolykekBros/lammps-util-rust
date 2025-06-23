use anyhow::{anyhow, Result};
use clap::{Args, Parser, Subcommand};
use colorgrad::preset::viridis;
use geomutil_util::{Point2, Point3};
use itertools::Itertools;
use lammps_util_rust::{get_runs_dirs, DumpFile, RunDir, XYZ};
use plotters::{
    chart::ChartBuilder,
    prelude::{BitMapBackend, IntoDrawingArea},
    style::WHITE,
};
use rayon::{prelude::*, ThreadPoolBuilder};
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
    x_count: usize,
    y_count: usize,
}

impl SurfaceValues {
    pub fn new(x_count: usize, y_count: usize) -> Self {
        Self {
            data: vec![f32::NAN; (x_count * y_count) as usize],
            x_count,
            y_count,
        }
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
    let x_count = (domain.width() / square_width as f32).ceil() as usize;
    let y_count = (domain.height() / square_width as f32).ceil() as usize;
    let mut values = SurfaceValues::new(x_count, y_count);
    xyz.into_iter()
        .map(|p| {
            Point3::from([
                (p.x as f32 - domain.lo().x) / square_width,
                (p.y as f32 - domain.lo().y) / square_width,
                p.z as f32 - zero_lvl,
            ])
        })
        .filter(|p| p.x >= 0.0 && p.x <= x_count as f32 && p.y >= 0.0 && p.y <= y_count as f32)
        .for_each(|p| {
            let x_i = (p.x as usize).min(x_count - 1);
            let y_i = (p.y as usize).min(y_count - 1);
            *values.at_mut(x_i, y_i) = values.at(x_i, y_i).max(p.z);
        });
    log::debug!(
        "intial NaNs: {}",
        values.data.iter().filter(|v| v.is_nan()).count()
    );
    (0..x_count)
        .flat_map(|x_i| (0..y_count).map(move |y_i| (x_i, y_i)))
        .filter(|(x_i, y_i)| values.at(*x_i, *y_i).is_nan())
        .collect::<Vec<_>>()
        .into_iter()
        .for_each(|(x_i, y_i)| {
            let (sum, cnt) = (-1..=1)
                .flat_map(|x_i_offset| (-1..=1).map(move |y_i_offset| (x_i_offset, y_i_offset)))
                .map(|(x_i_offset, y_i_offset)| {
                    *values.at_i64(x_i as i64 + x_i_offset, y_i as i64 + y_i_offset)
                })
                .filter(|value| !value.is_nan())
                .fold((0.0, 0), |(sum, cnt), value| (sum + value, cnt + 1));
            if cnt != 0 {
                *values.at_mut(x_i, y_i) = sum / cnt as f32;
            }
        });
    log::debug!(
        "after interpolation NaNs: {}",
        values.data.iter().filter(|v| v.is_nan()).count()
    );
    values
}

fn plot_surface_2d(output_dir: &Path, values: &SurfaceValues, domain: &Domain) -> Result<()> {
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
    heatmap(values, domain, &colorbar, chart_builder)?;
    Ok(())
}

fn write_surface_coords(path: &Path, domain: &Domain, values: &SurfaceValues) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "{} {}", domain.lo().x, domain.lo().y)?;
    writeln!(writer, "{} {}", domain.hi().x, domain.hi().y)?;
    for x_i in 0..values.x_count {
        write!(writer, "{}", values.at(x_i, 0))?;
        for y_i in 1..values.y_count {
            write!(writer, "\t{}", values.at(x_i, y_i))?;
        }
        write!(writer, "\n")?;
    }
    Ok(())
}

fn analyze_single_run(path: &Path, square_width: f64, zero_lvl: f64) -> Result<SurfaceValues> {
    let dump_final = DumpFile::read(&path.join("dump.final"), &[])?;
    let snapshot = dump_final.get_snapshots()[0];
    let domain = Domain::new(
        Point2::from([snapshot.sym_box.xlo as f32, snapshot.sym_box.ylo as f32]),
        Point2::from([snapshot.sym_box.xhi as f32, snapshot.sym_box.yhi as f32]),
    );
    let surface_coords_txt_path = path.join("surface_coords.txt");
    let coords = snapshot.get_coordinates();
    let values = get_surface_values(coords, &domain, square_width as f32, zero_lvl as f32);
    plot_surface_2d(path, &values, &domain)?;
    write_surface_coords(&surface_coords_txt_path, &domain, &values)?;
    Ok(values)
}

fn do_run_dir(
    run_dir: &RunDir,
    square_width: f64,
    zero_lvl: f64,
) -> Result<(usize, SurfaceValues)> {
    let values = analyze_single_run(&run_dir.path, square_width, zero_lvl)?;
    Ok((run_dir.num, values))
}

fn analyze_results_dir(
    dir: &Path,
    threads: usize,
    square_width: f64,
    zero_lvl: f64,
) -> Result<SurfaceValues> {
    let tp = ThreadPoolBuilder::new().num_threads(threads).build()?;
    let run_dirs = get_runs_dirs(dir)?;
    let mut results = tp.install(|| {
        run_dirs
            .into_par_iter()
            .map(|run_dir| do_run_dir(&run_dir, square_width, zero_lvl))
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

// fn process_results()

// fn main() -> Result<()> {
//     env_logger::init();
//     let cli = Cli::parse();

//     let dump_final_path = &cli.dump_final[0];
//     let values = cli.dump_final[1..]
//         .iter()
//         .try_fold(values, |values, dump_final_path| {
//             let values_new = parse_dump_final(dump_final_path, &domain, cli.width, cli.zero_lvl)?;
//             anyhow::Ok::<DMatrix<f64>>(values + values_new)
//         })?;
//     if let Some(output_dir) = cli.output_dir {
//         let values = values.scale(1.0 / cli.dump_final.len() as f64);
//         plot_surface_2d(&output_dir, &values, &domain)?;
//         let surface_coords_txt = output_dir.join("surface_coords.txt");
//         write_surface_coords(&surface_coords_txt, &domain, &values)?;
//     }
//     Ok(())
// }

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    let info = match &cli.command {
        Commands::Single(args) => analyze_single_run(&args.run_dir, cli.width, cli.zero_lvl)?,
        Commands::Multi(args) => {
            analyze_results_dir(&args.results_dir, args.threads, cli.width, cli.zero_lvl)?
        }
    };
    Ok(())
}

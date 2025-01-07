use anyhow::{anyhow, Result};
use clap::Parser;
use colorgrad::preset::viridis;
use itertools::izip;
use lammps_util_rust::{DumpFile, DumpSnapshot};
use nalgebra::{point, DMatrix, DVector, Point3};
use plotters::{
    chart::ChartBuilder,
    prelude::{BitMapBackend, IntoDrawingArea},
    style::WHITE,
};
use std::path::{Path, PathBuf};

mod heatmap;
use heatmap::{heatmap, Colorbar, Domain};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Input files of the final crystal
    #[arg()]
    dump_final: Vec<PathBuf>,

    /// Directory where to save resulting plots and tables
    #[arg(short, long)]
    output_dir: Option<PathBuf>,

    /// Width of the approximation square (A)
    #[arg(short, long, default_value_t = 5.43 / 2.0)]
    width: f64,

    /// Zero level of the crystal surface (A)
    #[arg(short, long)]
    zero_lvl: f64,
}

fn get_surface_values(xyz: &[Point3<f64>], domain: &Domain, square_width: f64) -> DMatrix<f64> {
    assert!(square_width > 0.0);
    let x_count = (domain.width() / square_width).ceil() as usize;
    let y_count = (domain.height() / square_width).ceil() as usize;
    let mut values = DMatrix::repeat(x_count, y_count, f64::NAN);
    xyz.iter()
        .copied()
        .filter(|p| {
            p.x >= domain.lo().x
                && p.x <= domain.hi().x
                && p.y >= domain.lo().y
                && p.y <= domain.hi().y
        })
        .for_each(|p| {
            let x_i = ((p.x - domain.lo().x) / square_width).floor() as usize;
            let y_i = ((p.y - domain.lo().y) / square_width).floor() as usize;
            let x_i = x_i.min(x_count - 1);
            let y_i = y_i.min(y_count - 1);
            values[(x_i, y_i)] = values[(x_i, y_i)].max(p.z);
        });
    let check_value = |vals: &DMatrix<f64>, x_i: i64, y_i: i64| {
        if x_i < 0 || y_i < 0 || x_i >= x_count as i64 || y_i >= y_count as i64 {
            f64::NAN
        } else {
            vals[(x_i as usize, y_i as usize)]
        }
    };
    println!("NaNs: {}", values.iter().filter(|v| v.is_nan()).count());
    (0..x_count as i64).for_each(|x_i| {
        (0..y_count as i64).for_each(|y_i| {
            if values[(x_i as usize, y_i as usize)].is_nan() {
                let (value_sum, value_cnt) = [
                    check_value(&values, x_i - 1, y_i - 1),
                    check_value(&values, x_i - 1, y_i),
                    check_value(&values, x_i - 1, y_i + 1),
                    check_value(&values, x_i, y_i - 1),
                    check_value(&values, x_i, y_i),
                    check_value(&values, x_i, y_i + 1),
                    check_value(&values, x_i + 1, y_i - 1),
                    check_value(&values, x_i + 1, y_i),
                    check_value(&values, x_i + 1, y_i + 1),
                ]
                .into_iter()
                .filter(|value| !value.is_nan())
                .fold((0.0, 0), |(sum, cnt), value| (sum + value, cnt + 1));
                values[(x_i as usize, y_i as usize)] = value_sum / value_cnt as f64;
            }
        });
    });
    values
}

fn plot_surface_2d(output_dir: &Path, values: &DMatrix<f64>, domain: &Domain) -> Result<()> {
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

fn parse_dump_final(
    snapshot: &DumpSnapshot,
    domain: &Domain,
    square_width: f64,
    zero_lvl: f64,
) -> DMatrix<f64> {
    let snapshot_x = DVector::from_column_slice(snapshot.get_property("x"));
    let snapshot_y = DVector::from_column_slice(snapshot.get_property("y"));
    let snapshot_z = DVector::from_column_slice(snapshot.get_property("z"));
    let xyz = izip![snapshot_x.iter(), snapshot_y.iter(), snapshot_z.iter()]
        .map(|(&x, &y, &z)| point![x, y, z])
        .collect::<Vec<_>>();
    let values = get_surface_values(&xyz, domain, square_width);
    println!("values.shape: {:?}", values.shape());
    values.add_scalar(-zero_lvl)
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let dump_final_path = &cli.dump_final[0];
    let dump_final = DumpFile::read(dump_final_path, &[])?;
    let snapshot = dump_final.get_snapshots()[0];
    let domain = Domain::new(
        point![snapshot.sym_box.xlo, snapshot.sym_box.ylo],
        point![snapshot.sym_box.xhi, snapshot.sym_box.yhi],
    );
    println!("{}", dump_final_path.to_string_lossy());
    let values = parse_dump_final(snapshot, &domain, cli.width, cli.zero_lvl);
    let values = cli.dump_final[1..]
        .iter()
        .try_fold(values, |values, dump_final_path| {
            let parent_dir = dump_final_path.parent().ok_or(anyhow!(
                "Unable to take parent directory: {}",
                dump_final_path.to_string_lossy()
            ))?;
            let _surface_coords_txt = parent_dir.join("surface_coords.txt");

            let dump_final = DumpFile::read(dump_final_path, &[])?;
            let snapshot = dump_final.get_snapshots()[0];
            println!("{}", dump_final_path.to_string_lossy());
            let values_new = parse_dump_final(snapshot, &domain, cli.width, cli.zero_lvl);

            plot_surface_2d(parent_dir, &values_new, &domain)?;
            anyhow::Ok::<DMatrix<f64>>(values + values_new)
        })?;
    if let Some(output_dir) = cli.output_dir {
        let values = values.scale(1.0 / cli.dump_final.len() as f64);
        plot_surface_2d(&output_dir, &values, &domain)?;
    }
    Ok(())
}

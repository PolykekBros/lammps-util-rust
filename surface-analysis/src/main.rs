use anyhow::Result;
use clap::Parser;
use colorgrad::preset::viridis;
use itertools::izip;
use lammps_util_rust::DumpFile;
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
    #[arg(value_name = "DUMP FINAL")]
    dump_final_file: PathBuf,

    #[arg(value_name = "OUTPUT DIR")]
    output_dir: PathBuf,

    #[arg(short, long, value_name = "SQUARE WIDTH (A)", default_value_t = 5.43 / 2.0)]
    width: f64,
}

fn get_surface_values(xyz: &[Point3<f64>], domain: &Domain, square_width: f64) -> DMatrix<f64> {
    assert!(square_width > 0.0);
    let x_count = (domain.width() / square_width).ceil() as usize + 1;
    let y_count = (domain.height() / square_width).ceil() as usize + 1;
    let mut values = DMatrix::repeat(y_count, x_count, f64::NEG_INFINITY);
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
            values[(x_i, y_i)] = values[(x_i, y_i)].max(p.z);
        });
    values
}

fn plot_surface_2d(output_dir: &Path, values: &DMatrix<f64>, domain: &Domain) -> Result<()> {
    let plot_path = output_dir.join("surface_2d.png");
    let drawing_area = BitMapBackend::new(&plot_path, (512, 512)).into_drawing_area();
    drawing_area.fill(&WHITE).unwrap();

    let pixel_width = drawing_area.dim_in_pixel().0;
    let (left, right) = drawing_area.split_horizontally(pixel_width - 60);
    let zero_lvl = 82.7813;
    let colorbar = Colorbar::new(zero_lvl - 20.0, zero_lvl + 10.0, viridis());
    colorbar.draw(ChartBuilder::on(&right));

    let mut chart_builder = ChartBuilder::on(&left);
    chart_builder.margin(10);
    heatmap(values, domain, &colorbar, chart_builder)?;
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let dump_final = DumpFile::read(&cli.dump_final_file, &[])?;
    let snapshot = dump_final.get_snapshots()[0];
    let snapshot_x = DVector::from_column_slice(snapshot.get_property("x"));
    let snapshot_y = DVector::from_column_slice(snapshot.get_property("y"));
    let snapshot_z = DVector::from_column_slice(snapshot.get_property("z"));
    let domain = Domain::new(
        point![snapshot_x.min(), snapshot_y.min()],
        point![snapshot_x.max(), snapshot_y.max()],
    );
    let xyz = izip![snapshot_x.iter(), snapshot_y.iter(), snapshot_z.iter()]
        .map(|(&x, &y, &z)| point![x, y, z])
        .collect::<Vec<_>>();
    let values = get_surface_values(&xyz, &domain, cli.width);
    plot_surface_2d(&cli.output_dir, &values, &domain)?;
    Ok(())
}

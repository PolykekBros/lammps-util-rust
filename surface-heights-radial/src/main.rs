use anyhow::{Context, Result, anyhow};
use clap::Parser;
use geomutil_util::{BoundingBox2, Point2, Point3};
use lammps_util_rust::IteratorAvg;
use std::{
    array,
    fs::File,
    io::{BufRead, BufReader},
    iter,
    path::{Path, PathBuf},
};

mod parser;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    surface_coords_path: PathBuf,
}

struct SurfaceData {
    bbox: BoundingBox2,
    data: Vec<Point3>,
}

impl SurfaceData {
    fn from_file(p: &Path) -> Result<SurfaceData> {
        let file = File::open(p).with_context(|| format!("file: {}:", p.display()))?;
        let mut token_iter = parser::Parser::new(BufReader::new(file).lines());
        let bbox = (|| {
            let bounds = token_iter
                .into_iter::<f32>()
                .take(4)
                .collect::<Result<Vec<_>, _>>()?;
            match bounds.len() {
                4 => Ok(BoundingBox2::from((
                    [bounds[0], bounds[1]].into(),
                    [bounds[2], bounds[3]].into(),
                ))),
                _ => Err(anyhow!("missing")),
            }
        })()
        .with_context(|| "parsing bounding box:")?;
        println!("bbox: {bbox:?}");
        let values = token_iter
            .into_iter::<f32>()
            .collect::<Result<Vec<_>, _>>()
            .with_context(|| "parsing surface data:")?;
        let dim = (values.len() as f32).sqrt().round() as usize;
        assert!(values.len() == dim * dim);
        println!("dims: {dim}x{dim}");
        let step = Point2::from([
            bbox.dimensions().x / dim as f32,
            bbox.dimensions().y / dim as f32,
        ]);
        let start = Point2::from([bbox.lower().x, bbox.lower().y]) + step / 2.0;
        let data = iter::zip(values, (0..dim).flat_map(|i| (0..dim).map(move |j| (i, j))))
            .map(|(z, (i, j))| {
                Point3::from([start.x + i as f32 * step.x, start.y + j as f32 * step.y, z])
            })
            .collect();
        Ok(SurfaceData { bbox, data })
    }

    fn radial_heights_distrib(&self, center: Point2) -> impl Iterator<Item = (usize, f32)> {
        const ANGLE_MAX: usize = 360;
        const ANGLE_STEP: usize = 10;
        const ANGLE_START: usize = ANGLE_STEP / 2;
        const ANGLE_N: usize = ANGLE_MAX / ANGLE_STEP;
        let mut sectors: [Vec<f32>; ANGLE_N] = array::from_fn(|_| Vec::default());
        self.data
            .iter()
            .filter(|point| point.z > 0.0)
            .for_each(|point| {
                let xy = Point2::from([point.x, point.y]) - center;
                let z = point.z;
                let i = ((xy.polar_angle().to_degrees() / ANGLE_STEP as f32).floor() as usize)
                    .clamp(0, ANGLE_N - 1);
                sectors[i].push(z);
            });
        let sectors: [f32; ANGLE_N] = array::from_fn(|i| {
            let sector = &sectors[i];
            sector.len() as f32 * sector.iter().copied().avg().unwrap_or_default()
        });
        sectors
            .into_iter()
            .enumerate()
            .map(|(i, x)| (ANGLE_START + i * ANGLE_STEP, x))
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let surface_data = SurfaceData::from_file(&cli.surface_coords_path)?;
    surface_data
        .radial_heights_distrib([0.0, 0.0].into())
        .for_each(|(angle, x)| println!("{angle}\t{x}"));
    Ok(())
}

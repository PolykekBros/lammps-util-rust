use anyhow::{Context, Result, anyhow, bail};
use clap::Parser;
use geomutil_util::{BoundingBox2, Point2};
use std::{
    fs::File,
    io::{BufRead, BufReader},
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
    data: Vec<Vec<f32>>,
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
                4 => Ok(([bounds[0], bounds[1]].into(), [bounds[2], bounds[3]].into()).into()),
                _ => Err(anyhow!("missing")),
            }
        })()
        .with_context(|| "parsing bounding box:")?;
        let values = token_iter
            .into_iter::<f32>()
            .collect::<Result<Vec<_>, _>>()
            .with_context(|| "parsing surface data:")?;
        let dim = (values.len() as f32).sqrt().round() as usize;
        let data = {
            let mut data = vec![Default::default(); dim];
            for (i, data) in data.iter_mut().enumerate() {
                *data = values
                    .get(i * dim..(i + 1) * dim)
                    .with_context(|| "missing surface data")?
                    .into();
            }
            data
        };
        Ok(SurfaceData { bbox, data })
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let surface_data = SurfaceData::from_file(&cli.surface_coords_path)?;
    Ok(())
}

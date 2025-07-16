use anyhow::{Context, Result, anyhow};
use clap::Parser;
use geomutil_util::Point2;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

mod parser;
use parser::TokenIterator;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    surface_coords_path: PathBuf,
}

struct BoundingBox {
    min: Point2,
    max: Point2,
}

impl BoundingBox {
    fn new(a: Point2, b: Point2) -> Self {
        Self {
            min: Point2::from([a.x.min(b.x), a.y.min(b.y)]),
            max: Point2::from([a.x.max(b.x), a.y.max(b.y)]),
        }
    }
}

struct SurfaceData {
    bbox: BoundingBox,
    data: Vec<Vec<f32>>,
}

impl SurfaceData {
    fn from_file(p: &Path) -> Result<SurfaceData> {
        let file = File::open(p).with_context(|| format!("file: {}:", p.to_string_lossy()))?;
        let token_iter = TokenIterator::new(BufReader::new(file).lines());
        Err(anyhow!(""))
    }
}

fn main() {
    let cli = Cli::parse();
    let surface_data = SurfaceData::from_file(&cli.surface_coords_path);
}

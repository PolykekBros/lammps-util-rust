mod clusterizer;
mod dump_file;
mod dump_snapshot;

pub use clusterizer::Clusterizer;
pub use dump_file::DumpFile;
pub use dump_snapshot::{DumpSnapshot, SymBox};

pub fn check_cutoff(a: (f64, f64, f64), b: (f64, f64, f64), cutoff: f64) -> bool {
    let d_x = a.0 - b.0;
    let d_y = a.1 - b.1;
    let d_z = a.2 - b.2;
    d_x * d_x + d_y * d_y + d_z * d_z <= cutoff * cutoff
}

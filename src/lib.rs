mod clusterizer;
mod dump_file;
mod dump_snapshot;
mod xyz;

pub use clusterizer::Clusterizer;
pub use dump_file::DumpFile;
pub use dump_snapshot::{DumpSnapshot, SymBox};
pub use xyz::{check_cutoff, XYZ};

pub fn range_f64(begin: f64, end: f64, count: usize) -> Vec<f64> {
    assert!(begin <= end);
    let step = (end - begin) / 1.max(count - 1) as f64;
    (0..count).map(|n| begin + n as f64 * step).collect()
}

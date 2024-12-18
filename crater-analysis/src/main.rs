use clap::Parser;
use core::error;
use lammps_util_rust::dump_file::DumpFile;
use lammps_util_rust::dump_snapshot::{DumpSnapshot, SymBox};
use std::collections::HashMap;
use std::io;
use std::path::Path;
use std::path::PathBuf;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[arg(value_name = "DUMP_INPUT")]
    dump_input_file: PathBuf,

    #[arg(value_name = "DUMP_FINAL")]
    dump_final_file: PathBuf,

    #[arg(value_name = "OUTPUT_DIR")]
    output_dir: PathBuf,

    #[arg(short, long, value_name = "MAX_DEPTH (A)", default_value_t = 50.0)]
    max_depth: f64,

    #[arg(short, long, value_name = "CUTOFF (A)", default_value_t = 3.0)]
    cutoff: f64,

    #[arg(short, long, value_name = "STRIPE_WIDTH (A)", default_value_t = 5.43 / 2.0)]
    width: f64,
}

struct Stripes<'a> {
    indices: Vec<Vec<usize>>,
    sym_box: SymBox,
    x: &'a [f64],
    y: &'a [f64],
    z: &'a [f64],
}

fn check_cutoff(a: (f64, f64, f64), b: (f64, f64, f64), cutoff: f64) -> bool {
    let d_x = a.0 - b.0;
    let d_y = a.1 - b.1;
    let d_z = a.2 - b.2;
    d_x * d_x + d_y * d_y + d_z * d_z <= cutoff * cutoff
}

impl<'a> Stripes<'a> {
    pub fn new(snap: &'a DumpSnapshot, zero_lvl: f64, max_depth: f64, width: f64) -> Self {
        assert!(max_depth >= 0.0);
        assert!(width >= 0.0);
        let x = snap.get_property("x");
        let y = snap.get_property("y");
        let z = snap.get_property("z");
        let count = (max_depth / width).ceil() as usize;
        let indices = (0..count).fold(Vec::new(), |mut indices, _| {
            indices.push(Vec::new());
            indices
        });
        let indices = z
            .iter()
            .copied()
            .enumerate()
            .fold(indices, |mut indices, (i, z)| {
                if zero_lvl - max_depth < z && z <= zero_lvl {
                    let j = ((zero_lvl - z) / width).floor() as usize;
                    indices[j].push(i);
                }
                indices
            });
        Self {
            indices,
            sym_box: snap.sym_box.clone(),
            x,
            y,
            z,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.indices.len()
    }

    #[inline]
    pub fn get_xyz(&self, i: usize) -> (f64, f64, f64) {
        (self.x[i], self.y[i], self.z[i])
    }

    pub fn get_missing_indices(&self, snap: &Stripes, i: usize, cutoff: f64) -> Vec<usize> {
        self.indices[i]
            .iter()
            .copied()
            .fold(Vec::new(), |mut indices, self_i| {
                let self_xyz = self.get_xyz(self_i);
                if self.sym_box.xlo + cutoff < self_xyz.0
                    && self_xyz.0 < self.sym_box.xhi - cutoff
                    && self.sym_box.ylo + cutoff < self_xyz.1
                    && self_xyz.1 < self.sym_box.yhi - cutoff
                {
                    let missing = snap.indices[i]
                        .iter()
                        .copied()
                        .fold(true, |missing, snap_i| {
                            let snap_xyz = snap.get_xyz(snap_i);
                            missing && !check_cutoff(self_xyz, snap_xyz, cutoff)
                        });
                    if missing {
                        indices.push(self_i);
                    }
                }
                indices
            })
    }
}

fn get_crater_candidate_indecies(
    snap_input: &DumpSnapshot,
    snap_final: &DumpSnapshot,
    max_depth: f64,
    cutoff: f64,
    width: f64,
) -> Vec<usize> {
    let zero_lvl = snap_input.get_zero_lvl();
    let stripes_input = Stripes::new(snap_input, zero_lvl, max_depth, width);
    let stripes_final = Stripes::new(snap_final, zero_lvl, max_depth, width);
    (0..stripes_input.len()).fold(Vec::new(), |mut indices, i| {
        let stripe_indices = stripes_input.get_missing_indices(&stripes_final, i, cutoff);
        indices.extend(stripe_indices);
        indices
    })
}

fn save_crater_candidates(
    snap_input: &DumpSnapshot,
    ids: &[usize],
    output_path: &Path,
) -> io::Result<()> {
    let keys: HashMap<String, usize> = snap_input
        .get_keys()
        .into_iter()
        .enumerate()
        .map(|(i, key)| (key.to_string(), i))
        .collect();
    let mut snap_output =
        DumpSnapshot::new(keys, snap_input.step, ids.len(), snap_input.sym_box.clone());
    for (new_i, i) in ids.iter().copied().enumerate() {
        for (j, _) in snap_input.get_keys().iter().enumerate() {
            snap_output.set_atom_value(j, new_i, snap_input.get_atom_value(j, i));
        }
    }
    let dump_file = DumpFile::new(vec![snap_output]);
    dump_file.save(output_path)
}

fn main() -> Result<(), Box<dyn error::Error>> {
    let cli = Cli::parse();
    let dump_input = DumpFile::read(&cli.dump_input_file, &[])?;
    let snapshot_input = dump_input.get_snapshots()[0];
    let dump_final = DumpFile::read(&cli.dump_final_file, &[])?;
    let snapshot_final = dump_final.get_snapshots()[0];
    let candidate_indicies = get_crater_candidate_indecies(
        snapshot_input,
        snapshot_final,
        cli.max_depth,
        cli.cutoff,
        cli.width,
    );
    save_crater_candidates(
        snapshot_input,
        &candidate_indicies,
        &cli.output_dir.join("dump.crater_candidates"),
    )?;
    // println!("{info}");
    Ok(())
}

// let crater_count = crater_indexes.len();
// let mut surface_count = 0;
// let mut z_avg = 0.0;
// let mut z_min = f64::INFINITY;
// for i in crater_indexes {
//     let z = snap_input.get_property("z")[i];
//     if z > -2.4 * 0.707 + zero_lvl {
//         surface_count += 1;
//     }
//     z_min = z_min.min(z - zero_lvl);
//     z_avg += z - zero_lvl;
// }
// z_avg /= crater_count as f64;
// let volume = crater_count as f64 * 20.1;
// let surface = surface_count as f64 * 7.3712;
// format!("{crater_count} {volume} {surface} {z_avg} {z_min}")

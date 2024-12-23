use clap::Parser;
use core::error;
use lammps_util_rust::{check_cutoff, Clusterizer, DumpFile, DumpSnapshot, SymBox, XYZ};
use std::{collections::HashMap, path::PathBuf};

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
    pub fn get_xyz(&self, i: usize) -> XYZ {
        XYZ::from([self.x[i], self.y[i], self.z[i]])
    }

    pub fn get_missing_indices(&self, snap: &Stripes, i: usize, cutoff: f64) -> Vec<usize> {
        self.indices[i]
            .iter()
            .copied()
            .fold(Vec::new(), |mut indices, self_i| {
                let self_xyz = self.get_xyz(self_i);
                if self.sym_box.xlo + cutoff < self_xyz.x()
                    && self_xyz.x() < self.sym_box.xhi - cutoff
                    && self.sym_box.ylo + cutoff < self_xyz.y()
                    && self_xyz.y() < self.sym_box.yhi - cutoff
                {
                    let neighbours =
                        snap.indices[i]
                            .iter()
                            .copied()
                            .fold(0, |neighbours, snap_i| {
                                let snap_xyz = snap.get_xyz(snap_i);
                                match check_cutoff(self_xyz, snap_xyz, cutoff) {
                                    true => neighbours + 1,
                                    false => neighbours,
                                }
                            });
                    if neighbours == 0 {
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

fn crater_candidates(snap_input: &DumpSnapshot, ids: &[usize]) -> DumpFile {
    let keys = snap_input.get_keys_map();
    let snap_output = ids.iter().copied().enumerate().fold(
        DumpSnapshot::new(
            keys.clone(),
            snap_input.step,
            ids.len(),
            snap_input.sym_box.clone(),
        ),
        |mut snap_output, (new_i, i)| {
            snap_input.get_keys().iter().enumerate().for_each(|(j, _)| {
                snap_output.set_atom_value(j, new_i, snap_input.get_atom_value(j, i));
            });
            snap_output
        },
    );
    DumpFile::new(vec![snap_output])
}

fn get_crater_info(snapshot: &DumpSnapshot, zero_lvl: f64) -> String {
    let cluster = snapshot.get_property("cluster");
    let cluster_cnt = cluster
        .iter()
        .copied()
        .fold(HashMap::new(), |mut indices, cluster| {
            indices
                .entry(cluster as u64)
                .and_modify(|cnt| *cnt += 1)
                .or_insert(1);
            indices
        });
    let (cluster_max, _) = cluster_cnt.iter().fold(
        (0, u64::MIN),
        |(cluster_max, max_cnt), (cluster, cnt)| match *cnt > max_cnt {
            true => (*cluster, *cnt),
            false => (cluster_max, max_cnt),
        },
    );
    let z = snapshot.get_property("z");
    let mut crater_count = 0;
    let mut surface_count = 0;
    let mut z_avg = 0.0;
    let mut z_min = f64::INFINITY;
    for i in 0..snapshot.atoms_count {
        let z = z[i];
        let cluster = cluster[i];
        if cluster as u64 == cluster_max {
            if z > -2.4 * 0.707 + zero_lvl {
                surface_count += 1;
            }
            crater_count += 1;
            z_min = z_min.min(z - zero_lvl);
            z_avg += z - zero_lvl;
        }
    }
    z_avg /= crater_count as f64;
    let volume = crater_count as f64 * 20.1;
    let surface = surface_count as f64 * 7.3712;
    format!("{crater_count} {volume} {surface} {z_avg} {z_min}")
}

fn main() -> Result<(), Box<dyn error::Error>> {
    let cli = Cli::parse();
    let dump_input = DumpFile::read(&cli.dump_input_file, &[])?;
    let snapshot_input = dump_input.get_snapshots()[0];
    let zero_lvl = snapshot_input.get_zero_lvl();
    let dump_final = DumpFile::read(&cli.dump_final_file, &[])?;
    let snapshot_final = dump_final.get_snapshots()[0];
    let candidate_indicies = get_crater_candidate_indecies(
        snapshot_input,
        snapshot_final,
        cli.max_depth,
        cli.cutoff,
        cli.width,
    );
    let dump_crater = crater_candidates(snapshot_input, &candidate_indicies);
    let dump_crater_clusterized =
        Clusterizer::new(dump_crater.get_snapshots()[0], 3.0).clusterize();
    dump_crater_clusterized.save(&cli.output_dir.join("dump.crater_candidates"))?;
    let info = get_crater_info(dump_crater_clusterized.get_snapshots()[0], zero_lvl);
    println!("{info}");
    Ok(())
}

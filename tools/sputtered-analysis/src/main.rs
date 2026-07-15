#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]
use anyhow::Result;
use clap::Parser;
use clap::Subcommand;
use lammps_util::{DumpFile, DumpSnapshot, MainWrapper, RunDir, Task, get_clusters, process_results_dir};
use std::{collections::HashMap, iter, path::PathBuf};

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Default analysis: count atoms by type per run
    Default(DefaultArgs),
    /// Print cluster element composition counts (e.g. Si1O2)
    ClusterComposition(ClusterCompositionArgs),
}

#[derive(Parser)]
struct DefaultArgs {
    /// Atom types "<type 1>,<type 2>,...,<type N>", ex. "Si,C,O"
    #[arg(short, long)]
    particles: String,

    /// Dump file name (default: dump.sputter)
    #[arg(long, default_value = "dump.sputter")]
    dump_file: String,

    /// Results directories (comma-separated or --dirs repeated)
    #[arg(short, long, value_delimiter = ',')]
    dirs: Vec<PathBuf>,

    /// Path to a single dump file (overrides --dirs)
    #[arg(long)]
    dump_path: Option<PathBuf>,
}

#[derive(Parser)]
struct ClusterCompositionArgs {
    /// Element to type mapping, ex. "Si:1,O:2"
    #[arg(long)]
    elements: String,

    /// Dump file name (default: dump.sputter)
    #[arg(long, default_value = "dump.sputter")]
    dump_file: String,

    /// Results directories (comma-separated or --dirs repeated)
    #[arg(short, long, value_delimiter = ',')]
    dirs: Vec<PathBuf>,

    /// Path to a single dump file (overrides --dirs)
    #[arg(long)]
    dump_path: Option<PathBuf>,
}

struct Atom {
    _id: usize,
    atype: usize,
    _coords: [f64; 3],
    velocity: [f64; 3],
    mass: f64,
}

impl Atom {
    const fn new(id: usize, atype: usize, coords: [f64; 3], velocity: [f64; 3], mass: f64) -> Self {
        Self {
            _id: id,
            atype,
            _coords: coords,
            velocity,
            mass,
        }
    }
}

#[derive(Default)]
struct Cluster {
    mass: f64,
    counts: Vec<usize>,
    momentum: [f64; 3],
    ek: f64,
    angle: f64,
}

impl Cluster {
    fn new(atoms: &[Atom], types_map: &HashMap<usize, String>) -> Self {
        let mut counts = vec![0; types_map.len()];
        let mut mass = 0.0;
        let mut momentum = [0.0; 3];
        for atom in atoms {
            mass += atom.mass;
            iter::zip(&mut momentum, atom.velocity).for_each(|(m, v)| *m += v * atom.mass);
            counts[atom.atype - 1] += 1;
        }
        let ek =
            2.0 * 5.1875 * 1e-5 * momentum.iter().map(|m| m.powi(2)).sum::<f64>() / (2.0 * mass);
        let angle = (momentum[2] / momentum[0].hypot(momentum[1])).atan();
        Self {
            mass,
            counts,
            momentum,
            ek,
            angle,
        }
    }
}

fn parse_types(s: &str) -> (HashMap<usize, String>, Vec<String>) {
    let type_names = s
        .split(',')
        .map(|s| s.trim().to_string())
        .collect::<Vec<_>>();
    let types_map = type_names.iter().map(String::from).enumerate().collect();
    (types_map, type_names)
}

fn parse_elements(s: &str) -> HashMap<usize, String> {
    s.split(',')
        .map(|pair| {
            let mut parts = pair.split(':');
            let name = parts.next().unwrap().trim().to_string();
            let type_id = parts.next().unwrap().trim().parse::<usize>().expect("invalid type id");
            (type_id, name)
        })
        .collect()
}

fn cluster_composition_string(
    counts: &[usize],
    elements_map: &HashMap<usize, String>,
) -> String {
    let mut parts = Vec::new();
    for (type_id, name) in elements_map.iter() {
        let idx = type_id - 1;
        if idx < counts.len() && counts[idx] > 0 {
            parts.push(format!("{name}{}", counts[idx]));
        }
    }
    parts.join("")
}

fn analyze_clusters(dump: &DumpSnapshot, types_map: &HashMap<usize, String>) -> Vec<Cluster> {
    let id = dump.get_property("id");
    let atype = dump.get_property("type");
    let x = dump.get_property("x");
    let y = dump.get_property("y");
    let z = dump.get_property("z");
    let vx = dump.get_property("vx");
    let vy = dump.get_property("vy");
    let vz = dump.get_property("vz");
    let mass = dump.get_property("mass");
    get_clusters(dump)
        .into_values()
        .map(|atoms_idx| {
            let atoms = atoms_idx
                .into_iter()
                .map(|atom_idx| {
                    Atom::new(
                        id[atom_idx] as usize,
                        atype[atom_idx] as usize,
                        [x[atom_idx], y[atom_idx], z[atom_idx]],
                        [vx[atom_idx], vy[atom_idx], vz[atom_idx]],
                        mass[atom_idx],
                    )
                })
                .collect::<Vec<_>>();
            Cluster::new(&atoms, types_map)
        })
        .collect()
}

fn do_single_dir(
    dir: &RunDir,
    types_map: &HashMap<usize, String>,
    dump_file: &str,
) -> Result<Vec<Cluster>> {
    let dump = DumpFile::read(dir.path.join(dump_file), &[])?;
    let clusters = analyze_clusters(&dump.get_snapshots()[0], types_map);
    Ok(clusters)
}

fn do_single_file(
    path: &PathBuf,
    types_map: &HashMap<usize, String>,
) -> Result<Vec<Cluster>> {
    let dump = DumpFile::read(path, &[])?;
    let clusters = analyze_clusters(&dump.get_snapshots()[0], types_map);
    Ok(clusters)
}

struct SputteredAnalysisTask {
    types_map: HashMap<usize, String>,
    type_names: Vec<String>,
    dump_file: String,
    dirs: Vec<PathBuf>,
    dump_path: Option<PathBuf>,
}

impl SputteredAnalysisTask {
    fn new(
        types_map: HashMap<usize, String>,
        type_names: Vec<String>,
        dump_file: String,
        dirs: Vec<PathBuf>,
        dump_path: Option<PathBuf>,
    ) -> Self {
        Self {
            types_map,
            type_names,
            dump_file,
            dirs,
            dump_path,
        }
    }
}

impl Task for SputteredAnalysisTask {
    type Output = ();
    fn run(&self) -> Result<Self::Output> {
        let all_clusters = if let Some(ref path) = self.dump_path {
            let clusters = do_single_file(path, &self.types_map)?;
            vec![(RunDir { path: path.clone(), num: 0 }, clusters)]
        } else if self.dirs.is_empty() {
            anyhow::bail!("no directories specified. Use --dirs or --dump-path.");
        } else {
            let mut results = Vec::new();
            for dir in &self.dirs {
                match process_results_dir(dir, |run_dir| {
                    let result = do_single_dir(run_dir, &self.types_map, &self.dump_file);
                    eprintln!(
                        "{} {:?}",
                        run_dir.path.to_string_lossy(),
                        &result.as_ref().map(Vec::len)
                    );
                    result
                }) {
                    Ok(r) => {
                        if r.is_empty() {
                            eprintln!(
                                "warning: no numeric run directories found in '{}'",
                                dir.display()
                            );
                        }
                        results.extend(r);
                    }
                    Err(e) => {
                        eprintln!("error processing '{}': {e}", dir.display());
                        return Err(e);
                    }
                }
            }
            results
        };

        println!("# № {} ∑", self.type_names.join(" "));
        for (run_dir, clusters) in all_clusters {
            let counts = clusters.into_iter().map(|cluster| cluster.counts).fold(
                vec![0; self.types_map.len()],
                |mut agg, counts| {
                    iter::zip(&mut agg, counts).for_each(|(a, b)| *a += b);
                    agg
                },
            );
            let sum = counts.iter().sum::<usize>();
            print!("{}\t", run_dir.num);
            for count in counts {
                print!("{count}\t");
            }
            println!("{sum}");
        }
        Ok(())
    }
}

struct ClusterCompositionTask {
    elements_map: HashMap<usize, String>,
    dump_file: String,
    dirs: Vec<PathBuf>,
    dump_path: Option<PathBuf>,
}

impl ClusterCompositionTask {
    fn new(
        elements_map: HashMap<usize, String>,
        dump_file: String,
        dirs: Vec<PathBuf>,
        dump_path: Option<PathBuf>,
    ) -> Self {
        Self {
            elements_map,
            dump_file,
            dirs,
            dump_path,
        }
    }
}

impl Task for ClusterCompositionTask {
    type Output = ();
    fn run(&self) -> Result<Self::Output> {
        let per_dir_counts: Vec<HashMap<String, usize>> = if let Some(ref path) = self.dump_path {
            let clusters = do_single_file(path, &self.elements_map)?;
            let mut counts = HashMap::new();
            for cluster in clusters {
                let comp = cluster_composition_string(&cluster.counts, &self.elements_map);
                *counts.entry(comp).or_default() += 1;
            }
            vec![counts]
        } else if self.dirs.is_empty() {
            anyhow::bail!("no directories specified. Use --dirs or --dump-path.");
        } else {
            let mut per_dir = Vec::new();
            for dir in &self.dirs {
                match process_results_dir(dir, |run_dir| {
                    let result = do_single_dir(run_dir, &self.elements_map, &self.dump_file);
                    eprintln!(
                        "{} {:?}",
                        run_dir.path.to_string_lossy(),
                        &result.as_ref().map(Vec::len)
                    );
                    result
                }) {
                    Ok(r) => {
                        if r.is_empty() {
                            eprintln!(
                                "warning: no numeric run directories found in '{}'",
                                dir.display()
                            );
                            continue;
                        }
                        let mut counts = HashMap::new();
                        for (_run_dir, clusters) in r {
                            for cluster in clusters {
                                let comp = cluster_composition_string(&cluster.counts, &self.elements_map);
                                *counts.entry(comp).or_default() += 1;
                            }
                        }
                        per_dir.push(counts);
                    }
                    Err(e) => {
                        eprintln!("error processing '{}': {e}", dir.display());
                        return Err(e);
                    }
                }
            }
            per_dir
        };

        if per_dir_counts.is_empty() {
            anyhow::bail!("no clusters found in any directory");
        }

        let n_dirs = per_dir_counts.len();
        let mut total_counts: HashMap<String, usize> = HashMap::new();
        for dir_counts in &per_dir_counts {
            for (comp, count) in dir_counts {
                *total_counts.entry(comp.clone()).or_default() += count;
            }
        }

        let mut entries: Vec<_> = total_counts.into_iter().collect();
        entries.sort_by(|(a, _), (b, _)| a.cmp(b));

        println!("# Cluster Composition Counts (averaged across {n_dirs} dir(s))");
        println!("# Composition\tTotal\tAvg");
        for (comp, total) in entries {
            let avg = total as f64 / n_dirs as f64;
            println!("{comp}\t{total}\t{avg:.2}");
        }

        Ok(())
    }
}

fn main() -> Result<()> {
    MainWrapper::<Cli>::default().run(|cli| match cli.command {
        Commands::Default(args) => {
            let (types_map, type_names) = parse_types(&args.particles);
            Box::new(SputteredAnalysisTask::new(
                types_map,
                type_names,
                args.dump_file,
                args.dirs,
                args.dump_path,
            ))
        }
        Commands::ClusterComposition(args) => {
            let elements_map = parse_elements(&args.elements);
            Box::new(ClusterCompositionTask::new(
                elements_map,
                args.dump_file,
                args.dirs,
                args.dump_path,
            ))
        }
    })
}

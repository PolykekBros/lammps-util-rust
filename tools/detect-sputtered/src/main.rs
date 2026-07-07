#![deny(
    clippy::all,
    clippy::correctness,
    clippy::suspicious,
    clippy::complexity,
    clippy::perf,
    clippy::style,
    clippy::pedantic
)]

use anyhow::Result;
use clap::{Args, Parser, Subcommand};
use lammps_files::{Dump, Snapshot};
use lammps_util::{
    clusterize_snapshot, copy_snapshot_with_indices, get_clusters, process_results_dir,
};
use log::info;
use std::marker::PhantomData;
use std::{
    collections::HashSet,
    path::{Path, PathBuf},
};

const DUMP_FINAL: &str = "dump.final";
const DUMP_SPUTTERED: &str = "dump.sputtered";
const DUMP_NON_SPUTTERED: &str = "dump.non_sputtered";

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Sputtered cluster threshold (Count)
    #[arg(short, long, default_value_t = 1000)]
    threshold: usize,

    /// Cluster neighbors cutoff (A)
    #[arg(short, long, default_value_t = 3.0)]
    cutoff: f64,
}

#[derive(Subcommand)]
enum Commands {
    /// Crater analysis for a single run dir
    Single(SingleCMD),

    /// Crater analysis for the whole results folder
    Multi(MultiCMD),
}

#[derive(Args)]
struct SingleCMD {
    run_dir: PathBuf,
}

#[derive(Args)]
struct MultiCMD {
    results_dir: PathBuf,
}

trait DirTask: Clone {
    type Output;
    fn run<P: AsRef<Path>>(&self, dir: P) -> Result<Self::Output>;
}

trait Task {
    type Output;
    fn run(&self) -> Result<Self::Output>;
}

#[derive(Debug, Clone)]
struct DetectSputteredTask {
    threshold: usize,
    cutoff: f64,
}

struct DetectSputteredTaskResult {
    pub snapshot_sputtered: Snapshot,
    pub snapshot_non_sputtered: Snapshot,
}

impl DetectSputteredTask {
    fn new(threshold: usize, cutoff: f64) -> Self {
        Self { threshold, cutoff }
    }
}

impl DirTask for DetectSputteredTask {
    type Output = DetectSputteredTaskResult;
    fn run<P: AsRef<Path>>(&self, dir: P) -> Result<Self::Output> {
        use std::ops::Sub;

        let dir = dir.as_ref();
        let dump_final = Dump::open(dir.join(DUMP_FINAL), &[])?;
        let snapshot_final = clusterize_snapshot(&dump_final.get_snapshots()[0], self.cutoff);
        let clusters = get_clusters(&snapshot_final);
        let indices: HashSet<usize> = clusters
            .values()
            .flat_map(|ids| ids.iter())
            .copied()
            .collect();
        let sputtered_indices: HashSet<usize> = clusters
            .values()
            .filter(|ids| ids.len() < self.threshold)
            .flat_map(|ids| ids.iter())
            .copied()
            .collect();
        let non_sputtered_indices: HashSet<usize> = indices.sub(&sputtered_indices);
        let snapshot_sputtered =
            copy_snapshot_with_indices(&snapshot_final, sputtered_indices.into_iter());
        let snapshot_non_sputtered =
            copy_snapshot_with_indices(&snapshot_final, non_sputtered_indices.into_iter());
        Ok(DetectSputteredTaskResult {
            snapshot_sputtered,
            snapshot_non_sputtered,
        })
    }
}

#[derive(Debug, Clone)]
struct DetectSputteredAndSaveTask<T> {
    parent: T,
}

impl<T> DetectSputteredAndSaveTask<T> {
    fn new(parent: T) -> Self {
        Self { parent }
    }
}

impl<T> DirTask for DetectSputteredAndSaveTask<T>
where
    T: DirTask<Output = DetectSputteredTaskResult>,
{
    type Output = ();
    fn run<P: AsRef<Path>>(&self, dir: P) -> Result<Self::Output> {
        let result = self.parent.run(&dir)?;
        let dir = dir.as_ref();
        let count = result.snapshot_sputtered.get_atoms_count();
        Dump::new(vec![result.snapshot_sputtered]).save(dir.join(DUMP_SPUTTERED))?;
        Dump::new(vec![result.snapshot_non_sputtered]).save(dir.join(DUMP_NON_SPUTTERED))?;
        info!("{} | {count}", dir.display());
        Ok(())
    }
}

struct SingleTask<T> {
    run_dir: PathBuf,
    task: T,
}

impl<T> SingleTask<T> {
    fn new(run_dir: PathBuf, task: T) -> Self {
        Self { run_dir, task }
    }
}

impl<T> Task for SingleTask<T>
where
    T: DirTask<Output = DetectSputteredTaskResult>,
{
    type Output = ();
    fn run(&self) -> Result<Self::Output> {
        DetectSputteredAndSaveTask::new(self.task.clone()).run(&self.run_dir)
    }
}

struct MultiTask<T> {
    results_dir: PathBuf,
    task: T,
}

impl<T> MultiTask<T> {
    fn new(results_dir: PathBuf, task: T) -> Self {
        Self { results_dir, task }
    }
}

impl<T> Task for MultiTask<T>
where
    T: DirTask<Output = DetectSputteredTaskResult> + Send + Sync,
{
    type Output = ();
    fn run(&self) -> Result<Self::Output> {
        process_results_dir(&self.results_dir, |run_dir| {
            SingleTask::new(run_dir.path.clone(), self.task.clone()).run()
        })
        .map(|_| ())
    }
}

struct MainWrapper<T> {
    cli: PhantomData<T>,
}

impl<T> Default for MainWrapper<T> {
    fn default() -> Self {
        Self { cli: PhantomData }
    }
}

#[allow(clippy::unused_self)]
impl<T: Parser> MainWrapper<T> {
    fn run<F: FnMut(T) -> Box<dyn Task<Output = ()>>>(self, f: F) -> Result<()> {
        env_logger::init();
        let cli = T::parse();
        let mut f = f;
        f(cli).run()
    }
}

fn main() -> Result<()> {
    MainWrapper::<Cli>::default().run(|cli| {
        let task = DetectSputteredTask::new(cli.threshold, cli.cutoff);
        match cli.command {
            Commands::Single(args) => Box::new(SingleTask::new(args.run_dir, task)),
            Commands::Multi(args) => Box::new(MultiTask::new(args.results_dir, task)),
        }
    })
}

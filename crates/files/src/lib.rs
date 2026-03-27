#![deny(
    clippy::all,
    clippy::correctness,
    clippy::suspicious,
    clippy::complexity,
    clippy::perf,
    clippy::style,
    clippy::pedantic
)]

pub mod dump;
pub mod error;
mod parser;
pub mod snapshot;

pub use dump::Dump;
pub use snapshot::Snapshot;

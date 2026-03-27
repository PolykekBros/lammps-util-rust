use std::{
    io::BufRead,
    collections::HashMap,
};

use geomutil_util::BoundingBox3;

use crate::{
    error::{Result, Error},
    parser::{Parser},
};

// Example dump file
// ITEM: TIMESTEP
// 100
// ITEM: NUMBER OF ATOMS
// 5
// ITEM: BOX BOUNDS pp pp pp
// 0.0000000000000000e+00 1.0000000000000000e+01
// 0.0000000000000000e+00 1.0000000000000000e+01
// 0.0000000000000000e+00 1.0000000000000000e+01
// ITEM: ATOMS id type x y z
// 1 1 1.25000 1.25000 1.25000
// 2 1 3.75000 1.25000 1.25000
// 3 2 1.25000 3.75000 1.25000
// 4 2 3.75000 3.75000 1.25000
// 5 1 5.00000 5.00000 5.00000

pub const HEADER_TIMESTEP: &str = "ITEM: TIMESTEP";
pub const HEADER_NUM_OF_ATOMS: &str = "ITEM: NUMBER OF ATOMS";
pub const HEADER_SYM_BOX: &str = "ITEM: BOX BOUNDS";
pub const HEADER_ATOMS: &str = "ITEM: ATOMS";

#[derive(Debug, Clone)]
pub struct SnapshotMeta {
    pub timestep: u64,
    pub atoms_count: usize,
    pub sym_box: SymBox,
    keys: HashMap<String, usize>,
}

fn parse_timestep() {
    todo!()
}

fn parse_atom_count() {
    todo!()
}

impl SnapshotMeta {
    fn parse<B: BufRead>(parser: &mut Parser<B>) -> Result<Self>
    {
        let line = parser.next().ok_or(Error::Timestep(parser.current()))??;
        if line.trim() != HEADER_TIMESTEP {
            return Err(Error::Timestep(parser.current()));
        }
        let timestep = parser.next().ok_or(Error::Timestep(parser.current()))??
            .trim().parse::<u64>().map_err(|_| Error::Timestep(parser.current()))?;
        let line = parser.next().ok_or(Error::NumOfAtoms(parser.current()))??;
        if line.trim() != HEADER_NUM_OF_ATOMS {
            return Err(Error::NumOfAtoms(parser.current()));
        }
        let atoms_count = parser.next().ok_or(Error::NumOfAtoms(parser.current()))??
            .trim().parse::<usize>().map_err(|_| Error::NumOfAtoms(parser.current()))?;

        todo!()
    }
}

pub enum Boundaries {
    
}

#[derive(Debug, Clone)]
pub struct SymBox {
    pub boundaries: String,
    pub bbox: BoundingBox3<f64>,
}

impl SymBox {
    fn parse() {
        todo()!
    }
}

pub struct Snapshot {
    meta, atoms
};


use geomutil_util::Point3;
use kd_tree::KdPoint;
use std::ops::{Deref, DerefMut};

use crate::SymBox;

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub struct XYZ {
    pub coords: Point3,
    pub index: usize,
}

impl Deref for XYZ {
    type Target = Point3;
    fn deref(&self) -> &Self::Target {
        &self.coords
    }
}

impl DerefMut for XYZ {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.coords
    }
}

impl KdPoint for XYZ {
    type Scalar = f32;
    type Dim = typenum::U3;
    fn at(&self, i: usize) -> f32 {
        self.coords.coords[i]
    }
}

impl XYZ {
    pub fn from(coords: impl Into<Point3>, index: usize) -> Self {
        Self {
            coords: coords.into(),
            index,
        }
    }

    pub fn check_cutoff(a: Self, b: Self, cutoff: f32) -> bool {
        a.distance_squared(*b) <= cutoff * cutoff
    }

    pub fn get_supercell_coords(coords: &mut Vec<Self>, sym_box: &SymBox, cutoff: f32) {
        let lo = sym_box.bbox.lower();
        let hi = sym_box.bbox.upper();
        let shift = sym_box.bbox.dimensions();
        let periods_and_shifts = (-1..=1)
            .flat_map(|px| (-1..=1).map(move |py| (px, py)))
            .flat_map(|(px, py)| (-1..=1).map(move |pz| [px, py, pz]))
            .filter(|periods| periods.iter().any(|&period| period != 0))
            .map(|periods| Point3::from(periods.map(|p| p as f32)))
            .map(|periods| (periods, shift * periods))
            .collect::<Vec<_>>();
        for atom_idx in 0..coords.len() {
            for (period, shift) in &periods_and_shifts {
                if (0..3).all(|i| match period[i] {
                    1.0 => coords[atom_idx][i] < lo[i] + cutoff,
                    -1.0 => coords[atom_idx][i] > hi[i] - cutoff,
                    _ => true,
                }) {
                    coords.push(XYZ::from(*coords[atom_idx] + *shift, coords.len()));
                }
            }
        }
    }
}

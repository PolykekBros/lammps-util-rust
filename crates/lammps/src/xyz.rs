use geomutil_util::Point3;
use kd_tree::KdPoint;
use std::{
    array,
    ops::{Deref, DerefMut},
};

use crate::SymBox;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XYZ {
    pub coords: Point3<f64>,
    pub index: usize,
    pub is_ghost: bool,
}

impl Deref for XYZ {
    type Target = Point3<f64>;
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
    type Scalar = f64;
    type Dim = typenum::U3;
    fn at(&self, i: usize) -> f64 {
        self.coords[i]
    }
}

impl XYZ {
    #[must_use]
    pub const fn new(coords: Point3<f64>, index: usize, is_ghost: bool) -> Self {
        Self {
            coords,
            index,
            is_ghost,
        }
    }

    #[must_use]
    pub fn check_cutoff(a: Self, b: Self, cutoff: f64) -> bool {
        a.distance_squared(*b) <= cutoff * cutoff
    }

    pub fn get_supercell_coords(coords: &mut Vec<Self>, sym_box: &SymBox, cutoff: f64) {
        let boundaries = sym_box.boundaries.split_whitespace().collect::<Vec<_>>();
        let periods: [_; 3] = array::from_fn(|i| match boundaries[i] {
            "pp" => vec![-1, 0, 1],
            _ => vec![0],
        });
        assert_eq!(boundaries.len(), 3);
        let lo = sym_box.bbox.lower;
        let hi = sym_box.bbox.upper;
        let shift = sym_box.bbox.dimensions();
        let periods_and_shifts = periods[0]
            .iter()
            .copied()
            .flat_map(|px| periods[1].iter().copied().map(move |py| (px, py)))
            .flat_map(|(px, py)| periods[2].iter().copied().map(move |pz| [px, py, pz]))
            .filter(|periods| periods.iter().any(|&period| period != 0))
            .map(|periods| Point3::from(periods.map(f64::from)))
            .map(|periods| (periods, shift * periods))
            .collect::<Vec<_>>();
        for atom_idx in 0..coords.len() {
            for (period, shift) in &periods_and_shifts {
                if (0..3).all(|i| match period[i] {
                    1.0 => coords[atom_idx][i] < lo[i] + cutoff,
                    -1.0 => coords[atom_idx][i] > hi[i] - cutoff,
                    _ => true,
                }) {
                    coords.push(Self::new(
                        *coords[atom_idx] + *shift,
                        coords[atom_idx].index,
                        true,
                    ));
                }
            }
        }
    }
}

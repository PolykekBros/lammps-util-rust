use kd_tree::KdPoint;
use nalgebra::Point3;
use std::{
    hash::{Hash, Hasher},
    ops::{Deref, DerefMut},
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XYZ(Point3<f64>, usize);

impl Eq for XYZ {}

impl Hash for XYZ {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.coords.iter().for_each(|n| {
            n.to_bits().hash(state);
        });
    }
}

impl Deref for XYZ {
    type Target = Point3<f64>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for XYZ {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
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
    pub fn from(xyz: [f64; 3], i: usize) -> Self {
        Self(xyz.into(), i)
    }

    pub fn index(&self) -> usize {
        self.1
    }
}

pub fn check_cutoff(a: XYZ, b: XYZ, cutoff: f64) -> bool {
    nalgebra::distance_squared(&a, &b) <= cutoff * cutoff
}

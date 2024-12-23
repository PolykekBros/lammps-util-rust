use kd_tree::KdPoint;
use std::hash::{Hash, Hasher};

#[derive(Debug, Clone, Copy)]
pub struct XYZ {
    coords: [f64; 3],
}

impl XYZ {
    pub fn from(coords: [f64; 3]) -> Self {
        Self { coords }
    }

    pub fn x(&self) -> f64 {
        self.coords[0]
    }
    pub fn y(&self) -> f64 {
        self.coords[1]
    }
    pub fn z(&self) -> f64 {
        self.coords[2]
    }
}

impl PartialEq for XYZ {
    fn eq(&self, other: &Self) -> bool {
        self.coords == other.coords
    }
}

impl Eq for XYZ {}

impl Hash for XYZ {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.coords.iter().for_each(|n| {
            n.to_bits().hash(state);
        });
    }
}

impl KdPoint for XYZ {
    type Scalar = f64;
    type Dim = typenum::U3;
    fn at(&self, i: usize) -> f64 {
        self.coords[i]
    }
}

pub fn check_cutoff(a: XYZ, b: XYZ, cutoff: f64) -> bool {
    let d_x = a.x() - b.x();
    let d_y = a.y() - b.y();
    let d_z = a.z() - b.z();
    d_x * d_x + d_y * d_y + d_z * d_z <= cutoff * cutoff
}

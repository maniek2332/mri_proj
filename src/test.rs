extern crate nalgebra as na;
use na::{DMat};
use approxl1_i0;

pub fn test_approx1l() {
    let a: DMat<f32> = DMat::new_random(4,4);
    println!("Hello, world!");
    approxl1_i0::compute(a);
}

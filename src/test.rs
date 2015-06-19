extern crate nalgebra as na;
use na::{DMat};
use approxl1_i0;
use matlab_fun;

pub fn test_approx1l() {
    println!("Starting test_approx1l");
    let a: DMat<f32> = DMat::new_ones(1,3);
    let b = approxl1_i0::compute(a);
	matlab_fun::print_matrix(&b);
}

pub fn test_besseli() {
	println!("Starting test_besseli");
	let a: DMat<f32> = DMat::new_ones(1,3);
	
	let b: DMat<f32> = matlab_fun::besseli(1, &a);
	matlab_fun::print_matrix(&b);
}

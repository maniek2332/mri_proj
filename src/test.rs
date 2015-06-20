extern crate nalgebra as na;
use na::{DMat};
use approxl1_i0;
use matlab_fun;
use correct_rice_gauss;
use em_ml_rice2D;


pub fn test_approx1l() {
    println!("Starting test_approx1l");
    let a: DMat<f32> = DMat::new_ones(1,3);
    let b = approxl1_i0::compute(&a);
	matlab_fun::print_matrix(&b);
}

pub fn test_besseli() {
	println!("Starting test_besseli");
	let a: DMat<f32> = DMat::new_ones(1,3);
	
	let b: DMat<f32> = matlab_fun::besseli(1, &a);
	matlab_fun::print_matrix(&b);
}

pub fn test_correct_rice_gauss() {
	println!("Starting test_correct_rice_gauss");
	
	let a: DMat<f32> = DMat::from_elem(1,3, 2.0);
	
	let b = correct_rice_gauss::compute(&a);
	matlab_fun::print_matrix(&b);
}

pub fn test_em_ml_rice2D() {
	println!("Starting test_em_ml_rice2D");
	let a : DMat<f32> = DMat::from_elem(3,3, 5.0); 
	let (b,c) = em_ml_rice2D::compute(&a, 10, 3);
	
	println!("b mat");
	matlab_fun::print_matrix(&b);
	
	println!("c mat");
	matlab_fun::print_matrix(&c);
	
}




extern crate nalgebra as na;
use na::{DMat};
use std::ops::*;
use matlab_fun;
use matrix_math;

pub fn compute(mat : &DMat<f32>, n : i32, ws : &[f32]){
	let mat64 = matlab_fun::cast_to_f64(&mat);
	
	let pr = matrix_math::prod(&ws);
	let ones: DMat<f32> = DMat::new_ones(ws[0], ws[1]);
	let mask = ones.div(pr);
}
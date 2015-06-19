extern crate nalgebra as na;
use na::{DMat};
use std::ops::*;

pub fn find_less_than(mat: &DMat<f32>, value: f32) -> Vec<(usize, usize)> {
   let rows_size = mat.nrows();
   let cols_size = mat.ncols();
   let mut result: Vec<(usize, usize)> = Vec::new();
 
   for r in 0..rows_size {
	for c in 0..cols_size {
		if mat[(r,c)] <  value {
			result.push((r,c));
		}
	}
   }
  
  return result;	
}

pub fn print_matrix(mat: &DMat<f32>) {
   let rows_size = mat.nrows();
   let cols_size = mat.ncols();
   for r in 0..rows_size {
	for c in 0..cols_size {
		println!("({}, {}) = {}",r,c,mat[(r,c)]);
	}
   }
}

pub fn besseli(nu : f32, z : &DMat<f32>) {

}
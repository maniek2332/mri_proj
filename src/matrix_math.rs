extern crate nalgebra as na;
use na::{DMat};
use std::ops::*;


pub fn div_matrix_by_matrix_each_value(mat1: &DMat<f32>, mat2: &DMat<f32>) -> DMat<f32> {
   assert_eq!(mat1.nrows(), mat2.nrows());
   assert_eq!(mat1.ncols(), mat2.ncols());
   
   let rows_size = mat1.nrows();
   let cols_size = mat1.ncols();
   let mut result: DMat<f32> = DMat::new_zeros(rows_size, cols_size);

   for r in 0..rows_size {
	for c in 0..cols_size {
		result[(r,c)] = mat1[(r,c)] / mat2[(r,c)];
	}
   }
  
  return result;
}

pub fn div_scalar_by_matrix(scalar: f32, mat: &DMat<f32>) -> DMat<f32>{
   let rows_size = mat.nrows();
   let cols_size = mat.ncols();
   let mut result: DMat<f32> = DMat::new_zeros(rows_size, cols_size);

   for r in 0..rows_size {
	for c in 0..cols_size {
		result[(r,c)] = scalar / mat[(r,c)];
	}
   }
  
  return result;
}


pub fn power_matrix_by_scalar(mat: &DMat<f32>, s: i32) -> DMat<f32> {
   let rows_size = mat.nrows();
   let cols_size = mat.ncols();
   let mut result: DMat<f32> = DMat::new_zeros(rows_size, cols_size);
 
   for r in 0..rows_size {
	for c in 0..cols_size {
		result[(r,c)] = pow(mat[(r,c)], s);
	}
   }
  
  return result;
}

pub fn pow(x: f32, s: i32) -> f32{
  let mut result: f32 = 1.0;

  for _ in 0..s {
	result = result * x;
  }
  return result;
}
extern crate nalgebra as na;
use na::{DMat};
use std::f32::consts::PI;
use std::cmp::Ordering;
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

pub fn mul_matrix_by_matrix_each_value(mat1: &DMat<f32>, mat2: &DMat<f32>) -> DMat<f32> {
   assert_eq!(mat1.nrows(), mat2.nrows());
   assert_eq!(mat1.ncols(), mat2.ncols());
   
   let rows_size = mat1.nrows();
   let cols_size = mat1.ncols();
   let mut result: DMat<f32> = DMat::new_zeros(rows_size, cols_size);

   for r in 0..rows_size {
	for c in 0..cols_size {
		result[(r,c)] = mat1[(r,c)] * mat2[(r,c)];
	}
   }
  
  return result;
}

pub fn prod(vec : &[f32]) -> f32 {
	let mut result = 1.0;
	
	for x in vec {
		result = result * x;
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

pub fn matrix_sqrt(mat : &DMat<f32>) -> DMat<f32> {
	let nrows = mat.nrows();
	let ncols = mat.ncols();
	
	let mut result : DMat<f32> = DMat::new_zeros(nrows, ncols);
	
	for r in 0..nrows {
		for c in 0..ncols {
			result[(r,c)] = mat[(r,c)].sqrt();
		}
	}
	
	return result;
}

pub fn mat_map<F>(mat_src: &na::DMat<f32>, mut f: F) -> na::DMat<f32> where F: FnMut(f32) -> f32 {
    let mut mat = na::DMat::new_zeros(mat_src.nrows(), mat_src.ncols());
    for i in 0..mat.nrows() {
        for j in 0..mat.ncols() {
            mat[(i, j)] = f(mat_src[(i, j)]);
        }
    }
    mat
}

pub fn mat_map_mut<F>(mut mat: na::DMat<f32>, mut f: F) -> na::DMat<f32> where F: FnMut(f32) -> f32 {
    for i in 0..mat.nrows() {
        for j in 0..mat.ncols() {
            mat[(i, j)] = f(mat[(i, j)]);
        }
    }
    mat
}

pub fn mat_max<N : Clone + Copy + PartialOrd>(mat: &na::DMat<N>) -> N {
    let mut mat_vec = mat.clone().to_vec();
    mat_vec.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let max_val = *mat_vec.last().unwrap();
    max_val
}

pub fn sum_vec<T>(s: &Vec<T>, init: &T) -> T
    where T : Copy + Add<T, Output=T>
{
    s.iter().fold(*init, |acc, &item| acc + item)    
}

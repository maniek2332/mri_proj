extern crate nalgebra as na;
use na::{DMat};
use std::ops::*;

pub fn find_less_than(mat: &DMat<f64>, value: f64) -> Vec<(usize, usize)> {
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

pub fn find_equal(mat: &DMat<f64>, value: f64) -> Vec<(usize, usize)> {
   let rows_size = mat.nrows();
   let cols_size = mat.ncols();
   let mut result: Vec<(usize, usize)> = Vec::new();
 
   for r in 0..rows_size {
	for c in 0..cols_size {
		if mat[(r,c)] ==  value {
			result.push((r,c));
		}
	}
   }
  
  return result;	
}

pub fn print_matrix(mat: &DMat<f64>) {
   let rows_size = mat.nrows();
   let cols_size = mat.ncols();
   for r in 0..rows_size {
	for c in 0..cols_size {
		println!("({}, {}) = {:.*}",r,c,6,mat[(r,c)]);
	}
   }
}

//to zwróci macierz 1xN gdzie n to ilosc indeksow
pub fn select_cells(mat : &DMat<f64>, indexes: &Vec<(usize, usize)>) -> DMat<f64> {
   let mut result: DMat<f64> = DMat::new_zeros(1, indexes.len());
 
   for i in 0..indexes.len() {
		result[(0, i)] = mat[indexes[i]];
   }
  
  return result;
}


pub fn besseli(nu : usize, mat : &DMat<f64>) -> DMat<f64>{
   let rows_size = mat.nrows();
   let cols_size = mat.ncols();
   let mut result: DMat<f64> = DMat::new_zeros(rows_size, cols_size);

   for r in 0..rows_size {
	for c in 0..cols_size {
		result[(r,c)] = besseli_for_value(nu, mat[(r,c)]);
	}
   }
  
  return result;
}

fn besseli_for_value(nu : usize, x : f64) -> f64 {
	let powers;
	let denominators;
	
	if nu == 0 {
		powers = next_powers(x, 0, 8, 2);
		denominators = [1.0, 4.0, 64.0, 2304.0, 147456.0]; 
	} else {
		powers = next_powers(x, 1, 9, 2);
		denominators = [2.0, 16.0, 384.0, 18432.0, 737280.0]; 
	}
	
	return compute_besseli(&powers, &denominators);
}

pub fn next_powers(x: f64, from: i32, to: i32, step: i32) -> Vec<f64>{
  let mut vec: Vec<f64> = Vec::new();
  let mut result : f64 = 1.0;	
  let mut c = 1;
  let end = to + 1;
  if from == 0 {
	vec.push(1.0);
	c = step;
  }
  
  for i in 1..end {
	result = result * x;
	if i == c {
		vec.push(result);
		c = c + step;
	}	
  }
  
  return vec;
}

fn compute_besseli(powers: &Vec<f64>, denominators: &[f64]) -> f64 {
	let mut result : f64 = 0.0;
	
	for i in  0..powers.len() {
		result = result + powers[i]/denominators[i];
	}
	
	return result;
}

pub fn cast_to_f64(mat : &DMat<f64>) -> DMat<f64> {
	let rows_size = mat.nrows();
	let cols_size = mat.ncols();
	let mut result : DMat<f64> = DMat::new_zeros(rows_size, cols_size);

	for r in 0..rows_size {
		for c in 0..cols_size {
			result[(r,c)] = mat[(r,c)] as f64;
		}
	}
	return result;
}

pub fn biggest_of_values(mat : &DMat<f64>, value: f64) -> DMat<f64>{
	let rows_size = mat.nrows();
	let cols_size = mat.ncols();
	let mut result : DMat<f64> = DMat::new_zeros(rows_size, cols_size);

	for r in 0..rows_size {
		for c in 0..cols_size {
			if mat[(r,c)] < value {
				result[(r,c)] = value;		
			} else {
				result[(r,c)] = mat[(r,c)];
			}
		}
	}
	return result;
}

pub fn abs(mat: &DMat<f64>) -> DMat<f64> {
	let rows_size = mat.nrows();
	let cols_size = mat.ncols();
	let mut result : DMat<f64> = DMat::new_zeros(rows_size, cols_size);

	for r in 0..rows_size {
		for c in 0..cols_size {
			result[(r,c)] = mat[(r,c)].abs();
		}
	}
	return result;
}

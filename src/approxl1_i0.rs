extern crate nalgebra as na;
use na::{DMat, DVec};
use std::ops::*;
use matrix_math;
use matlab_fun;
use na::Iterable;

pub fn compute(z: &DMat<f32>) -> DMat<f32>{
	let cont: DVec<f32> = check(1.5, &z);
	let z_clone = z.clone();

	let z8 = z_clone.mul(8.0);
	let z8_to_2 = matrix_math::power_matrix_by_scalar(&z8, 2);
	let z8_to_3 = matrix_math::power_matrix_by_scalar(&z8, 3);

	let mn = compute_mn(&z8, &z8_to_2, &z8_to_3);
	let md = compute_md(&z8, &z8_to_2, &z8_to_3);


	let mut m = matrix_math::div_matrix_by_matrix_each_value(&mn, &md);
	let mut k = Vec::new();
	if check_cont(&cont, 1.0) {
		k = matlab_fun::find_less_than(&z, 1.5);
		let zk = matlab_fun::select_cells(&z, &k);
		let besseli1 = matlab_fun::besseli(1, &zk);
		let besseli0 = matlab_fun::besseli(0, &zk);
		let besseli_div = matrix_math::div_matrix_by_matrix_each_value(&besseli1, &besseli0);

		for i in 0..k.len() {
			m[k[i]] = besseli_div[(0,i)];
		}
	}
  
	let zeros: Vec<(usize, usize)> = matlab_fun::find_equal(&z, 0.0);
	for i in 0..zeros.len() {
		m[zeros[i]] = 0.0;
	}
	return m;
}

fn check(treshold: f32, z: &DMat<f32>) -> DVec<f32> {
	let rows_size = z.nrows();
	let cols_size = z.ncols();
	let mut  result: DVec<f32> = DVec::new_zeros(cols_size);
	
	for c in 0..cols_size {
		for r in 0..rows_size {
			if z[(r,c)] < treshold {
				result[c] = result[c]  + 1.0;
			}
		}
   }
  return result;
}

fn check_cont(cont : &DVec<f32>, value: f32) -> bool {
	for x in cont.iter() {
		if x <= &value {
			return false;
		}
	}
	return true;
}

fn compute_mn(z8: &DMat<f32>, z8_to_2: &DMat<f32>, z8_to_3: &DMat<f32>) -> DMat<f32> {
   let z8_clone = z8.clone();
  //Mn=1-a-b-c
  
  //a = 3./z8
  let a = matrix_math::div_scalar_by_matrix(3.0, &z8_clone);
  //b = 15./2./(z8).^2
  let b = matrix_math::div_scalar_by_matrix(7.5, &z8_to_2);


  //c =52.5./(z8).^3
  let c = matrix_math::div_scalar_by_matrix(52.5, &z8_to_3);
  
  let ones = DMat::new_ones(z8_clone.nrows(), z8_clone.ncols());
  //return z8.clone();
  return ones.sub(a).sub(b).sub(c);
}

fn compute_md(z8: &DMat<f32>, z8_to_2: &DMat<f32>, z8_to_3: &DMat<f32>) -> DMat<f32> {
	let z8_clone = z8.clone();
	//Md =1+a+b+c;
	
	//a = 1/z8
	let a = matrix_math::div_scalar_by_matrix(1.0, &z8_clone);
	
	//b = 4.5/(z8)^2
	let b = matrix_math::div_scalar_by_matrix(4.5, &z8_to_2);
	
	//c = 37.5/(z8)^3
	let c = matrix_math::div_scalar_by_matrix(37.5, &z8_to_3);
	
	let ones = DMat::new_ones(z8_clone.nrows(), z8_clone.ncols());

	return ones.add(a).add(b).add(c);
}


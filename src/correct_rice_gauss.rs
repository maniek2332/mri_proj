extern crate nalgebra as na;
use na::{DMat};
use std::ops::*;
use matlab_fun;
use matrix_math;

pub fn compute(snr : &DMat<f32>) -> DMat<f32> {
	let coefs : [f32;9] = [-0.289549906258443, -0.0388922575606330, 0.409867108141953, -0.355237628488567,	0.149328280945610, -0.0357861117942093,	0.00497952893859122, -0.000374756374477592, 0.0000118020229140092];
	
	//Fc=Coefs(1)+Coefs(2).*a1+Coefs(3).*a1.^2+Coefs(4).*a1.^3+Coefs(5).*a1.^4+Coefs(6).*a1.^5+Coefs(7).*a1.^6+Coefs(8).*a1.^7+Coefs(9).*a1.^8;
	
	let rows_size = snr.nrows();
	let cols_size = snr.ncols();
	let mut fc : DMat<f32> = DMat::new_zeros(rows_size, cols_size);

   for r in 0..rows_size {
	for c in 0..cols_size {
	    let powers = matlab_fun::next_powers(snr[(r,c)], 0, 8, 1); 
		fc[(r,c)] = coefs[0] + coefs[1] * powers[1] + coefs[2] * powers[2] + coefs[3] * powers[3] + coefs[4] * powers[4]
		+ coefs[5] * powers[5] + coefs[6] * powers[6] + coefs[7] * powers[7] + coefs[8] * powers[8];
	}
   }
   
   //matlab_fun::print_matrix(&fc);
	let is_less_than = is_less_then(snr, 7.0);
	
	let fc_final = matrix_math::mul_matrix_by_matrix_each_value(&fc, &is_less_than);
	
	return fc_final;
}

fn is_less_then(mat: &DMat<f32>, value: f32) -> DMat<f32> {
	let rows = mat.nrows();
	let cols = mat.ncols();

	let mut result : DMat<f32> = DMat::new_zeros(rows, cols);

	for r in 0..rows {
		for c in 0..cols {
			if mat[(r,c)] <  value {
				result[(r,c)] = 1.0;
			}
		}
	}
   
	return result;
}
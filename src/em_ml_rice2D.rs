extern crate nalgebra as na;
use na::{DMat};
use std::ops::*;
use matlab_fun;
use matrix_math;
use filter2b;
use approxl1_i0;

pub fn compute(mat : &DMat<f64>, n : i32, ws : usize) -> (DMat<f64>, DMat<f64>){
	//jezeli błąd będzie za duży spróbować użyć f64
	//szkoda że to wymaga przerobienia połowy programu...
	//let mat64 = matlab_fun::cast_to_f64(&mat);
	
	let tt: f64 = (ws * ws) as f64;
	let mask_value: f64 = 1.0/tt;
	let mask_size: usize = ws;
	
	let mat_to_2 = matrix_math::power_matrix_by_scalar(&mat, 2);
	let mat_to_4 = matrix_math::power_matrix_by_scalar(&mat, 4);
	
	let mat_to_2_filter = filter2b::filter_img(&mat_to_2, mask_size, mask_value);
	let mat_to_4_filter = filter2b::filter_img(&mat_to_4, mask_size, mask_value);
	
	let mut a_k = compute_a_k(&mat_to_2_filter, &mat_to_4_filter);
	let mut sigma_k2 = compute_sigma_k(&a_k, &mat_to_2_filter);
	
	for _ in 0..n {
		a_k = compute_next_a_k(&a_k, &mat, &sigma_k2, mask_size, mask_value);
		sigma_k2 = compute_next_signal_k(&mat, &a_k, mask_size, mask_value);
	}

	let sigma_n = matrix_math::matrix_sqrt(&sigma_k2);
	
	return (a_k, sigma_n);	
}

fn compute_a_k(mat_to_2_filter : &DMat<f64>, mat_to_4_filter: &DMat<f64>) -> DMat<f64> {

	let mat_to_2_filter_to_2 = matrix_math::power_matrix_by_scalar(&mat_to_2_filter, 2);
	let a = mat_to_2_filter_to_2.mul(2.0);
	
	let b = a.sub(mat_to_4_filter.clone());
	
	let max_b = matlab_fun::biggest_of_values(&b, 0.0);
	
	let max_b_sqrt = matrix_math::matrix_sqrt(&max_b);
	
	return matrix_math::matrix_sqrt(&max_b_sqrt);
}

fn compute_sigma_k(a_k: &DMat<f64>, mat_to_2_filter: &DMat<f64>) -> DMat<f64> {
	let a_k_to_2 = matrix_math::power_matrix_by_scalar(&a_k, 2);
	let c = mat_to_2_filter.clone().sub(a_k_to_2);
	
	let d = matlab_fun::biggest_of_values(&c, 0.01);
	
	return d.mul(0.5);	
}

fn compute_next_a_k(a_k: &DMat<f64>, mat : &DMat<f64>, sigma_k: &DMat<f64>, mask_size: usize, mask_value: f64) -> DMat<f64> {	
	let a_k_times_mat = matrix_math::mul_matrix_by_matrix_each_value(&a_k, &mat);
	let ak_div_sigma = matrix_math::div_matrix_by_matrix_each_value(&a_k_times_mat, &sigma_k);
	let approx = approxl1_i0::compute(&ak_div_sigma);	
	let approx_times_mat = matrix_math::mul_matrix_by_matrix_each_value(&approx, &mat); 
	let filtered = filter2b::filter_img(&approx_times_mat, mask_size, mask_value);	
	let new_ak =  matlab_fun::biggest_of_values(&filtered, 0.0);
	
	return new_ak;
}

fn compute_next_signal_k(mat : &DMat<f64>, a_k: &DMat<f64>, mask_size: usize, mask_value: f64) -> DMat<f64> {
	let abs_mat = matlab_fun::abs(&mat);
	let abs_to_2 = matrix_math::power_matrix_by_scalar(&abs_mat, 2);
	let filtered = filter2b::filter_img(&abs_to_2, mask_size, mask_value);
	let filtered_times_half = filtered.mul(0.5);
	let a_k_to_2 = matrix_math::power_matrix_by_scalar(&a_k, 2);
	let a_k_to_2_half = a_k_to_2.div(2.0);
	let temp = filtered_times_half.sub(a_k_to_2_half);
	
	return matlab_fun::biggest_of_values(&temp, 0.01);
}

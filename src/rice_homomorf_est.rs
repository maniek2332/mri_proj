extern crate nalgebra as na;
use na::{DMat};
use std::ops::*;
use matlab_fun;
use matrix_math;
use filter2b;
use em_ml_rice2D;
use lpf;
use correct_rice_gauss;

pub fn compute_for_uknown_snr(mat : &DMat<f64>, lpf: f64, modo: i32) -> (DMat<f64>, DMat<f64>){
	let (m2, sigma_n) = em_ml_rice2D::compute(&mat,10, 3);
	let snr = matrix_math::div_matrix_by_matrix_each_value(&m2, &sigma_n);
	
	return compute(&mat, &snr, lpf, modo);
}

pub fn compute(mat : &DMat<f64>, snr: &DMat<f64>, _lpf: f64, modo: i32) -> (DMat<f64>, DMat<f64>){
	let (m2, sigma_n) = em_ml_rice2D::compute(&mat,10, 3);
	
	let m1= filter2b::filter_img(&mat, 5, 0.04);
    println!("M1:\n{:?}\n", m1);
	
	let sub = mat.clone().sub(m1.clone());
	let rn = matlab_fun::abs(&sub);
	let vioz = value_or_zeros(&rn, 0.001);
	let rn_vioz = rn.add(vioz);
	
	let l_rn = matrix_math::mat_map(&rn_vioz, |x| x.ln());
	let lpf2 = lpf::lpf2(l_rn, _lpf);
	let mapa2 = matrix_math::mat_map(&lpf2, |x| x.exp());
	
	//2./sqrt(2).*exp(-psi(1)./2) = 1.8874
	let const_val: f64 = 1.8874;
	let mapa_g = matrix_math::mat_map(&mapa2, |x| x * const_val);
	
	let mut rn_abs = DMat::new_ones(1,1);
	
	if modo == 1 {
		let temp = mat.clone().sub(m1.clone());
		rn_abs = matlab_fun::abs(&temp);
	} else if modo==2{
		let temp = mat.clone().sub(m2);
		rn_abs = matlab_fun::abs(&temp);
	} else {
		rn_abs = matlab_fun::abs(&mat);
	}
	
	let vioz_2 = value_or_zeros(&rn_abs, 0.001);
	let rn_abs_vioz_2 = rn_abs.add(vioz_2);
	let l_rn_abs = matrix_math::mat_map(&rn_abs_vioz_2, |x| x.ln());
	let lpf2_abs = lpf::lpf2(l_rn_abs, _lpf);
	
	let fc1= correct_rice_gauss::compute(&snr);
	let lpf1 = lpf2_abs.sub(fc1);
	let lpf1_2 = lpf::lpf2(lpf1,_lpf + 2.0);		
		
	let mapa1 = matrix_math::mat_map(&lpf1_2, |x| x.exp());

	let mapa_r = matrix_math::mat_map(&mapa1, |x| x * const_val);
	
	return (mapa_r, mapa_g);
}

fn value_or_zeros(mat : &DMat<f64>, value: f64) -> DMat<f64> {
   let rows_size = mat.nrows();
   let cols_size = mat.ncols();
   let mut result: DMat<f64> = DMat::new_zeros(rows_size, cols_size);
 
   for r in 0..rows_size {
	for c in 0..cols_size {
		if(mat[(r,c)] ==  0.0 ) {
			result[(r,c)] = value;
		}
	}
   }
  
  return result;	
}

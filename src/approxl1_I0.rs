extern crate nalgebra as na;
use na::{DMat};
use std::ops::*;
use matrix_math;

pub fn compute(z: DMat<f32>) {
  //let cont = check(1.5, z.as_vec());
  
  let z8 = z.mul(8.0);
  let z8_to_2 = matrix_math::power_matrix_by_scalar(&z8, 2);
  let z8_to_3 = matrix_math::power_matrix_by_scalar(&z8, 3);
  
  let mn = compute_mn(&z8, &z8_to_2, &z8_to_3);
  let md = compute_md(&z8, &z8_to_2, &z8_to_3);
  
  let m = matrix_math::div_matrix_by_matrix_each_value(&mn, &md);
}

fn check(treshold: f32, z: &[f32]) -> bool {
  for x in z {
    if x < &treshold {
      return true;
    }
  }

  return false;
}

fn compute_mn(z8: &DMat<f32>, z8_to_2: &DMat<f32>, z8_to_3: &DMat<f32>) -> DMat<f32> {
   let z8_clone = z8.clone();
  //Mn=1-a-b-c
  
  //a = 3./z8
  println!("z {} ", z8.index((0,0)));
  let a = matrix_math::div_scalar_by_matrix(3.0, &z8_clone);
  println!("a {} ", a.index((0,0)));
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

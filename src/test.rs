extern crate nalgebra as na;
use na::{DMat};
use na::Transpose;
use approxl1_i0;
use matlab_fun;
use correct_rice_gauss;
use em_ml_rice2D;
use rice_homomorf_est;
use lpf::{lpf2, dct_inplace, idct_inplace, dct2, idct2};


pub fn test_approx1l() {
    println!("Starting test_approx1l");
    let a: DMat<f32> = DMat::new_ones(1,3);
    let b = approxl1_i0::compute(&a);
	matlab_fun::print_matrix(&b);
}

pub fn test_besseli() {
	println!("Starting test_besseli");
	let a: DMat<f32> = DMat::new_ones(1,3);
	
	let b: DMat<f32> = matlab_fun::besseli(1, &a);
	matlab_fun::print_matrix(&b);
}

pub fn test_correct_rice_gauss() {
	println!("Starting test_correct_rice_gauss");
	
	let a: DMat<f32> = DMat::from_elem(1,3, 2.0);
	
	let b = correct_rice_gauss::compute(&a);
	matlab_fun::print_matrix(&b);
}

pub fn test_em_ml_rice2D() {
	println!("Starting test_em_ml_rice2D");
	let a : DMat<f32> = DMat::from_elem(3,3, 5.0); 
	let (b,c) = em_ml_rice2D::compute(&a, 10, 3);
	
	println!("b mat");
	matlab_fun::print_matrix(&b);
	
	println!("c mat");
	matlab_fun::print_matrix(&c);	
}

pub fn test_rice_homomorf_est() {
	println!("Starting test_rice_homomorf_est");
	let a : DMat<f32> = DMat::from_elem(4,4, 5.0); 
	
	let (b,c) = rice_homomorf_est::compute_for_uknown_snr(&a, 3.4, 2);
	
	println!("b mat");
	matlab_fun::print_matrix(&b);
	
	println!("c mat");
	matlab_fun::print_matrix(&c);
}

pub fn test_lpf() {
    let mut test_img = na::DMat::from_elem(6, 6, 25.);
    let test_img_lpf = lpf2(test_img, 4.8);
    println!("TEST LPF (1):\n{:?}\n", test_img_lpf);
    let mut test_img2 = na::DMat::from_col_vec(2, 2,
                                               &[25., 12., 8.33, 25.]);
    test_img2.transpose_mut();
    let test_img2_lpf = lpf2(test_img2, 4.8);
    println!("TEST LPF (2):\n{:?}\n", test_img2_lpf);

}

pub fn test_dct() {
    let mut test_vec = na::DVec::from_slice(6, &[25., 25., 25., 25., 25., 25.]);
    println!("DCT VEC  (1): {:?}", test_vec);
    dct_inplace(&mut test_vec);
    println!("DCT TEST (1): {:?}", test_vec);
    let mut test_vec2 = na::DVec::from_slice(6, &[25., 25., 8.33, 25., 12.5, 0.]);
    println!("DCT VEC  (2): {:?}", test_vec2);
    dct_inplace(&mut test_vec2);
    println!("DCT TEST (2): {:?}", test_vec2);
}

pub fn test_idct() {
    let mut test_vec = na::DVec::from_slice(6, &[25., 25., 25., 25., 25., 25.]);
    println!("IDCT VEC  (1): {:?}", test_vec);
    idct_inplace(&mut test_vec);
    println!("IDCT TEST (1): {:?}", test_vec);
    let mut test_vec2 = na::DVec::from_slice(6, &[25., 25., 8.33, 25., 12.5, 0.]);
    println!("IDCT VEC  (2): {:?}", test_vec2);
    idct_inplace(&mut test_vec2);
    println!("IDCT TEST (2): {:?}", test_vec2);
}

pub fn test_dct2() {
    let mut test_mat = na::DMat::from_elem(3, 3, 25.);
    println!("DCT MAT  (1):\n{:?}\n", test_mat);
    println!("DCT TEST (1):\n{:?}\n", dct2(&test_mat));
    let mut test_mat2 = na::DMat::from_row_vec(3, 3, &[3., 0., 3., 2., 2., 1., 0.5, 0.7, 5.]);
    println!("DCT MAT  (2):\n{:?}\n", test_mat2);
    println!("DCT TEST (2):\n{:?}\n", dct2(&test_mat2));
}
pub fn test_idct2() {
    let mut test_mat = na::DMat::from_elem(3, 3, 25.);
    println!("IDCT MAT  (1):\n{:?}\n", test_mat);
    println!("IDCT TEST (1):\n{:?}\n", idct2(&test_mat));
    let mut test_mat2 = na::DMat::from_row_vec(3, 3, &[3., 0., 3., 2., 2., 1., 0.5, 0.7, 5.]);
    println!("IDCT MAT  (2):\n{:?}\n", test_mat2);
    println!("IDCT TEST (2):\n{:?}\n", idct2(&test_mat2));
}

pub fn test_dct2_and_idct2() {
    let mut test_mat2 = na::DMat::from_row_vec(3, 3, &[3., 0., 3., 2., 2., 1., 0.5, 0.7, 5.]);
    println!("DCT IDCT MAT  (1):\n{:?}\n", test_mat2);
    println!("DCT IDCT TEST (1):\n{:?}\n", idct2(&dct2(&test_mat2)));
}

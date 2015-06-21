extern crate nalgebra as na;
use std::f64::consts::PI;
use matrix_math::{sum_vec, mat_map_mut, mat_max};
use na::{RowSlice, ColSlice};
use na::Transpose;

pub fn dct_inplace(y : &mut na::DVec<f64>) {
    let x = y.clone();
    for k in 0..x.len() {
        let mut acc = 0f64;
        for n in 0..x.len() {
            acc += x[n] * (PI * (k as f64) * (2f64 * (n as f64) + 1f64) / (2f64 * (x.len() as f64))).cos();
        }
        y[k] = acc * 2f64;
        if k == 0 {
            y[k] = y[k] * (1f64/(4f64 * x.len() as f64)).sqrt();
        } else {
            y[k] = y[k] * (1f64/(2f64 * x.len() as f64)).sqrt();
        }
    }
}

pub fn idct_inplace(y : &mut na::DVec<f64>) {
    let x = y.clone();

    for k in 0..x.len() {
        let mut acc = 0f64;
        for n in 1..x.len() {
            acc += x[n] * (PI * (k as f64 + 0.5f64) * n as f64 / x.len() as f64).cos();
        }
        y[k] = x[0] / (x.len() as f64).sqrt() + acc * (2f64 / x.len() as f64).sqrt();
    }
}

pub fn dct2(mat_src : &na::DMat<f64>) -> na::DMat<f64>  {
    let mut mat = mat_src.clone();
    //println!("MID DCT2 MATRIX (1):\n{:?}\n", mat);
    for i in 0..mat.nrows() {
        let mut slice = mat.row_slice(i, 0, mat.ncols());
        dct_inplace(&mut slice);
        for j in 0..mat.ncols() {
            mat[(i, j)] = slice[j];
        }
    }
    for i in 0..mat.ncols() {
        let mut slice = mat.col_slice(i, 0, mat.nrows());
        //println!("DCT SLICE COL (1): {:?}", slice);
        dct_inplace(&mut slice);
        //println!("DCT SLICE COL (2): {:?}", slice);
        for j in 0..mat.nrows() {
            mat[(j, i)] = slice[j];
        }
        //println!("DCT MAT SLICE UPDATE ():\n{:?}\n", mat);
    }
    mat
}

pub fn idct2(mat_src : &na::DMat<f64>) -> na::DMat<f64>  {
    let mut mat = mat_src.clone();
    for i in 0..mat.nrows() {
        let mut slice = mat.row_slice(i, 0, mat.ncols());
        idct_inplace(&mut slice);
        for j in 0..mat.ncols() {
            mat[(i, j)] = slice[j];
        }
    }
    for i in 0..mat.ncols() {
        let mut slice = mat.col_slice(i, 0, mat.nrows());
        idct_inplace(&mut slice);
        for j in 0..mat.nrows() {
            mat[(j, i)] = slice[j];
        }
    }
    mat
}

fn gaussian_filter(rows : usize, cols : usize, sigma : f64) -> na::DMat<f64> {
    assert_eq!(rows % 2, 0);
    assert_eq!(cols % 2, 0);
    let mut mat = na::DMat::new_zeros(rows, cols);
    for y in 0..rows {
        let yv = y as f64 + 0.5f64;
        for x in 0..cols {
            let xv = x as f64 + 0.5f64;
            mat[(y, x)] = (
                (-(xv * xv + yv * yv) / (2f64 * sigma * sigma)).exp()
                );
        }
    }
    let mat_sum = sum_vec(&mat.clone().to_vec(), &0.);
    mat = mat_map_mut(mat, |x| x / (mat_sum * 4.0f64));
    mat
}

pub fn lpf2(mat : na::DMat<f64>, sigma : f64) -> na::DMat<f64> {
    let mut gauss_mat = gaussian_filter(mat.nrows(), mat.ncols(), sigma * 2f64);
    let gauss_max = mat_max(&gauss_mat);
    gauss_mat = mat_map_mut(gauss_mat, |x| x / gauss_max);
    //println!("TEST IMG:\n{:?}\n", mat);
    //println!("GAUSS MAT:\n{:?}\n", gauss_mat);
    let mut mat_dct = dct2(&mat);
    //println!("DCT MAT:\n{:?}\n", mat_dct);
    for i in 0..mat.nrows() {
        for j in 0..mat.ncols() {
            mat_dct[(i, j)] = mat_dct[(i, j)] * gauss_mat[(i, j)];
        }
    }
    //println!("DCT MAT (with Gauss):\n{:?}\n", mat_dct);
    let mat_lpf = idct2(&mat_dct);
    mat_lpf
}


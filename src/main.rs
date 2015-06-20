extern crate nalgebra as na;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Write;
use std::cmp::{min, max};
use std::cmp::Ordering;
use std::f32::consts::PI;
use na::{RowSlice, ColSlice};

mod approxl1_i0;
mod matrix_math;
mod matlab_fun;
mod test;
mod correct_rice_gauss;
mod em_ml_rice2D;
mod filter2b;
mod config;
mod rice_homomorf_est;

fn load_image(path : String) -> na::DMat<f32> {
    let mut mat = na::DMat::new_zeros(256, 256);

    let file = match File::open(path) {
        Ok(file) => file,
        Err(..)  => panic!("panic"),
    };
    let br = BufReader::new(&file);

    let mut row = 0;
    for line in br.lines().map(|line| line.unwrap()) {
        let mut col = 0;
        for num_str in line.split(",") {
            let num = num_str.parse::<f32>().unwrap();
            mat[(row, col)] = num;
            col += 1;
        }
        row += 1;
    }
    mat
}

fn save_image(img : &na::DMat<f32>, path : String) {
    let mut file = File::create(path).unwrap();
    for i in 0..img.nrows() {
        for j in 0..img.ncols() - 1 {
            write!(file, "{},", img[(i, j)]);
        }
        write!(file, "{}\n", img[(i, img.ncols() - 1)]);
    }
}

fn dct_inplace(y : &mut na::DVec<f32>) {
    let x = y.clone();
    for k in 0..x.len() {
        let mut acc = 0f32;
        for n in 0..x.len() {
            acc += x[n] * (PI * (k as f32) * (2f32 * (n as f32) + 1f32) / (2f32 * (x.len() as f32))).cos();
        }
        y[k] = acc * 2f32;
        if k == 0 {
            y[k] = y[k] * (1f32/(4f32 * x.len() as f32)).sqrt();
        } else {
            y[k] = y[k] * (1f32/(2f32 * x.len() as f32)).sqrt();
        }
    }
}

fn idct_inplace(y : &mut na::DVec<f32>) {
    let x = y.clone();

    for k in 0..x.len() {
        let mut acc = 0f32;
        for n in 1..x.len() {
            acc += x[n] * (PI * (k as f32 + 0.5f32) * n as f32 * x.len() as f32).cos();
        }
        y[k] = x[0] / (x.len() as f32).sqrt() + acc * (2f32 / x.len() as f32).sqrt();
    }
}

fn dct2(mat_src : &na::DMat<f32>) -> na::DMat<f32>  {
    let mut mat = mat_src.clone();
    for i in 0..mat.nrows() {
        let mut slice = mat.row_slice(i, 0, mat.ncols());
        dct_inplace(&mut slice);
        for j in 0..mat.ncols() {
            mat[(i, j)] = slice[j];
        }
    }
    for i in 0..mat.ncols() {
        let mut slice = mat.col_slice(i, 0, mat.nrows());
        dct_inplace(&mut slice);
        for j in 0..mat.nrows() {
            mat[(j, i)] = slice[j];
        }
    }
    mat
}

fn idct2(mat_src : &na::DMat<f32>) -> na::DMat<f32>  {
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

fn gaussian_filter(rows : usize, cols : usize, sigma : f32) -> na::DMat<f32> {
    let mut mat = na::DMat::new_zeros(rows, cols);
    for y in 0..rows {
        for x in 0..cols {
            mat[(y, x)] = (
                (-(((x * x) + (y * y)) as f32) / (2. * sigma * sigma)).exp()
                );
        }
    }
    mat
}

fn mat_map<F>(mat_src: &na::DMat<f32>, mut f: F) -> na::DMat<f32> where F: FnMut(f32) -> f32 {
    let mut mat = na::DMat::new_zeros(mat_src.nrows(), mat_src.ncols());
    for i in 0..mat.nrows() {
        for j in 0..mat.ncols() {
            mat[(i, j)] = f(mat_src[(i, j)]);
        }
    }
    mat
}

fn mat_map_mut<F>(mut mat: na::DMat<f32>, mut f: F) -> na::DMat<f32> where F: FnMut(f32) -> f32 {
    for i in 0..mat.nrows() {
        for j in 0..mat.ncols() {
            mat[(i, j)] = f(mat[(i, j)]);
        }
    }
    mat
}

fn mat_max<N : Clone + Copy + PartialOrd>(mat: &na::DMat<N>) -> N {
    let mut mat_vec = mat.clone().to_vec();
    mat_vec.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let max_val = *mat_vec.last().unwrap();
    max_val
}

fn lpf2(mat : na::DMat<f32>, sigma : f32) -> na::DMat<f32> {
    let mut gauss_mat = gaussian_filter(mat.nrows(), mat.ncols(), sigma * 2.);
    let gauss_max = mat_max(&gauss_mat);
    //gauss_mat = mat_map_mut(gauss_mat, |x| x / gauss_max);
    println!("GAUSS MAT:\n{:?}", gauss_mat);
    let mut mat_dct = dct2(&mat);
    mat_dct = mat_dct * gauss_mat;
    let mat_lpf = idct2(&mat_dct);
    mat_lpf
}

fn main() {
	test::test_em_ml_rice2D();
    let img = load_image("test.csv".to_string());
    let img_f = filter2b::filter_img(&img, 5, 1. / 25.);
    println!("IMG1: {:?}", img[(0, 0)]);
    println!("IMG2: {:?}", img_f[(0, 0)]);
    save_image(&img_f, "output.csv".to_string());
    let mut test_img = na::DMat::from_elem(5, 5, 25.);
    //test_img[(5, 5)] = 5.;
    //let img_dct = dct2(&test_img);
    //let img_idct = idct2(&test_img);
    //println!("IMG:\n{:?}\n", test_img);
    //println!("IMG DCT:\n{:?}\n", img_dct);
    //println!("IMG IDCT:\n{:?}\n", img_idct);
    //let mat_gauss = gaussian_filter(7, 7, 0.84);
    //println!("GAUSS:\n{:?}\n", mat_gauss);
    let test_img_lpf = lpf2(test_img, 4.8);
    println!("TEST LPF:\n{:?}\n", test_img_lpf);
    //let mut test_vec = na::DVec::from_elem(5, 25.);
    ////test_vec[1] = 23f32;
    //println!("TEST VEC: {:?}", test_vec);
    //dct_inplace(&mut test_vec);
    //println!("VEC DCT: {:?}", test_vec);
    //idct_inplace(&mut test_vec);
    //println!("VEC IDCT: {:?}", test_vec);
}

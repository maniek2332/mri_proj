extern crate nalgebra as na;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Write;
use std::cmp::{min, max};
use std::f64::consts::PI;
use std::ops::Add;
use std::vec::Vec;
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
mod lpf;

use lpf::lpf2;

fn load_image(path : String) -> na::DMat<f64> {
    let mut mat = na::DMat::new_zeros(256, 256);

    let file = match File::open(path) {
        Ok(file) => file,
        Err(..)  => panic!("Cannot open file"),
    };
    let br = BufReader::new(&file);

    let mut row = 0;
    for line in br.lines().map(|line| line.unwrap()) {
        let mut col = 0;
        for num_str in line.split(",") {
            let num = num_str.parse::<f64>().unwrap();
            mat[(row, col)] = num;
            col += 1;
        }
        row += 1;
    }
    mat
}

fn save_image(img : &na::DMat<f64>, path : String) {
    let mut file = File::create(path).unwrap();
    for i in 0..img.nrows() {
        for j in 0..img.ncols() - 1 {
            write!(file, "{},", img[(i, j)]);
        }
        write!(file, "{}\n", img[(i, img.ncols() - 1)]);
    }
}

fn run() {
    let ref config = config::CONFIG;
    let mr_noisy = load_image(config.input_filename.to_string());
    let mr_snr = load_image(config.input_snr.to_string());
    println!("Starting algorithm...");
    let (mr_rice_map, mr_gauss_map) = rice_homomorf_est::compute(&mr_noisy, &mr_snr, 3.4, 2);
    //let (mr_rice_map, mr_gauss_map) = rice_homomorf_est::compute_for_uknown_snr(&mr_noisy, 3.4, 2);
    //let (mr_rice_map, mr_gauss_map) = rice_homomorf_est::compute(&mr_noisy, &mr_snr, 3.4, 1);
    //let (mr_rice_map, mr_gauss_map) = rice_homomorf_est::compute_for_uknown_snr(&mr_noisy, 3.4, 1);
    println!("Completed, saving results");
    save_image(&mr_rice_map, config.output_filename_Rician.to_string());
    save_image(&mr_gauss_map, config.output_filename_Gaussian.to_string());
}

fn main() {
    run();

    //test::test_dct();
    //test::test_idct();
    //test::test_dct2();
    //test::test_idct2();
    //test::test_dct2_and_idct2();
    //test::test_lpf();
    //test::test_rice_homomorf_est();
    //let img = load_image("test.csv".to_string());
    //let img_f = filter2b::filter_img(&img, 5, 1. / 25.);
    //println!("IMG1: {:?}", img[(0, 0)]);
    //println!("IMG2: {:?}", img_f[(0, 0)]);
    //save_image(&img_f, "output.csv".to_string());
    //let mut test_img = na::DMat::from_elem(6, 6, 25.);
    ////test_img[(5, 5)] = 5.;
    ////let img_dct = dct2(&test_img);
    ////let img_idct = idct2(&test_img);
    ////println!("IMG:\n{:?}\n", test_img);
    ////println!("IMG DCT:\n{:?}\n", img_dct);
    ////println!("IMG IDCT:\n{:?}\n", img_idct);
    ////let mat_gauss = gaussian_filter(7, 7, 0.84);
    ////println!("GAUSS:\n{:?}\n", mat_gauss);
    //let test_img_lpf = lpf2(test_img, 4.8);
    //println!("TEST LPF:\n{:?}\n", test_img_lpf);
    //let mut test_vec = na::DVec::from_elem(5, 25.);
    ////test_vec[1] = 23f64;
    //println!("TEST VEC: {:?}", test_vec);
    //dct_inplace(&mut test_vec);
    //println!("VEC DCT: {:?}", test_vec);
    //idct_inplace(&mut test_vec);
    //println!("VEC IDCT: {:?}", test_vec);
}

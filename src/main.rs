extern crate nalgebra as na;

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
mod utils;

use utils::{load_image, save_image};

fn run() {
    let ref config = config::load_config();
    let mr_noisy = load_image(config.input_filename.to_string());
    let mr_snr = load_image(config.input_snr.to_string());
    println!("Starting algorithm...");
    if config.use_snr == 1 && config.ex_filter_type == 2 {
        let (mr_rice_map, mr_gauss_map) = rice_homomorf_est::compute(&mr_noisy, &mr_snr, config.lpf_f, 2, &config);
        println!("Completed, saving results");
        save_image(&mr_rice_map, config.output_filename_Rician.to_string());
        save_image(&mr_gauss_map, config.output_filename_Gaussian.to_string());
    } else if config.use_snr == 0 && config.ex_filter_type == 2 {
        let (mr_rice_map, mr_gauss_map) = rice_homomorf_est::compute_for_uknown_snr(&mr_noisy, config.lpf_f, 2, &config);
        println!("Completed, saving results");
        save_image(&mr_rice_map, config.output_filename_Rician.to_string());
        save_image(&mr_gauss_map, config.output_filename_Gaussian.to_string());
    } else if config.use_snr == 1 && config.ex_filter_type == 1 {
        let (mr_rice_map, mr_gauss_map) = rice_homomorf_est::compute(&mr_noisy, &mr_snr, config.lpf_f, 1, &config);
        println!("Completed, saving results");
        save_image(&mr_rice_map, config.output_filename_Rician.to_string());
        save_image(&mr_gauss_map, config.output_filename_Gaussian.to_string());
    } else if config.use_snr == 0 && config.ex_filter_type == 1 {
        let (mr_rice_map, mr_gauss_map) = rice_homomorf_est::compute_for_uknown_snr(&mr_noisy, config.lpf_f, 1, &config);
        println!("Completed, saving results");
        save_image(&mr_rice_map, config.output_filename_Rician.to_string());
        save_image(&mr_gauss_map, config.output_filename_Gaussian.to_string());
    } else {
        panic!("Config error");
    }
}

fn main() {
    run();
    //test::test_filter();
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

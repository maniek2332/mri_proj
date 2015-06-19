extern crate nalgebra as na;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::cmp::{min, max};

mod approxl1_i0;
mod matrix_math;
mod matlab_fun;
mod test;

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

fn filter(image_src : &na::DMat<f32>, mask_size : usize, mask_val : f32) -> na::DMat<f32> {
    let mut image = image_src.clone();
    assert_eq!(image.nrows(), image.ncols());
    let range = mask_size as i32 / 2;
    let mask = na::DMat::from_elem(mask_size, mask_size, mask_val);
    let size = image.nrows() as i32;

    let pos = |p| min(size - 1, max(0, p)) as usize;
    let mask_pos = |p| (p + range) as usize;

    for x in 0..size {
        for y in 0..size {
            let mut acc = 0f32;
            for xd in -range..range+1 {
                for yd in -range..range+1 {
                    acc += image[(pos(x + xd), pos(y + yd))] as f32 *
                        mask[(mask_pos(xd), mask_pos(yd))] as f32;
                }
            }
            image[(x as usize, y as usize)] = acc;
        }
    }
    image
}

fn main() {
	test::test_approx1l();
    let img = load_image("test.csv".to_string());
    let img_f = filter(&img, 5, 1. / 25.);
    println!("IMG1: {:?}", img[(0, 0)]);
    println!("IMG2: {:?}", img_f[(0, 0)]);
    //println!("Mat: {:?}", mat);
}

extern crate nalgebra as na;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Write;
use std::cmp::{min, max};
use std::f32::consts::PI;

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

fn filter_img(image_src : &na::DMat<f32>, mask_size : usize, mask_val : f32) -> na::DMat<f32> {
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
                    acc += image_src[(pos(x + xd), pos(y + yd))] as f32 *
                        mask[(mask_pos(xd), mask_pos(yd))] as f32;
                }
            }
            image[(x as usize, y as usize)] = acc;
        }
    }
    image
}

fn dct2(image_src : &na::DMat<f32>) -> na::DMat<f32> {
    let mut image = na::DMat::new_zeros(image_src.nrows(), image_src.ncols());
    for u in 0..image.nrows() {
        for v in 0..image.ncols() {
            let mut acc = 0f32;
            for i in 0..image.nrows() {
                for j in 0..image.ncols() {
                    acc += (
                        image_src[(i, j)] *
                        (PI / (image.nrows() as f32) * ((i as f32) + 1. / 2.) * (u as f32)).cos() *
                        (PI / (image.ncols() as f32) * ((j as f32) + 1. / 2.) * (v as f32)).cos()
                        );
                }
            }
            image[(u, v)] = acc;
        }
    }
    image
}

fn idct2(image_src : &na::DMat<f32>) -> na::DMat<f32> {
    let mut image = na::DMat::new_zeros(image_src.nrows(), image_src.ncols());
    for u in 0..image.nrows() {
        for v in 0..image.ncols() {
            image[(u, v)] = 1./4. * image_src[(0, 0)];
            for i in 1..image.nrows() {
                image[(u, v)] += 1./2. * image_src[(i, 0)];
            }
            for j in 1..image.ncols() {
                image[(u, v)] += 1./2. * image_src[(0, j)];
            }

            for i in 1..image.nrows() {
                for j in 1..image.ncols() {
                    image[(u, v)] += (
                        image_src[(i, j)] *
                        (PI / (image.nrows() as f32) * ((u as f32) + 1./2.) * (i as f32)).cos() *
                        (PI / (image.ncols() as f32) * ((v as f32) + 1./2.) * (j as f32)).cos()
                        );
                }
            }
            image[(u, v)] *= 2./(image.nrows() as f32) * 2./(image.ncols() as f32);
        }
    }
    image
}

fn main() {
    let img = load_image("test.csv".to_string());
    let img_f = filter_img(&img, 5, 1. / 25.);
    println!("IMG1: {:?}", img[(0, 0)]);
    println!("IMG2: {:?}", img_f[(0, 0)]);
    save_image(&img_f, "output.csv".to_string());
    let mut test_img = na::DMat::from_elem(6, 6, 255.);
    test_img[(5, 5)] = 5.;
    let img_dct = dct2(&test_img);
    let img_idct = idct2(&img_dct);
    println!("IMG:\n{:?}\n", test_img);
    println!("IMG DCT:\n{:?}\n", img_dct);
    println!("IMG IDCT:\n{:?}\n", img_idct);
}

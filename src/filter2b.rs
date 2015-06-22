extern crate nalgebra as na;
use na::{DMat};
use std::cmp::{min, max};
use std::ops::*;

pub fn filter_img(image_src : &na::DMat<f64>, mask_size : usize, mask_val : f64) -> na::DMat<f64> {
    let mut image = image_src.clone();
    assert_eq!(image.nrows(), image.ncols());
    let range = mask_size as i32 / 2;
    let mask = na::DMat::from_elem(mask_size, mask_size, mask_val);
    let size = image.nrows() as i32;

    let pos = |p| min(size - 1, max(0, p)) as usize;
    //let mask_pos = |p| (p + range) as usize;

    for x in 0..size {
        for y in 0..size {
            let mut acc = 0f64;
            for xd in -range..range+1 {
                for yd in -range..range+1 {
                    let val = image_src[(pos(x + xd), pos(y + yd))];
                    //println!("Parts: {}, {}", image_src[(pos(x + xd), pos(y + yd))], mask[(mask_pos(xd), mask_pos(yd))]);
                    //println!("...at: ({}, {}), ({}, {})", pos(x + xd), pos(y + yd), mask_pos(xd), mask_pos(yd));
                    //println!("...val: {}", val);
                    acc += val;
                }
            }
            acc = acc * mask_val;
            image[(x as usize, y as usize)] = acc;
        }
    }
    image
}

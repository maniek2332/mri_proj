extern crate nalgebra as na;
use na::{DMat};
use std::cmp::{min, max};
use std::ops::*;

pub fn filter_img(image_src : &na::DMat<f32>, mask_size : usize, mask_val : f32) -> na::DMat<f32> {
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
extern crate nalgebra as na;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Write;

pub fn load_image(path : String) -> na::DMat<f64> {
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
            let num = num_str.trim().parse::<f64>().unwrap();
            mat[(row, col)] = num;
            col += 1;
        }
        row += 1;
    }
    mat
}

pub fn save_image(img : &na::DMat<f64>, path : String) {
    let mut file = File::create(path).unwrap();
    for i in 0..img.nrows() {
        for j in 0..img.ncols() - 1 {
            write!(file, "{},", img[(i, j)]);
        }
        write!(file, "{}\n", img[(i, img.ncols() - 1)]);
    }
}

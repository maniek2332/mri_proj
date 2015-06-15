extern crate nalgebra as na;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;

fn main() {
    let mut mat = na::DMat::new_zeros(256, 256);

    let file = match File::open("test.csv") {
        Ok(file) => file,
        Err(..)  => panic!("panic"),
    };
    let mut br = BufReader::new(&file);

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
    println!("Hello, world! {}", mat[(0, 0)]);
    //println!("Mat: {:?}", mat);
}

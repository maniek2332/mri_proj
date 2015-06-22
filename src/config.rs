use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;

pub struct Config {
	pub use_snr: i32,
	pub ex_filter_type: i32,
	pub ex_window_size: usize,
	pub ex_iterations: i32,
	pub lpf_f: f64,
	pub lpf_f_SNR: f64,
	pub lpf_f_Rice: f64,
    pub input_filename: &'static str,
    pub input_snr: &'static str,
	pub output_filename_Gaussian: &'static str,
	pub output_filename_Rician: &'static str
}

pub fn load_config() -> Config {
    let mut cfg = Config {
        use_snr: 0,
        ex_filter_type: 1,
        ex_window_size: 5,
        ex_iterations: 10,
        lpf_f: 3.4,
        lpf_f_SNR: 1.2,
        lpf_f_Rice: 5.4,
        input_filename: "MR_noisy.csv",
        input_snr: "MR_SNR.csv",
        output_filename_Gaussian: "MR_Gaussian_Map.csv",
        output_filename_Rician: "MR_Rician_Map.csv"
    };

    let file = match File::open("config.ini") {
        Ok(file) => file,
        Err(..)  => panic!("Cannot open config file"),
    };
    let br = BufReader::new(&file);

    for line in br.lines().map(|line| line.unwrap()) {
        let mut col = 0;
        let mut words = line.split(" ");
        if words.clone().count() < 3 {
            continue;
        }
        let var_name = words.next().unwrap();
        let eq = words.next().unwrap();
        let var_val_str = words.next().unwrap();
        if eq != "=" {
            continue;
        }
        //println!("{} {} {}", var_name, eq, var_val_str);
        if var_name == "use_snr" {
            cfg.use_snr = var_val_str.trim().parse::<i32>().unwrap();
        }
        if var_name == "ex_filter_type" {
            cfg.ex_filter_type = var_val_str.trim().parse::<i32>().unwrap();
        }
        if var_name == "ex_window_size" {
            cfg.ex_window_size = var_val_str.trim().parse::<usize>().unwrap();
        }
        if var_name == "ex_iterations" {
            cfg.ex_iterations = var_val_str.trim().parse::<i32>().unwrap();
        }
        if var_name == "lpf_f" {
            cfg.lpf_f = var_val_str.trim().parse::<f64>().unwrap();
        }
        if var_name == "lpf_f_SNR" {
            cfg.lpf_f_SNR = var_val_str.trim().parse::<f64>().unwrap();
        }
        if var_name == "lpf_f_Rice" {
            cfg.lpf_f_Rice = var_val_str.trim().parse::<f64>().unwrap();
        }
        //if var_name == "input_filename" {
            //cfg.input_filename = var_val_str.clone().trim();
        //}
        //if var_name == "input_snr" {
            //cfg.input_snr = var_val_str.clone().trim();
        //}
        //if var_name == "output_filename_Gaussian" {
            //cfg.output_filename_Gaussian = var_val_str.clone().trim();
        //}
        //if var_name == "output_filename_Rician" {
            //cfg.output_filename_Rician = var_val_str.clone().trim();
        //}
    }
    cfg
}


pub static CONFIG : Config = Config {
	ex_filter_type: 1,
	ex_window_size: 5,
	ex_iterations: 10,
	lpf_f: 3.4,
	lpf_f_SNR: 1.2,
	lpf_f_Rice: 5.4,
    input_filename: "MR_noisy.csv",
	output_filename_Gaussian: "MR_Gaussian_Map.csv",
	output_filename_Rician: "MR_Rician_Map.csv"
};

pub struct Config {
	pub ex_filter_type: i32,
	pub ex_window_size: i32,
	pub ex_iterations: i32,
	pub lpf_f: f64,
	pub lpf_f_SNR: f64,
	pub lpf_f_Rice: f64,
    pub input_filename: &'static str,
	pub output_filename_Gaussian: &'static str,
	pub output_filename_Rician: &'static str
}

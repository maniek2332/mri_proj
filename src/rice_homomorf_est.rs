extern crate nalgebra as na;
use na::{DMat};
use std::ops::*;
use matlab_fun;
use matrix_math;
use filter2b;
use em_ml_rice2D;

pub fn compute_for_uknown_snr(mat : &DMat<f32>, lpf: f32, modo: i32) {
	let (m2, sigma_n) = em_ml_rice2D::compute(&mat,10, 3);
	let snr = matrix_math::div_matrix_by_matrix_each_value(&m2, sigma_n);
	
	retrun compute(&mat, &snr, lpf, modo);
}

pub fn compute(mat : &DMat<f32>, snr: &DMat<f32>, lpf: f32, modo: i32) {
	let const eq: f32 = 0.5772156649015328606;
	let (m2, sigma_n) = em_ml_rice2D::compute(&mat,10, 3);
	
	let m1= filter2b::filter_img(&mat, 5, 0.04);
}

Rn=abs(In-M1);
lRn=log(Rn.*(Rn~=0)+0.001.*(Rn==0));
LPF2=lpf((lRn),LPF);
Mapa2=exp(LPF2);
MapaG=Mapa2.*2./sqrt(2).*exp(-psi(1)./2);

%Rician-------------------------
if Modo==1
    LocalMean=M1;
elseif Modo==2
    LocalMean=M2;
else
    LocalMean=0;
end

Rn=abs(In-LocalMean);
lRn=log(Rn.*(Rn~=0)+0.001.*(Rn==0));
LPF2=lpf((lRn),LPF);
Fc1=correct_rice_gauss(SNR);
LPF1=LPF2-Fc1;
LPF1=lpf((LPF1),LPF+2,2);
Mapa1=exp(LPF1);
MapaR=Mapa1.*2./sqrt(2).*exp(-psi(1)./2);
//use std::f32;
use std::f64::consts::PI as PI;
const SIGMA_T: f64 = 6.65e-25; // cm ^2
const STEF_BOLTZ: f64 = 5.6704e-5;
const c: f64 = 3e10; //cm s^-1????
const BOLTZ_CONST: f64 = 1.3806504e-16; //erg/K
const PLANCK_CONST: f64 = 6.63e-27; //ergs
const MASS_ELECTRON : f64= 9.109e-28; //g
//alpha=1.
//d_1=2.32e14
//theta=N.pi/2.

fn main(){
	println!("yolo world");
	let B = 1.8; // Gauss
	let v_b = 2.8e6 * B;
	let dl = 2230.0 * 3.08568e24;
	let z = 0.409; //redshift
	let freq_blr = 2.47e15;
	let d_1 = 2.32e14;
	let alpha = 1.;

	//Parameters
	let n1 = 1.2;
	let n2 = 3.2;
	let gamma_break = 500; //gamma for break frequency where electron distribution drops 900
	let R = 3.0e15; //cm size scale of blob
	let no = 1e4;
	let Ld = 1.7e45; //erg s^-1
	let bulk_lorentz = 4;
	//doppler_factor=3
	let gamma_max = 7.5e4;
	let R_diss = 3.0e16;
	let theta2 = 3. * PI; //3.*N.pi/180.;
	let gamma_min = 1.;
	let gamma_step = 0.1;

	//loggammastep=0.01 #0.01
	let beta_factor = (1. - (1./bulk_lorentz.pow(2))).sqrt(); //N.sqrt(1.-(1./(bulk_lorentz**2.)))
	let doppler_factor = 1./(bulk_lorentz * (1. - beta_factor * (3. * PI).cos()));
	let fBLR = 0.1;
	let Rblr = 1.0e17 * (Ld/1.0e45).sqrt();
	let optical_depth = SIGMA_T * R * no;
	//let logvmin=N.round(N.log10(vb),2)#8.41 #6.41 #8.41
	//vmin=10**(logvmin)
	//logblrvmin = N.log10(vblr)
	//logvcompmax=30. #N.log10(vb)+4*loggammamax#26.41
	//vmax=10**(logvcompmax)
	let freq_step = 0.1;
	//Radiation Energy Densities
	let uB = B.pow(2)/(8. * PI);

	let num_bins: u64 = (gamma_max - gamma_min / gamma_step).round();
	//MAKE SYNCHROTRON SPECTRUM
	let emissivity_synchrotron: [f64; num_bins] = [0; num_bins]; // needs to be vector
	let emissivity_compton: [f64; num_bins] = [0; num_bins]; // num_bins not known compile time

	fn gamma_value(bin_number: usize, gamma_min: f64, gamma_max: f64, num_bins: u64) -> f64{
		gamma_min + (gamma_max - gamma_min) * (bin_number as f64/ num_bins as f64)
	}

	fn n_gamma(gamma_value: f64, no: f64, n1: f64, n2: f64, gamma: f64, gamma_break: f64) -> f64{
		no * gamma_value.pow(-n1) * (1. + (gamma / gamma_break)).pow(n1 - n2)
	}

	for bin in 0..num_bins{
		let gamma_val = gamma_value(bin, gamma_min, gamma_max, num_bins);
	    let k = SIGMA_T * c * uB /(6. * PI * v_b);
	    let freq_synchrotron = gamma_val.pow(2) * v_b; //need to make sure is cloning gamma_val
	    emissivity_synchrotron[bin] = k * gamma_val * n_gamma(gamma_val,
			no, n1, n2, gamma_val, gamma_break);
	}

	//SSA CUT_OFF
	let freq_peak = ((2. * d_1 * no * (B.pow(1.5 + alpha)) * R)/((PI / 2.).sin() * 1.)).pow(2./(5.+2.*alpha)); //the 1 is optical depth at v(ssa)
	/*peak=int(round((N.log10(vpeak)-logvmin)/vstep))
	kssa=js[peak]/(vpeak**(5./2.))
	for ip in range (peak):
	    js[ip]=0.
	    #js[ip]=kssa*(freq[ip]**(5./2.))
	*/
}

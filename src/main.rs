//use std::f32;
extern crate gnuplot;
use gnuplot::{Figure, Caption, Color};
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
	let B: f64 = 1.8; // Gauss
	let v_b = 2.8e6 * B;
	let v_max: f64 = 1.0e14;
	let v_step: f64 = 0.1;
	let dl: f64 = 2230.0 * 3.08568e24;
	let z: f64 = 0.409; //redshift
	let freq_blr: f64 = 2.47e15;
	let d_1: f64 = 2.32e14;
	let alpha: f64 = 1.;

	//Parameters
	let n1: f64 = 1.2;
	let n2: f64 = 3.2;
	let gamma_break: f64 = 500.0; //gamma for break frequency where electron distribution drops 900
	let R = 3.0e15; //cm size scale of blob
	let no: f64 = 1.0e4;
	let Ld: f64 = 1.7e45; //erg s^-1
	let bulk_lorentz: f64 = 4.0;
	//doppler_factor=3
	let gamma_max: f64 = 7.5e4;
	let R_diss: f64 = 3.0e16;
	let theta2 = 3. * PI; //3.*N.pi/180.;
	let gamma_min: f64 = 1.;
	let gamma_step: f64 = 0.1;

	//loggammastep=0.01 #0.01
	let beta_factor = (1. - (1./bulk_lorentz.powf(2.0))).sqrt(); //N.sqrt(1.-(1./(bulk_lorentz**2.)))
	let doppler_factor = 1./(bulk_lorentz * (1. - beta_factor * (3. * PI).cos()));
	let fBLR: f64 = 0.1;
	let Rblr = 1.0e17 * (Ld/1.0e45).sqrt();
	let optical_depth = SIGMA_T * R * no;
	//let logvmin=N.round(N.log10(vb),2)#8.41 #6.41 #8.41
	//vmin=10**(logvmin)
	//logblrvmin = N.log10(vblr)
	//logvcompmax=30. #N.log10(vb)+4*loggammamax#26.41
	//vmax=10**(logvcompmax)
	let freq_step: f64 = 0.1;
	//Radiation Energy Densities
	let uB = B.powf(2.0)/(8. * PI);

	let num_bins: usize = (gamma_max - gamma_min / gamma_step).round() as usize;
	let v_step: f64 = (v_max - v_b / num_bins as f64);
	//MAKE SYNCHROTRON SPECTRUM
	let mut emissivity_synchrotron = vec![(0.0, 0.0); num_bins as usize]; // needs to be vector
	let mut emissivity_compton = vec![(0.0, 0.0); num_bins as usize]; // num_bins not known compile time

	fn gamma_value(bin_number: usize, gamma_min: f64, gamma_max: f64, num_bins: usize) -> f64{
		gamma_min + (gamma_max - gamma_min) * (bin_number as f64/ num_bins as f64)
	}

	fn n_gamma(gamma_value: f64, no: f64, n1: f64, n2: f64, gamma: f64, gamma_break: f64) -> f64{
		no * gamma_value.powf(-n1) * (1. + (gamma / gamma_break)).powf(n1 - n2)
	}

	//fn emissivity_to_luminosity(emissivity: Vec<f64>) -> Vec<f64>{}

	for bin in 0..num_bins{
		let gamma_val = gamma_value(bin, gamma_min, gamma_max, num_bins);
	    let k = SIGMA_T * c * uB /(6. * PI * v_b);
	    let freq_synchrotron = gamma_val.powf(2.0) * v_b; //need to make sure is cloning gamma_val
	    emissivity_synchrotron[bin] = (freq_synchrotron, k * gamma_val * n_gamma(gamma_val,
			no, n1, n2, gamma_val, gamma_break));
	}

	let emiss_to_flux = |x: Vec<(f64, f64)>| {x.iter().map(|y| y.1 * 4. * PI * R.powf(3.) / (3. * dl.powf(2.)));
		x
	};

	//println!("{:?}", emissivity_synchrotron);
	let flux_synchrotron = emiss_to_flux(emissivity_synchrotron);
	//println!("{:?}", flux_synchrotron);
	//SSA CUT_OFF
	let freq_peak = ((2. * d_1 * no * (B.powf(1.5 + alpha)) * R)/((PI / 2.).sin() * 1.)).powf(2./(5.+2.*alpha)); //the 1 is optical depth at v(ssa)

	fn bin_it(tuples: Vec<(f64, f64)>, min_bin: f64, max_bin: f64, step: f64) -> Vec<(f64, f64)>{
		println!("{}  max bin", max_bin);
		let mut binned_results = vec!();
		let mut current_bin = min_bin - step;
		for (freq, flux) in tuples{
			//println!("{}, {}", freq, current_bin);
			match freq{
				f if f >= max_bin => {
					println!("breaking");
						break
				},
				f if f >= current_bin + step => { // assuming cannot jump multiple bins
					println!("new bin");
					while f < current_bin + step{
						current_bin += step
					};
					binned_results.push((current_bin + step / 2., flux))
				},
				_ => {
					println!("add to bin");
					let mut bin = binned_results.iter().find(|&&x| x.0 == current_bin + step / 2.).unwrap();
					//println!(" flux {:?}", bin.1)
					bin = &(bin.0, bin.1 + flux) //unsafe if cant find
				}
			}
			//println!("binned res {:?}", binned_results);
		}
		binned_results
	};

	// x: log frequency, y: log frequency * log flux
	let result_tuple = |tuples: Vec<(f64, f64)>| -> (Vec<f64>, Vec<f64>){
		(tuples.iter().map(|x| x.0).collect::<Vec<f64>>(),
		 tuples.iter().map(|x| x.1 * x.0).collect::<Vec<f64>>())
	};
	let binned_results = bin_it(flux_synchrotron, v_b, v_max, v_step);
	let (frequencies, fluxes) = result_tuple(binned_results);
	//let y = [3u32, 4, 5];
	let mut fg = Figure::new();
	fg.axes2d()
	.lines(&frequencies, &fluxes, &[Caption("A line"), Color("black")]);
	fg.show();


	/*peak=int(round((N.log10(vpeak)-logvmin)/vstep))
	kssa=js[peak]/(vpeak**(5./2.))
	for ip in range (peak):
	    js[ip]=0.
	    #js[ip]=kssa*(freq[ip]**(5./2.))
	*/
}

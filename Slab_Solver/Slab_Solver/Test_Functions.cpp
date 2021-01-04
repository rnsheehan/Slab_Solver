#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::vec_mat_test()
{
	// Check to ensure that the vec mat product is correct

	int rs = 3, cs = 3;

	std::vector<std::vector<std::complex<double>>> M; 
	std::vector<std::complex<double>> X(rs, zero); 
	std::vector<std::complex<double>> B(rs, zero);

	M = vecut::zero_cmat(rs, cs); 

	M[0][0].real(4); M[0][1].real(2); M[0][2].real(-1); 
	M[1][0].real(0); M[1][1].real(7); M[1][2].real(5); 
	M[2][0].real(3); M[2][1].real(-5); M[2][2].real(-6);

	X[0].real(4); X[1].real(-2); X[2].real(1); 

	B = vecut::cmat_cvec_product(M, X); 

	for (int i = 0; i < rs; i++)
		std::cout << B[i] << "\n"; 
}

void testing::slab_wg_neff_calc()
{
	// Run a three layer slab waveguide neff calculation and check the results
	// R. Sheehan 2 - 5 - 2018

	double W = 1.5; 
	double WL = 1.55; 
	double Nc = 3.38; 
	double Ns = 3.17; 
	double Ncl = 1.0; 

	slab_tl_neff sl_obj; 

	sl_obj.set_params(W, WL, Nc, Ns, Ncl); 

	sl_obj.neff_search(TE); 

	std::cout << "Complete\n"; 
}

void testing::slab_wg_mode_calc()
{
	// Run a three layer slab waveguide mode calculation and check the results
	// R. Sheehan 2 - 5 - 2018

	double W, WL, Nc, Ns, Ncl; 

	W = 2.0;  WL = 1.55;
	Nc = 3.38; Ns = 1.0; Ncl = 1.0;

	//W = 1.7;  WL = 6.158;
	//Nc = 3.54102488; Ns = 1.75275048;

	/*W = 2.0; WL = 1.55; 
	Nc = 3.27562; Ns = 3.27446;*/

	//W = 2.5; WL = 1.55; 
	//Nc = 3.293; Ns = 3.29152;

	/*W = 0.3; WL = 1.535; 
	Nc = 2.00704; Ns = 1.44446; Ncl = 1.0;*/ 

	/*W = 0.3; WL = 1.55;
	Nc = 2.00704; Ns = 1.44428; Ncl = 1.0;*/

	slab_tl_mode sl_obj;

	sl_obj.set_params(W, WL, Nc, Ns, Ncl);

	sl_obj.compute_neff(TE); 

	std::cout << "TE Fund. Mode Field value at origin: ";
	std::cout << sl_obj.TE_TM(0.0, 0, TE) << "\n\n"; 

	sl_obj.compute_neff(TM);

	std::cout << "TM Fund Mode Field value at origin: ";
	std::cout << sl_obj.TE_TM(0.0, 0, TM) << "\n\n";

	int N = 501; 
	double Lz = 3.5*W; 
	std::string stor = ""; 

	sl_obj.output_all_stats(stor);

	sl_obj.output_modes(TE, N, Lz, stor); 

	std::cout << "Complete\n";
}

void testing::fl_slab_wg_neff_calc()
{
	// Compute the effective index in a four layer slab waveguide
	// R. Sheehan 5 - 10 - 2018

	double W, Wr, WL, Nc, Ns, Nr, Ncl; 

	// Case A => Field Oscillating in Core and Ridge
	// For there to be a solution one has to have ns <= ncl < nr < nc
	//W = 1.0; Wr = 1.0; WL = 1.55; 
	//W = 0.5; Wr = 0.5; WL = 1.55;
	//W = 1.0; Wr = 0.5; WL = 1.55;
	//Nc = 3.38; Ns = 1.0; Nr = 3.17; Ncl = 1.0;

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)
	//W = 0.5; Wr = 0.5; WL = 1.55; // n_{TM} = 3.274464
	//W = 0.5; Wr = 1.0; WL = 1.55; // n_{TM} = 3.27562
	//W = 0.6; Wr = 0.4; WL = 1.55; // n_{TM} = 3.29152
	//W = 0.6; Wr = 0.9; WL = 1.55; // n_{TM} = 3.293
	//Ns = 3.17; Nc = 3.38; Nr = 3.17; Ncl = 1.0;
	
	W = 0.45; Wr = 0.3; WL = 1.55; // no modes present case B
	Ns = 3.17; Nc = 3.38; Nr = 3.341; Ncl = 1.0; // no modes present case B

	//W = 0.6; Wr = 0.5; WL = 1.55;
	//Ns = 3.17; Nc = 3.38; Nr = 3.19; Ncl = 1.0; 

	slab_fl_mode_A fl_obj;

	//slab_fl_mode_B fl_obj;

	fl_obj.set_params(W, Wr, WL, Nc, Ns, Ncl, Nr); 

	fl_obj.compute_neff(TM); 
}

void testing::fl_slab_wg_mode_calc()
{
	// Compute the mode profile in a four layer slab waveguide
	// R. Sheehan 15 - 12 - 2020

	bool pol = TE; 

	double W, Wr, WL, Nc, Ns, Nr, Ncl;

	// Case A => Field Oscillating in Core and Ridge
	// For there to be a solution one has to have ns <= ncl < nr < nc
	W = 1.0; Wr = 1.5; WL = 1.55; 
	//W = 0.5; Wr = 0.5; WL = 1.55;
	//W = 1.0; Wr = 0.5; WL = 1.55;
	Nc = 3.38; Ns = 3.0; Nr = 3.17; Ncl = 1.0;

	// Case B: Field Oscillating in Core Only
	// For there to be a solution one has to have ncl < nm < nc, where nm = Max(nr,ns)
	//W = 0.5; Wr = 0.5; WL = 1.55; // n_{TM} = 3.274464
	//W = 0.5; Wr = 1.0; WL = 1.55; // n_{TM} = 3.27562
	//W = 0.6; Wr = 0.4; WL = 1.55; // n_{TM} = 3.29152
	//W = 0.6; Wr = 0.9; WL = 1.55; // n_{TM} = 3.293
	//Ns = 3.17; Nc = 3.38; Nr = 3.17; Ncl = 1.0;

	//W = 0.45; Wr = 0.3; WL = 1.55; // no modes present case B
	//Ns = 3.17; Nc = 3.38; Nr = 3.341; Ncl = 1.0; // no modes present case B

	//W = 0.6; Wr = 0.5; WL = 1.55;
	//Ns = 3.17; Nc = 3.38; Nr = 3.19; Ncl = 1.0; 

	slab_fl_mode_A fl_obj;

	//slab_fl_mode_B fl_obj;

	fl_obj.set_params(W, Wr, WL, Nc, Ns, Ncl, Nr);

	fl_obj.bracket_roots(pol, true); 

	fl_obj.compute_neff(pol); 

	// compute field value at x = 0.0
	std::cout << "\nField value at x = 0: " << fl_obj.TE_TM(0.0, 0, pol) << "\n";

	int N = 201;
	double Lz = 3.5 * (W + Wr);
	std::string stor = "";
	
	fl_obj.output_modes(pol, N, Lz, stor);
}

void testing::coupled_slab_wg_calc()
{
	// Compute the coupling coefficient for a pair of identical slab waveguides
	// R. Sheehan 1 - 10 - 2018

	// W: waveguide width, D: waveguide pitch, WL: wavelength, Nc: core RI, Ns: substrate RI
	double W, D, WL, Nc, Ns;

	// Code produces answer that matches with Okamoto example values
	// see Okamoto page 178, kappa = 0.39 mm^{-1}, L_{c} = 4 mm
	//W= 6; D = 12; WL = 1.5171; Nc = 1.46; Ns = 1.455;

	// Other examples

	// Si wire of height H = 0.22 um has neff^{TE} = 2.81 along vertical for nc = 3.45, ns = 1.45 and ncl = 1.0 at l = 1.55
	// Coupling coefficient between two Si wires of width W = 450 nm, D = 300 nm is kappa = 2.4 um^{-1}, L_{c} = 653 nm
	// Calculation is correct when TE polarisation is assumed in slab
	W = 0.45; D = 0.3; WL = 1.55; Nc = 2.81; Ns = 1.0;

	double kappa;

	//W = 4.0; WL = 6.09996577; Nc = 3.377661931; Ns = 3.148794556;

	//W = 2; WL = 1.535; Nc = 1.63701; Ns = 1;
	//W = 1.0; WL = 1.535; Nc = 1.475254387; Ns = 1;

	coupled_slab_tl_neff sl_obj;

	std::string filename = "Coupling_Coeff_W_" + template_funcs::toString(W, 2) + ".txt"; 

	std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

	D = 0.3; 

	for (int i = 0; i < 6; i++) {

		sl_obj.set_params(D, W, WL, Nc, Ns);

		kappa = sl_obj.compute_coupling_coeff(TM);

		std::cout << "Separation = " << D << " um\n"; // Is separation the same as pitch? No

		std::cout << "Coupling coefficient = " << kappa << " um^{-1}\n";

		std::cout << "Coupling length = " << (PI / (2.0*kappa)) << " um\n";

		std::cout << "\n\n";

		write << std::setprecision(10) << D << " , " << kappa << "\n"; 

		D += 0.1; 

	}

	write.close(); 
}

void testing::coupled_slab()
{
	// test the operation of the coupled slab calculator
	// R. Sheehan 4 - 9 - 2020

	// Si wire of height H = 0.22 um has neff^{TE} = 2.81 along vertical for nc = 3.45, ns = 1.45 and ncl = 1.0 at l = 1.55
	// Coupling coefficient between two Si wires of width W = 450 nm, D = 300 nm is kappa = 2.4 um^{-1}, L_{c} = 653 nm
	// Calculation is correct when TE polarisation is assumed in slab
	// Citation?
	//double Wa = 0.45, Wb = 0.45, l = 1.55, ncorea = 2.81, ncoreb = 2.81, nsub = 1.0, pitch = 0.75; 

	//double Wa = 1.5, Wb = 2.0, l = 1.55, ncorea = 3.38, ncoreb = 3.4, nsub = 3.17, pitch = 1.85;

	// Test Problem from Chuang
	// Strong coupling between pair of coupled Ti-diffused LiNbO3 WG
	// C = 0.168, kab = kba ~ 0.0027 for pitch = 3.9 um
	double Wa = 2.0, Wb = 2.0, l = 1.06, ncorea = 2.2, ncoreb = 2.2, nsub = 2.19, pitch = 3.9;
	double dn = +0.0;
	ncorea += 0.5 * dn; ncoreb -= 0.5 * dn;

	std::cout << "Domains\n"; 
	std::cout << "-a <= x <= +a: " << -0.5 * Wa << " <= x <= " << 0.5 * Wa << "\n"; 
	std::cout << "D-b <= x <= D+b: " << pitch-0.5 * Wb << " <= x <= " << pitch+0.5 * Wb << "\n\n"; 

	coupled_slabs wg_pair; 

	wg_pair.set_params(Wa, Wb, l, ncorea, ncoreb, nsub); 

	wg_pair.output_modes(pitch); 

	wg_pair.compute_coefficients(pitch, true); 

	wg_pair.propagate(1200, 10, 0, 1, false); 
}

void testing::coupled_slab_benchmark()
{
	// Benchmark test for the coupled slab class
	// Calculation is based on data from the paper
	// S-L. Chuang, ``Application of the strongly coupled-mode theory to integrated optical devices'', IEEE J. Quant. Electron., 23 (5), 1987
	// It involves the calculation of the overlap integrals and coupling coefficients for two strongly-coupled Ti LiNbO3 waveguides separated by
	// fixed distance, the RI in each waveguide is variable
	// R. Sheehan 9 - 10 - 2020

	// Test Problem from Chuang
	// Strong coupling between pair of coupled Ti-diffused LiNbO3 WG
	// C = 0.168, kab = kba ~ 0.0027 for pitch = 3.9 um
	double Wa = 2.0, Wb = 2.0, l = 1.06, ncorea = 2.2, ncoreb = 2.2, nsub = 2.19, pitch = 3.9;

	double dn_start = -0.004; 
	double dn_stop = 0.0045; 
	double dn = +0.0005;

	std::string filename = "Chuang_Benchmark.txt"; 
	std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

	while (dn_start < dn_stop) {
		ncorea = 2.2; 
		ncoreb = 2.2;
		ncorea += 0.5 * dn_start; ncoreb -= 0.5 * dn_start;

		coupled_slabs wg_pair;

		wg_pair.set_params(Wa, Wb, l, ncorea, ncoreb, nsub);

		wg_pair.compute_coefficients(pitch, false);

		write << std::setprecision(10) << dn_start << " , " << wg_pair.get_CAB() << " , " << wg_pair.get_kab() << " , " << wg_pair.get_kba() << " , " << wg_pair.get_async() << "\n"; 

		dn_start += dn; 
	}

	write.close(); 
	
}
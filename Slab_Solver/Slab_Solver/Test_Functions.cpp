#ifndef ATTACH_H
#include "Attach.h"
#endif

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

	double W = 0.22;
	double WL = 1.55;
	double Nc = 3.45;
	double Ns = 1.45;
	double Ncl = 1.0;

	slab_tl_mode sl_obj;

	sl_obj.set_params(W, WL, Nc, Ns, Ncl);

	sl_obj.compute_neff(TE); 

	int N = 501; 
	double Lz = 3.0; 
	std::string stor = ""; 

	sl_obj.output_all_stats(stor);

	sl_obj.output_modes(TE, N, Lz, stor); 

	std::cout << "Complete\n";
}

void testing::coupled_slab_wg_calc()
{
	// Compute the coupling coefficient for a pair of identical slab waveguides
	// R. Sheehan 1 - 10 - 2018

	double W, D, WL, Nc, Ns; 

	// Code produces answer that matches with Okamoto example values
	// see Okamoto page 178, kappa = 0.39 mm^{-1}, L_{c} = 4 mm
	//W= 6; D = 12; WL = 1.5171; Nc = 1.46; Ns = 1.455;

	// Other examples

	// Si wire of height H = 0.22 um has neff^{TE} = 2.81 along vertical for nc = 3.45, ns = 1.45 and ncl = 1.0 at l = 1.55
	// Coupling coefficient between two Si wires of width W = 450 nm, D = 300 nm is kappa = 2.4 um^{-1}, L_{c} = 653 nm
	W = 0.45; D = 0.3; WL = 1.55; Nc = 2.81; Ns = 1.0;

	coupled_slab_tl_neff sl_obj;

	sl_obj.set_params(D, W, WL, Nc, Ns);

	double kappa; 
	
	kappa = sl_obj.compute_coupling_coeff(TE);

	std::cout << "Coupling coefficient = " << kappa << " um^{-1}\n"; 

	std::cout << "Coupling length = " << (PI / (2.0*kappa)) << " um\n"; 

	std::cout << "Complete\n";
}
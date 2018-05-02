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

	double W = 1.5;
	double WL = 1.55;
	double Nc = 3.38;
	double Ns = 3.17;
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
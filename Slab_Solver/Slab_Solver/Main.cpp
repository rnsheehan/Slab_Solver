
#ifndef ATTACH_H
#include "Attach.h"
#endif

// This project implements the solvers for three and four layers slabs
// I want to implement the solvers in a more inherited manner
// R. Sheehan 2 - 5 - 2018

int main()
{

	//testing::vec_mat_test(); 

	//testing::slab_wg_mode_calc(); 

	//testing::fl_slab_wg_neff_calc(); 

	//testing::fl_slab_wg_mode_calc_a(); 

	testing::fl_slab_wg_mode_calc_b(); 

	//testing::coupled_slab_wg_calc(); 

	//testing::coupled_slab(); 
	
	//testing::coupled_slab_benchmark(); 

	std::cout << "Press return to close\n";
	std::cin.get(); 
	return 0; 
}
#ifndef ATTACH_H
#include "Attach.h"
#endif

// This project implements the solvers for three and four layers slabs
// I want to implement the solvers in a more inherited manner
// R. Sheehan 2 - 5 - 2018

int main()
{

	//testing::slab_wg_mode_calc(); 

	testing::coupled_slab_wg_calc(); 

	std::cout << "Press return to close\n";
	std::cin.get(); 
	return 0; 
}
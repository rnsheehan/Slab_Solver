#ifndef TEST_FUNCTIONS_H
#define TEST_FUNCTIONS_H

// Namespace containing functions used to check that calculations are being done correctly
// R. Sheehan 2 - 5 - 2018

namespace testing {

	void vec_mat_test();
	
	void slab_wg_neff_calc(); 

	void slab_wg_mode_calc(); 

	void fl_slab_wg_neff_calc(); 

	void fl_slab_wg_mode_calc_a(); 

	void fl_slab_wg_mode_calc_b(); 

	void coupled_slab_wg_calc(); 

	void coupled_slab(); 

	void coupled_slab_benchmark();

}

#endif

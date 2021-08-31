#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::energy_unit_conversion()
{
	// what sort of values do we get when talking about eV and J
	// R. Sheehan 19 - 8 - 2021

	for (int i = 1; i < 11; i++) {
		std::cout << i << " eV = " << template_funcs::convert_ev_J((double)(i)) << " J\n"; 
	}
}

void testing::potential_step()
{
	// test the potential step calculation
	// R. Sheehan 19 - 8 - 2021

	pot_step classical; 

	classical.set_params(M_ELECTRON_KG, 1, 1.1); 

	classical.compute_wavefunction("Step_Solution.txt"); 
}

void testing::infinite_well()
{
	// test the implementation of the code for the infinite square well
	// R. Sheehan 24 - 10 -  2016

	double length = 2.0; // length scale is assumed to be nm
	double pos = 0.7; 

	inf_well the_well(length, M_ELECTRON_KG); 

	for(int i = 1; i < 5; i++){
		std::cout<<"Energy of level "<<i<<": "<< the_well.energy_eigenvalue(i) <<"\n"; 
	}

	std::cout<<"\n";

	std::cout<<"Probability of being located at position x = 0.7: "<<template_funcs::DSQR( the_well.energy_eigenfunction(1, pos) )<<"\n";
}
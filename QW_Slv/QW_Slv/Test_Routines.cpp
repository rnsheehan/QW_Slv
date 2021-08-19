#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::infinite_well()
{
	// test the implementation of the code for the infinite square well
	// R. Sheehan 24 - 10 -  2016

	double length = 2.0; // length scale is assumed to be nm
	double pos = 0.7; 

	inf_well the_well(length, M_ELECTRON_KG); 

	for(int i = 1; i < 5; i++){
		std::cout<<"Energy of level "<<i<<": "<< the_well.energy_eigenvalue(i, false) <<"\n"; 
	}

	std::cout<<"\n";

	std::cout<<"Probability of being located at position x = 0.7: "<<template_funcs::DSQR( the_well.energy_eigenfunction(1, pos) )<<"\n";
}
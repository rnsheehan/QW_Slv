#ifndef ATTACH_H
#include "Attach.h"
#endif

// stop fucking about with the choice of units and go with eV only!! R. Sheehan 11 - 1 - 2018

int main(int argc, char* argv[])
{
	//testing::energy_unit_conversion(); 

	testing::potential_step(); 

	//testing::infinite_well(); 

	std::cout<<"Press enter to close\n"; 
	std::cin.get(); 

	return 0; 
}
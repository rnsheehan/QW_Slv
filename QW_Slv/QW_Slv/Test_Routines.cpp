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

	classical.set_params(M_ELECTRON_KG, 2, 1, true); 

	classical.compute_wavefunction("Step_Solution_E_gr_V.txt"); 
	//classical.compute_wavefunction("Step_Solution_E_le_V.txt"); 
}

void testing::potential_step_ratio()
{
	// Look at the variation of R, T as ratio of E to V changes
	// R. Sheehan 31 - 8 - 2021

	try {

		std::string filename = "Step_Probabilities.txt";

		std::ofstream write;

		write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

		if (write.is_open()) {
			// perform the calculations
			// energies in units of eV
			double step_height = 1.0;
			double particle_energy = 1.001;
			double delta_energy = 0.001;
			double ratio; 
			double particle_mass = M_ELECTRON_KG;

			pot_step the_step;

			for (int i = 1; i <= 1000; i++) {
				ratio = particle_energy / step_height; 

				the_step.set_params(particle_mass, particle_energy, step_height); // compute T, R probabilities

				write << std::setprecision(10) << ratio << " , " << the_step.get_T() << " , " << the_step.get_R() << "\n"; 

				particle_energy += delta_energy; // update particle energy
			}

			write.close(); 
		}
		else {
			std::string reason = "Error: void testing::potential_step_E_gr_V()\n";
			reason += "Could not open file: " + filename + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void testing::potential_barrier()
{
	// test the implementation of the code for the potential barrier solution
	// R. Sheehan 31 - 8 - 2021

	pot_barr classical;

	classical.set_params(M_ELECTRON_KG, 1, 1.1, 1, true); 

	classical.compute_wavefunction("Barrier_Solution.txt");
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
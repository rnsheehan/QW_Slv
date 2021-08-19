#ifndef ATTACH_H
#include "Attach.h"
#endif

pot_step::pot_step()
{
	// default constructor
	m = E = V = p1 = p2 = T = R = A = B = C = D = 0.0; 
	CC = DD = zero; 
}

pot_step::pot_step(double particle_mass, double particle_energy, double step_height)
{
	// primary constructor
	set_params(particle_mass, particle_energy, step_height); 
}

void pot_step::set_params(double particle_mass, double particle_energy, double step_height)
{
	// assign values to the parameters for the potential step calculation

	try {
		bool c1 = particle_mass > 0.0 ? true : false;
		bool c2 = particle_energy > 0.0 ? true : false;
		bool c3 = step_height > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			m = particle_mass; 
			E = particle_energy; 
			V = step_height; 
			p1 = sqrt(2.0 * m * E); 
			
			if (E > V) {
				//particle energy greater than step height
				p2 = sqrt(2.0 * m * (E - V));				
				double psum = p1 + p2; 
				double pdiff = p2 - p1; 
				B = 1;
				C = pdiff / psum;
				A = (2.0 * p1) / psum; 
				T = (4.0 * p1 * p2) / (template_funcs::DSQR(psum)); 
				R = (template_funcs::DSQR(pdiff)) / (template_funcs::DSQR(psum)); 
				std::cout << "T = " << T << ", R = " << R << ", T+R = " << T + R << "\n"; 
			}
			else {
				//particle energy less than step height
				double delta = V - E; 
				p2 = sqrt(2.0 * m * delta);
				B = 1; 
				D = sqrt(delta / E); 
				CC = (B / 2) * (one + eye * D); 
				DD = (B / 2) * (one - eye * D); 
				R = 1.0; 
				T = 0.0; 
			}
		}
		else {
			std::string reason = "Error: void pot_step::set_params(double particle_mass, double particle_energy, double step_height)\n";
			if (!c1) reason += "particle_mass is not positive\n";
			if (!c2) reason += "particle_energy is not positive\n";
			if (!c3) reason += "step_height is not positive\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}
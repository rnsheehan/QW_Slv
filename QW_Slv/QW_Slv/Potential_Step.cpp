#ifndef ATTACH_H
#include "Attach.h"
#endif

pot_step::pot_step()
{
	// default constructor
	m = E = V = p1 = p2 = T = R = A = B = C = D = 0.0; 
	CC = DD = zero; 
	params_defined = high = false; 
}

pot_step::pot_step(double particle_mass, double particle_energy, double step_height)
{
	// primary constructor
	set_params(particle_mass, particle_energy, step_height); 
}

void pot_step::set_params(double particle_mass, double particle_energy, double step_height, bool loud)
{
	// assign values to the parameters for the potential step calculation
	// particle mass in units of kg
	// energies in units of eV

	try {
		bool c1 = particle_mass > 0.0 ? true : false;
		bool c2 = particle_energy > 0.0 ? true : false;
		bool c3 = step_height > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			m = particle_mass; 
			E = template_funcs::convert_ev_J(particle_energy); 
			V = template_funcs::convert_ev_J(step_height);
			p1 = sqrt(2.0 * m * E); 
			
			if (E > V) {
				//particle energy greater than step height
				high = true; 
				p2 = sqrt(2.0 * m * (E - V));				
				double psum = p1 + p2; 
				double pdiff = p2 - p1; 
				B = 1;
				C = ( pdiff * B ) / psum;
				A = (2.0 * p1 * B) / psum; 
				T = (4.0 * p1 * p2) / (template_funcs::DSQR(psum)); 
				R = (template_funcs::DSQR(pdiff)) / (template_funcs::DSQR(psum)); 
				
				t1 = (eye * p1 * 1.0e-9) / H_BAR_J; // scale length to nm
				t2 = (eye * p2 * 1.0e-9) / H_BAR_J; // scale length to nm

				if (loud) {
					std::cout << "B = " << B << " , C = " << C << " , A = " << A << "\n";
					std::cout << "p1 = " << p1 << " , p2 = " << p2 << "\n";
					std::cout << "t1 = " << t1 << " , t2 = " << t2 << "\n";
					std::cout << "T = " << T << ", R = " << R << ", T+R = " << T + R << "\n";
				}
			}
			else {
				//particle energy less than step height
				high = false; 
				double delta = V - E; 
				p2 = sqrt(2.0 * m * delta);
				B = 1; 
				D = sqrt(delta / E); 
				CC = (B / 2) * (one + eye * D); 
				DD = (B / 2) * (one - eye * D); 
				R = 1.0; 
				T = 0.0; 
				
				t1 = (eye * p1 * 1.0e-9) / H_BAR_J;
				t2 = (-1.0e-9 * p2) / H_BAR_J;


				if (loud) {
					std::cout << "CC = " << CC << " , DD = " << DD << " , B = " << B << "\n";
					std::cout << "p1 = " << p1 << " , p2 = " << p2 << "\n";
					std::cout << "t1 = " << t1 << " , t2 = " << t2 << "\n";
					std::cout << "T = " << T << ", R = " << R << ", T+R = " << T + R << "\n";
				}
			}

			params_defined = true; 
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

std::complex<double> pot_step::wavefunction(double position)
{
	// compute the value of the particle wavefunction for the potential step problem
	// R. Sheehan 20 - 8 - 2021

	try {
	
		if (params_defined) {
			if (high) {
				// Case E > V
				if (position < 0.0) {
					// incoming sine wave
					return ( ( B * exp(t1 * position) ) + ( C * exp(-1.0 * t1 * position) ) ); 
				}
				else {
					// outgoing sine wave
					return ( A * exp(t2 * position) ); 
				}
			}
			else {
				// Case E < V
				if (position < 0.0) {
					// incoming sine wave
					return (CC * exp(t1 * position) + DD * exp(-1.0 * t1 * position));
				}
				else {
					// outgoing exponential decay
					return ( B * exp(t2 * position) ); 
				}
			}
		}
		else { 
			std::string reason = "Error: double pot_step::wavefunction(double position)\n"; 
			reason += "No parameters defined for pot_step class\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void pot_step::compute_wavefunction(std::string filename)
{
	// send the computed wavefunction to a file
	// R. Sheehan 20 - 8 - 2021

	try {
		if (params_defined && filename != empty_str) {
			std::ofstream write;

			write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

			if (write.is_open()) {

				int nn = 501; 
				double x0 = -3.0, dx = (2.0*fabs(x0)) / static_cast<double>(nn-1);
				std::complex<double> psi; 

				for (int i = 0; i < nn; i++) {
					//write << std::setprecision(10) << x0 << " , " << abs(wavefunction(x0)) << "\n"; 
					psi = wavefunction(x0); 
					write << std::setprecision(10) << x0 << " , " << psi.real() << " , " << psi.imag() << " , " << template_funcs::DSQR(psi.real()) + template_funcs::DSQR(psi.imag()) << "\n";
					x0 += dx; 
				}
				
				write.close(); 
			}
			else {
				std::string reason = "Error: void pot_step::compute_wavefunction(std::string filename)\n";
				reason += "Could not open file: " + filename + "\n"; 
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason = "Error: void pot_step::compute_wavefunction(std::string filename)\n";
			if(!params_defined) reason += "No parameters defined for pot_step class\n";
			if (filename == empty_str) reason += "Invalid filename\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}
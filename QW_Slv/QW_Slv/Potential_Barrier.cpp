#ifndef ATTACH_H
#include "Attach.h"
#endif

// Class definitions for the potential barrier class
// R. Sheehan 31 - 8 - 2021

pot_barr::pot_barr()
{
	// Default constructor
	m = W = E = V = p1 = p2 = T = R = AA = 0.0;
	
	BB = CC = DD = EE = t1 = t2 = zero; 
}

pot_barr::pot_barr(double particle_mass, double particle_energy, double barr_height, double barr_width)
{
	// Primary constructor
	set_params(particle_mass, particle_energy, barr_height, barr_width);
}

void pot_barr::set_params(double particle_mass, double particle_energy, double barr_height, double barr_width, bool loud)
{
	// assign values to the parameters for the potential barrier calculation
	// particle mass in units of kg
	// energies input in units of eV but converted to J for sake of calculation
	// distances assumed to be in units of nm

	try {
		bool c1 = particle_mass > 0.0 ? true : false;
		bool c2 = particle_energy > 0.0 ? true : false;
		bool c3 = barr_height > particle_energy ? true : false;
		bool c4 = barr_width > 0.0 ? true : false;
		bool c10 = c1 && c2 && c3 && c4;

		if (c10) {
			m = particle_mass;
			W = barr_width; 
			E = template_funcs::convert_ev_J(particle_energy);
			V = template_funcs::convert_ev_J(barr_height);
			p1 = sqrt(2.0 * m * E);
			p2 = sqrt(2.0 * m * (V - E));
			
			AA = 1.0; 
			BB = (0.5 * AA) * (1.0 + eye * (p1 / p2) ) * exp( ( ( eye * p1 ) - p2) * ((W * 1.0e-9) / H_BAR_J) );
			CC = (0.5 * AA) * (1.0 - eye * (p1 / p2) ) * exp( ( ( eye * p1 ) + p2) * ((W * 1.0e-9) / H_BAR_J) );
			DD = (0.5 * CC) * (1.0 + eye * (p2 / p1)) + (0.5 * BB) * (1.0 - eye * (p2 / p1)); 
			EE = (0.5 * CC) * (1.0 - eye * (p2 / p1)) + (0.5 * BB) * (1.0 + eye * (p2 / p1)); 

			t1 = (eye * p1 * 1.0e-9) / H_BAR_J; // momentum outside the barrier, scale length to nm
			t2 = (1.0e-9 * p2) / H_BAR_J; // momentum inside the barrier, scale length to nm

			T = template_funcs::DSQR(abs(AA)) / template_funcs::DSQR(abs(DD)); 
			R = template_funcs::DSQR(abs(EE)) / template_funcs::DSQR(abs(DD)); 

			if (loud) {
				std::cout << "DD = " << DD << " , EE = " << EE << "\nBB = " << BB << " , CC = " << CC << "\nAA = " << AA << "\n";
				std::cout << "p1 = " << p1 << " , p2 = " << p2 << "\n";
				std::cout << "t1 = " << t1 << " , t2 = " << t2 << "\n";
				std::cout << "T = " << T << ", R = " << R << ", T + R = " << T + R << "\n";
				std::cout << "Scale FOM: (p2 a) / hbar = " << (p2 * W * 1.0e-9) / H_BAR_J << "\n"; 
				std::cout << "Scale FOM: (p2 a) / hbar << 1 => |CC| ~ |BB| => Thin Barrier\n";
				std::cout << "Scale FOM: (p2 a) / hbar >> 1 => |CC| >> |BB| => Thick Barrier\n";
			}

			params_defined = true;
		}
		else {
			std::string reason = "Error: void pot_step::set_params(double particle_mass, double particle_energy, double step_height)\n";
			if (!c1) reason += "particle_mass is not positive\n";
			if (!c2) reason += "particle_energy is not positive\n";
			if (!c3) reason += "barr_height is less than particle_energy\n";
			if (!c4) reason += "barr_width is not positive\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

std::complex<double> pot_barr::wavefunction(double position)
{
	// Compute the value of the wavefunction for the potential barrier problem
	// R. Sheehan 31 - 8 - 2021

	try {
		if (params_defined) {
			if (position < 0.0) {
				// wavefunction before barrier
				return ( DD * exp(t1 * position) + EE * exp(-1.0*t1 * position) );
			}
			else if (position > W) {
				// wavefunction after barrier
				return ( AA * exp(t1 * position) );
			}
			else {
				// wavefunction inside barrier
				return ( BB * exp(t2 * position) + CC * exp(-1.0 * t2 * position) );
			}
		}
		else {
			std::string reason = "Error: std::complex<double> pot_barr::wavefunction(double position)\n";
			reason += "No parameters defined for pot_barr class\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void pot_barr::compute_wavefunction(std::string filename)
{
	// send the computed wavefunction to a file
	// R. Sheehan 20 - 8 - 2021

	try {
		if (params_defined && filename != empty_str) {
			std::ofstream write;

			write.open(filename.c_str(), std::ios_base::out | std::ios_base::trunc);

			if (write.is_open()) {

				int nn = 1001;
				double x0 = -3.0, dx = (2.0 * fabs(x0)) / static_cast<double>(nn - 1);
				std::complex<double> psi;

				for (int i = 0; i < nn; i++) {
					psi = wavefunction(x0);
					write << std::setprecision(10) << x0 << " , " << psi.real() << " , " << psi.imag() << " , " << template_funcs::DSQR(psi.real()) + template_funcs::DSQR(psi.imag()) << "\n";
					//write << std::setprecision(10) << x0 << " , " << psi.real() << " , " << psi.imag() << "\n";
					x0 += dx;
				}

				write.close();
			}
			else {
				std::string reason = "Error: void pot_barr::compute_wavefunction(std::string filename)\n";
				reason += "Could not open file: " + filename + "\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason = "Error: void pot_barr::compute_wavefunction(std::string filename)\n";
			if (!params_defined) reason += "No parameters defined for pot_barr class\n";
			if (filename == empty_str) reason += "Invalid filename\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}
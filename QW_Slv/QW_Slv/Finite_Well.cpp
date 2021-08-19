#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the methods associated with the finite well class
// R. Sheehan 24 - 10 - 2016

fin_well::fin_well()
{
	// Default Constructor

	n_states = 0; 

	M_ratio = 1.0; 
	A = L = Lhalf = M_well = M_barr = x_c = E_range = 0.0; 
}

fin_well::fin_well(double length, double mass_well, double mass_barrier, double barrier_height, double centre_position)
{
	// Constructor

	// length is the length of the well expressed in nm
	// mass is the mass of the particle in the well expressed as a multiple of the electron mass, which is in kg
	// barrier_height is the depth of the quantum well expressed in eV
	// centre_position of the well determines where the centre of the well is located and is expressed in nm

	set_well_params(length, mass_well, mass_barrier, barrier_height, centre_position); 
}

void fin_well::set_well_params(double length, double mass_well, double mass_barrier, double barrier_height, double centre_position)
{
	// method for assigning values to the well parameters

	// length is the length of the well expressed in nm
	// mass is the mass of the particle in the well expressed as a multiple of the electron mass, which is in kg
	// barrier_height is the depth of the quantum well expressed in eV
	// centre_position of the well determines where the centre of the well is located and is expressed in nm

	try{
	
		bool c1 = length > 0.0 ? true : false;
		bool c2 = mass_well > 0.0 ? true : false;
		bool c2a = mass_barrier > 0.0 ? true : false;
		bool c3 = barrier_height > 0.0 ? true : false;

		if(c1 && c2 && c3){

			L = length; 
			Lhalf = 0.5*L; 
			A = sqrt(2.0/L); // normalistion constant for energy eigenfunctions, it can be any complex number whose absolute value is sqrt(2.0/L)
			M_well = mass_well; 
			M_barr = mass_barrier;
			M_ratio = (M_barr / M_well); 
			x_c = centre_position; 
			well_depth = barrier_height; 
			E_range = 2.0 * M_well * well_depth * template_funcs::DSQR(Lhalf/H_BAR_eV); // range over which energies are sought
		}
		else{
			std::string reason = "Error: void fin_well::set_well_params(double length, double mass, double barrier_height, double centre_position)\n"; 
			if(!c1) reason += "length is not positive\n"; 
			if(!c2) reason += "mass in well is not positive\n"; 
			if(!c2a) reason += "mass in barrier is not positive\n"; 
			if(!c3) reason += "barrier height is not positive\n"; 
			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double fin_well::energy_eigenvalue(int n)
{
	// return the n^{th} energy eigenvalue in units of eV

	try{
	
		if(n > -1 && n < n_states){
			return energy_levels[n];
		}
		else{
			std::string reason = "Error: double fin_well::energy_eigenvalue(int n)\n"; 
			reason += "Value of n must be in range of allowed values\n"; 
			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what();
		return 0.0; 
	}
}

double fin_well::K(double beta)
{
	// wavenumber in the well

	try{
		
		return ( sqrt(2.0*M_well*beta)/H_BAR_eV ); 

	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what();
		return 0.0; 
	}
}

double fin_well::Alpha(double beta)
{
	// Decay constant in barrier

	try{
		return ( sqrt( 2.0*M_barr*(well_depth - beta) )/H_BAR_eV ); 
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what();
		return 0.0; 
	}
}

double fin_well::energy_eigenequation(double beta, bool state)
{
	// return the value of the energy eigenequation

	try{
	
		if(state){
			return even_eigenequation(beta); 
		}
		else{
			return odd_eigenequation(beta); 
		}

	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what();
		return 0.0; 
	}
}

double fin_well::even_eigenequation(double beta)
{
	// method that returns the value of the eigenequation for the even states

	try{
		double kval = K(beta); 

		return ( Alpha(beta) - M_ratio*kval*tan(kval*Lhalf) ); 
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what();
		return 0.0; 
	}
}

double fin_well::odd_eigenequation(double beta)
{
	// method that returns the value of the eigenequation for the odd states
	// Note cot(x) = tan( (pi/2) - x )

	try{
		double kval = K(beta); 

		return ( Alpha(beta) + M_ratio*kval*tan(PI_2 - kval*Lhalf) ); 
	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what();
		return 0.0; 
	}
}
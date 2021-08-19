#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definition of the methods declared in Infinite_Well.h

inf_well::inf_well()
{
	// Default constructor
	A = L = Lhalf = M = x_c = En_const = kn_const = bndry = 0.0; 
}

inf_well::inf_well(double length, double mass, double centre_position)
{
	// Constructor

	set_well_params(length, mass, centre_position); 
}

void inf_well::set_well_params(double length, double mass, double centre_position)
{
	try{
	
		bool c1 = length > 0.0 ? true : false;
		bool c2 = mass > 0.0 ? true : false;

		if(c1 && c2){

			L = length; 
			Lhalf = 0.5*L; 
			A = sqrt(2.0/L); // normalistion constant for energy eigenfunctions, it can be any complex number whose absolute value is sqrt(2.0/L)
			M = mass; 
			x_c = centre_position; 
			bndry = x_c + Lhalf; 
			En_const = template_funcs::DSQR(PLANCK_CONST_eV) / (8.0 * M * template_funcs::DSQR(L) ); // h^{2} / 8 m L^{2}
			kn_const = PI / L; // \pi / L
		}
		else{
			std::string reason = "Error: void inf_well::set_well_params(double length, double mass, double centre_position)\n"; 
			if(!c1) "length is not positive\n"; 
			if(!c2) "mass is not positive\n"; 
			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

double inf_well::energy_eigenvalue(int n)
{
	// return the n^{th} energy eigenvalue in units of eV

	try{
	
		if(n > 0){
			return (template_funcs::DSQR(n) * En_const); // n^{2} h^{2} / 8 m L^{2}
		}
		else{
			std::string reason = "Error: double inf_well::energy_eigenvalue(int n)\n"; 
			reason += "Value of n must be greater than zero\n"; 
			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what();
		return 0.0; 
	}
}

double inf_well::energy_eigenfunction(int n, double position)
{
	// compute the value of the n^{th} energy eigenfunction
	// position length scale is in nano-metres

	try{
	
		bool c1 = n > 0 ? true : false; 
		bool c2 = fabs(position) < bndry ? true : false; 

		if(c1 && c2){
			double arg = n * kn_const * (position - bndry); // ( n \pi / L ) * ( x - x_{c} + L/2 )
			return A * sin( arg ); 
		}
		else{
			std::string reason = "Error: double inf_well::energy_eigenfunction(int n, double position)\n"; 
			if(!c1) reason += "Value of n must be greater than zero\n"; 
			if(!c2) reason += "Position must be located inside the well\n"; 
			throw std::invalid_argument(reason); 
		}

	}
	catch(std::invalid_argument &e){
		std::cerr<<e.what();
		return 0.0; 
	}
}
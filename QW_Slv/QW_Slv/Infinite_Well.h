#ifndef INFINITE_WELL_H
#define INFINITE_WELL_H

// Implementation of a class that computes the solution of an infinite square well
// For calculation details see Bransden and Joachain.
// R. Sheehan 24 - 10 - 2016

class inf_well{

public: 
	inf_well(); 

	inf_well(double length, double mass, double centre_position = 0.0); 

	void set_well_params(double length, double mass, double centre_position = 0.0); 

	double energy_eigenvalue(int n); // return the energy associated with the n^{th} energy level

	double energy_eigenfunction(int n, double position); // return the value of the normalised wavefunction at some position inside the well

private:
	double L; // length of the well
	double Lhalf; // half the length of the well
	double M; // mass of the particle in the well
	double A; // wavefunction normalisation constant
	double x_c; // position of the centre of the well
	double En_const; // Constant associated with energy eigenvalues
	double kn_const; // Constant associated with energy eigenfunctions
	double bndry; // position of the well boundary
}; 

#endif
#ifndef FINITE_WELL_H
#define FINITE_WELL_H

// Implementation of the code required to compute the solutions in a finite quantum well
// R. Sheehan 24 - 10 - 2016

// The natural scale for energy is eV, scale energies accordingly
// The natural scale for length is nm, you can scale lengths accordingly

class fin_well{
public:
	fin_well(); 

	fin_well(double length, double mass_well, double mass_barrier, double barrier_height, double centre_position = 0.0); 

	void set_well_params(double length, double mass_well, double mass_barrier, double barrier_height, double centre_position = 0.0); 

	double energy_eigenvalue(int n); // return the energy associated with the n^{th} energy level

	double energy_eigenfunction(int n, double position); // return the value of the normalised wavefunction at some position inside the well

private:
	double energy_eigenequation(double beta, bool state);

	double even_eigenequation(double beta);

	double odd_eigenequation(double beta); 

	double K(double beta); 

	double Alpha(double beta); 

	void solve_energy_eigenequation(); 

private:
	int n_states; // num. of bound states in the well

	double L; // length of the well expressed in nm
	double Lhalf; // half the length of the well expressed in nm
	double M_well; // mass of the particle in the well expressed as a multiple of electron mass, which is in kg
	double M_barr; // mass of the particle in the well expressed as a multiple of electron mass, which is in kg
	double M_ratio; // ratio of mass in barrier to mass in well
	double A; // wavefunction normalisation constant
	double x_c; // position of the centre of the well expressed in nm
	double well_depth; // height of energy barrier expressed in eV
	double E_range; // range over which energies are sought expressed in eV

	std::vector<double> energy_levels;  
}; 

#endif
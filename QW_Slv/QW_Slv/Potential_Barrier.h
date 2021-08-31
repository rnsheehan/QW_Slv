#ifndef POTENTIAL_BARRIER_H
#define POTENTIAL_BARRIER_H

// Code for computing the solution in a potential barrier
// Notation taken from "Quantum Theory" by David Bohm
// R. Sheehan 31 - 8 - 2021

class pot_barr {
public:
	pot_barr(); 
	pot_barr(double particle_mass, double particle_energy, double barr_height, double barr_width);

	void set_params(double particle_mass, double particle_energy, double barr_height, double barr_width, bool loud = false);

	std::complex<double> wavefunction(double position);

	void compute_wavefunction(std::string filename);

	// getters
	inline double get_E() { return E; }
	inline double get_V() { return V; }
	inline double get_T() { return T; }
	inline double get_R() { return R; }

private:
	bool params_defined; // boolean to decide if parameters have been assigned to the class
	double m; // particle mass in units of kg
	double W; // barrier width in units of nm
	double E; // particle energy in units of J
	double V; // step height in units of J
	double p1; // momentum before / after barrier
	double p2; // momentum within barrier
	double T; // transmission probability
	double R; // reflection probability

	double AA; // solution constant
	std::complex<double> BB; // solution constant
	std::complex<double> CC; // solution constant
	std::complex<double> DD; // solution constant
	std::complex<double> EE; // solution constant
	std::complex<double> t1; // solution constant
	std::complex<double> t2; // solution constant
};

#endif

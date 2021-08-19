#ifndef POTENTIAL_STEP_H
#define POTENTIAL_STEP_H

// Implementation of code required to compute solution for potential step
// R. Sheehan 19 - 8 - 2021

class pot_step {
public:
	pot_step();
	pot_step(double particle_mass, double particle_energy, double step_height);

	void set_params(double particle_mass, double particle_energy, double step_height); 

private:
	double m; // particle mass in units of kg
	double E; // particle energy in units of eV
	double V; // step height in units of eV
	double p1; // momentum before step
	double p2; // momentim after step
	double T; // transmission probability
	double R; // reflection probability
	double A; // solution constant
	double B; // solution constant
	double C; // solution constant
	double D; // solution constant
	std::complex<double> CC; // solution constant
	std::complex<double> DD; // solution constant
};

#endif
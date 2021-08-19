#ifndef CHEBYSHEV_APPROXIMATION_H
#define CHEBYSHEV_APPROXIMATION_H

// Implementation of a namespace for methods that implement the Chebyshev polynomial approximation to an arbitrary function
// This will be used to implement the Bessel functions of non-integer order and the Airy functions
// R. Sheehan 17 - 10 - 2016

namespace cheb_appr{
	
	void chebft(double a, double b, double c[], int n, double (*func)(double)); 

	double chebev(double a, double b, double c[], int m, double x); 

	void chder(double a, double b, double c[], double cder[], int n); 

	void chint(double a, double b, double c[], double cint[], int n); 

}

#endif
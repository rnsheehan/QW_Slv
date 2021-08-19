#ifndef ATTACH_H
#include "Attach.h"
#endif

void cheb_appr::chebft(double a, double b, double c[], int n, double (*func)(double))
{
	// Chebyshev fit to a function func(x) on the interval [a, b]
	// Given a function func, lower and upper limits of the interval [a,b] and a maximum degree n
	// this routine computes the n coefficients c[0..n-1] such that func(c) \approx -c_{0}/2 + \sum_{k=0}^{n-1}c_{k}T_{k}(y)
	// where y = \frac{x-0.5*(b+a)}{0.5*(b-a)}
	// Routine is to be used with moderately large n, 30 or 50, the array is then subsequently truncated when used in the function evaluation
	// by choosing m < n such that c_{m} and subsequent elements are negligible

	try{

		bool c1 = ( a < b ? true : false); 
		bool c2 = ( fabs(b - a) > EPS ? true : false);
		bool c3 = ( n > 29 ? true : false);

		if(c1 && c2 && c3){

			int k,j;
			double fac,bpa,bma,*f;

			//f=vector(0,n-1);
			f = new (double [n]); 
			bma=0.5*(b-a);
			bpa=0.5*(b+a);
			for (k=0;k<n;k++) {
				double y=cos(PI*(k+0.5)/n);
				f[k]=(*func)(y*bma+bpa);
			}
			fac=2.0/n;
			for (j=0;j<n;j++) {
				double sum=0.0;
				for (k=0;k<n;k++)
					sum += f[k]*cos(PI*j*(k+0.5)/n);
				c[j]=fac*sum;
			}
			//free_vector(f,0,n-1);
			delete[] f; 
		}
		else{
			std::string reason = "Error: void cheb_appr::chebft(double a, double b, double c[], int n, double (*func)(double))\n"; 
			if(!c1 || !c2) reason += "Domain endpoints are not correctly defined\na = " + template_funcs::toString(a, 3) + ", b = "+ template_funcs::toString(b, 3) + "\n"; 
			if(!c3) reason += "Number of expansion coefficients n = " + template_funcs::toString(n) + " is not sufficient\n"; 
			throw std::invalid_argument(reason); 
		}	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
	catch(std::bad_alloc &ba){
		std::string reason = "Error: void cheb_appr::chebft(double a, double b, double c[], int n, double (*func)(double))\n";
		reason += ba.what(); 
		useful_funcs::exit_failure_output(reason); 
		exit(EXIT_FAILURE); 
	}
}

double cheb_appr::chebev(double a, double b, double c[], int m, double x)
{
	// Evaluation of Chebyshev approximation to func(x) on [a, b]
	// c[0..n-1] is the array of expansion coefficients computed using chebft
	// The Chebyshev polynomial -c_{0}/2 + \sum_{k=0}^{n-1}c_{k}T_{k}(y) is evaluated at the point y = \frac{x-0.5*(b+a)}{0.5*(b-a)}
	// and the result is returned as a function value

	try{

		bool c1 = ( a < b ? true : false); 
		bool c2 = ( fabs(b - a) > EPS ? true : false);
		bool c3 = ( m > 1 ? true : false);
		bool c4 = ((x-a)*(x-b) < 0.0 ? true : false); 

		if(c1 && c2 && c3 & c4){

			double d=0.0,dd=0.0,sv,y,y2;
			int j;

			//if ((x-a)*(x-b) > 0.0) std::cerr<<"x not in range in routine chebev\n";

			y2=2.0*(y=(2.0*x-a-b)/(b-a));
			for (j=m-1;j>=1;j--) {
				sv=d;
				d=y2*d-dd+c[j];
				dd=sv;
			}
			return y*d-dd+0.5*c[0];
	
		}
		else{
			std::string reason = "Error: double cheb_appr::chebev(double a, double b, double c[], int m, double x)\n"; 
			if(!c1 || !c2) reason += "Domain endpoints are not correctly defined\na = " + template_funcs::toString(a, 3) + ", b = "+ template_funcs::toString(b, 3) + "\n"; 
			if(!c3) reason += "Number of expansion coefficients m = " + template_funcs::toString(m) + " is not sufficient\n"; 
			if(!c4) reason += "x not in range in routine chebev\n"; 
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void cheb_appr::chder(double a, double b, double c[], double cder[], int n)
{
	// Chebyshev polynomial approximation to the derivative of a function func(x) on [a, b]
	// Given a, b, c[0..n-1] as output from routine chebft, compute cder[0..n-1], the 
	// Chebyshev expansion coefficients of the derivative of the function func(x) on [a, b]

	try{

		bool c1 = ( a < b ? true : false); 
		bool c2 = ( fabs(b - a) > EPS ? true : false);
		bool c3 = ( n > 29 ? true : false);

		if(c1 && c2 && c3){

			int j;
			double con;

			cder[n-1]=0.0;
			cder[n-2]=2*(n-1)*c[n-1];
			for (j=n-3;j>=0;j--)
				cder[j]=cder[j+2]+2*(j+1)*c[j+1];
			con=2.0/(b-a);
			for (j=0;j<n;j++)
				cder[j] *= con;
	
		}
		else{
			std::string reason = "Error: void cheb_appr::chebft(double a, double b, double c[], int n, double (*func)(double))\n"; 
			if(!c1 || !c2) reason += "Domain endpoints are not correctly defined\na = " + template_funcs::toString(a, 3) + ", b = "+ template_funcs::toString(b, 3) + "\n"; 
			if(!c3) reason += "Number of expansion coefficients n = " + template_funcs::toString(n) + " is not sufficient\n"; 
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}

void cheb_appr::chint(double a, double b, double c[], double cint[], int n)
{
	// Chebyshev polynomial approximation to the integral of a function func(x) on [a, b]
	// Given a, b, c[0..n-1] as output from routine chebft, compute cint[0..n-1], the 
	// Chebyshev expansion coefficients of the integral of the function func(x) on [a, b]
	// The constant of integration is set so that the integral vanishes at point a

	try{

		bool c1 = ( a < b ? true : false); 
		bool c2 = ( fabs(b - a) > EPS ? true : false);
		bool c3 = ( n > 29 ? true : false);

		if(c1 && c2 && c3){

			int j;
			double sum=0.0,fac=1.0,con;

			con=0.25*(b-a);
			for (j=1;j<=n-2;j++) {
				cint[j]=con*(c[j-1]-c[j+1])/j;
				sum += fac*cint[j];
				fac = -fac;
			}
			cint[n-1]=con*c[n-2]/(n-1);
			sum += fac*cint[n-1];
			cint[0]=2.0*sum;
	
		}
		else{
			std::string reason = "Error: void cheb_appr::chebft(double a, double b, double c[], int n, double (*func)(double))\n"; 
			if(!c1 || !c2) reason += "Domain endpoints are not correctly defined\na = " + template_funcs::toString(a, 3) + ", b = "+ template_funcs::toString(b, 3) + "\n"; 
			if(!c3) reason += "Number of expansion coefficients n = " + template_funcs::toString(n) + " is not sufficient\n"; 
			throw std::invalid_argument(reason); 
		}
	
	}
	catch(std::invalid_argument &e){
		useful_funcs::exit_failure_output(e.what()); 
		exit(EXIT_FAILURE); 
	}
}
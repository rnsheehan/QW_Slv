#ifndef ATTACH_H
#define ATTACH_H

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <string>
#include <sstream>
#include <fstream>

// need these for directory manipulation
#include <direct.h>
#include <errno.h>

#include <cmath>
#include <complex> 
#include <vector>

// Constants
static const double EPS=(1.0e-16);
static const double FPMIN=(1.0e-30);

static const double p=(atan(1.0)); // pi / 4
static const double Two_PI=(8.0*p); // 2 pi
static const double PI=(4.0*p); // pi
static const double PI_2=(2.0*p); // pi / 2
static const double PI_3=((4.0/3.0)*p); // pi / 3
static const double PI_4=(p); // pi / 4
static const double PI_5=((4.0/5.0)*p); // pi / 5
static const double PI_6=((2.0/3.0)*p); // pi / 6

static const double M_ELECTRON_KG = 9.109384e-31; // mass of electron in kg
static const double PLANCK_CONST_J = 6.6262e-34; // Planck's constant in Js
static const double PLANCK_CONST_eV = 4.1357e-15; // Planck's constant in eVs
static const double H_BAR_J = 1.05457e-34; // hbar = Planck's constant over 2 pi in Js
static const double H_BAR_eV = 6.5822e-16; // hbar = Planck's constant over 2 pi in eVs
//static const double One_H_BAR = 9.482534e+33; // 1/hbar

static const std::complex<double> zero(0.0, 0.0);
static const std::complex<double> eye(0.0, 1.0);
static const std::complex<double> one(1.0, 0.0);

#include "Templates.h"
#include "Useful.h"

#include "Potential_Step.h"
#include "Infinite_Well.h"
#include "Finite_Well.h"

#include "Test_Routines.h"
//#include "Chebyshev_Approximation.h"
//#include "Special_Functions.h"

#endif
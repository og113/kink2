/*
	declarations for function pointers
	
	N.B. this can only be included by the .cc file with main()
*/

#ifndef __FNPTRS_H_INCLUDED__
#define __FNPTRS_H_INCLUDED__

#include <complex>
#include "potentials.h"

using namespace std;

typedef complex<double> comp;

/*-------------------------------------------------------------------------------------------------------------------------
potential functions and their derivatives
	- (*V)
	- (*Vd)
	- (*dV)
	- (*dVd)
	- (*ddV)
	- (*ddVd)
	
-------------------------------------------------------------------------------------------------------------------------*/

//function pointer V
comp (*V) (const comp& phi, const params_for_V& parameters); 

//function pointer Vd, for gsl
double (*Vd) (double phi,  params_for_V parameters);

//function pointer dV
comp (*dV) (const comp& phi, const params_for_V& parameters);

//function pointer dVd, for gsl
double (*dVd) (double phi,  params_for_V parameters);

////function pointer ddV
comp (*ddV) (const comp& phi, const params_for_V& parameters);

//function pointer ddVd, for gsl
double (*ddVd) (double phi,  params_for_V parameters);

#endif // __FNPTRS_H_INCLUDED__

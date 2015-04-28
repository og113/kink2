/*-------------------------------------------------------------------------------------------------------------------------
 	declarations for functions to compute sphaleron
 -------------------------------------------------------------------------------------------------------------------------*/
 
#ifndef __SPHALERON_FNS_H_INCLUDED__
#define __SPHALERON_FNS_H_INCLUDED__

using namespace std;

typedef unsigned int uint;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. ode functions
	2. pde functions
	3. energy integrand
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. ode functions
		- void params
		- func
		- jac
-------------------------------------------------------------------------------------------------------------------------*/

// void params
struct void_params {};

// func
int func (double t, const double y[], double f[], void *params);

// jac
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

/*-------------------------------------------------------------------------------------------------------------------------
	2. pde functions
		- force parameters
		- funcPDE
		- jacPDE
-------------------------------------------------------------------------------------------------------------------------*/
	
// force parameters
struct force_params{double x0; double x1; unsigned int Nx; double sigma;};

// a function to give 'force' for time evolution
int funcPDE (double t, const double y[], double f[], void *params);
	
// jacobian for PDE solver	
int jacPDE (double t, const double y[], double *dfdy, double dfdt[], void *params);

/*-------------------------------------------------------------------------------------------------------------------------
	3. energy integrand
		- E_params
		- E_integrand
-------------------------------------------------------------------------------------------------------------------------*/

struct E_params {double Y_0;};
	
// energy integrand
double E_integrand (double x, void * parameters);

#endif // __SPHALERON_FNS_H_INCLUDED__

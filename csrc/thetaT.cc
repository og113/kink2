/*-------------------------------------------------------------------------------------------------------------------------
	definitions for functions specifically useful to the theta/T method
-------------------------------------------------------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include "gsl_extras.h"
#include "potentials.h"
#include "thetaT.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - collection of specific functions
	2 - calculation of secondary parameters
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. specific functions F, FDF etc
		- wrapped functions
		- VdV
		- dVddV
		- struct ec_params
		- ec (energy change: V(minima[1])-V(minima[0])-dE)
		- S1 integrand
		- rho integrand
-------------------------------------------------------------------------------------------------------------------------*/

//function pointer Vd, for gsl - "static" for functions just makes them local to the file (n.b. "static" is different for variables)
static double (*Vd_local) (const double& phi, const params_for_V& parameters);
static double (*dVd_local) (const double& phi, const params_for_V& parameters);
static double (*ddVd_local) (const double& phi, const params_for_V& parameters);

// wrapped functions
#define WRAP1(FN) double FN##_wrapped1(double x, void* parameters) { return FN(x, *((params_for_V*)parameters)); }
WRAP1(Vd_local)
WRAP1(dVd_local)
WRAP1(ddVd_local)
#undef WRAP1

//V FDF gsl function
void VdV (double x, void * parameters, double * f, double* df) 
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	*f =  Vd_local_wrapped1(x,params);
	*df = dVd_local_wrapped1(x,params);
	}
	
//dV FDF gsl functions
void dVddV (double x, void * parameters, double * f, double* df) 
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	*f =  dVd_local_wrapped1(x,params);
	*df = ddVd_local_wrapped1(x,params);
	}

//energy change gsl function : V(minima[1])-V(minima[0])-dE
double ec (double epsi, void * parameters) {
	struct ec_params * paramsIn = (struct ec_params *)parameters;
	struct params_for_V paramsOut;
	paramsOut.epsi = epsi;
	paramsOut.aa = (paramsIn->aa);
	double minima1 = (paramsIn->minima1);
	double minima0 = (paramsIn->minima0);
	double de = (paramsIn->de);
	return Vd_local(minima0,paramsOut) - Vd_local(minima1,paramsOut) - de;
}
	
//S1 integrand
double s1Integrand (double x, void * parameters) {
	struct params_for_V * params = (struct params_for_V *)parameters;
	return pow(2.0*Vd_local(x,*params),0.5);
}

//rho integrand
double rhoIntegrand (double x, void * parameters) {
	struct params_for_V * params = (struct params_for_V *)parameters;
	return pow(2.0*Vd_local(x,*params),-0.5);
}
	
/*-------------------------------------------------------------------------------------------------------------------------
	2. calculation of secondary parameters
		- epsilonFn
-------------------------------------------------------------------------------------------------------------------------*/

//program to find epsilon given gsl functions df and dE
void epsilonFn (gsl_function * xF, gsl_function * xEC, const double * xdE, double * xEpsilon, vector<double>* xMinima)
	{
	double closenessdE = 1.0e-14;
	vector<double> dE_test(1);	dE_test[0] = 1.0;
	double newdE = *xdE;
	struct params_for_V * Fparameters = (struct params_for_V *) (*xF).params;
	struct ec_params * ECparameters = (struct ec_params *) (*xEC).params;
	unsigned int counter = 0;
	unsigned int maxCounter = 1e4;
	while (dE_test.back()>closenessdE)
		{
		//find roots of ec(epsilon)=0
		*xEpsilon = brentRootFinder(xEC,*xEpsilon,*xEpsilon/2.0,*xEpsilon*2.0);
		//assign new value of epsilon to xF
		(*Fparameters).epsi = *xEpsilon;
		(*xF).params = Fparameters;
		//finding new roots of dV(phi)=0
		(*xMinima)[0] = brentMinimum(xF,-1.0,-3.0,0.0);
		(*xMinima)[1] = brentMinimum(xF,1.2,0.5,3.0);
		//assign new roots to xECDF
		(*ECparameters).minima0 = (*xMinima)[0];
		(*ECparameters).minima1 = (*xMinima)[1];
		(*xEC).params = ECparameters;
		//evaluating new dE
		newdE = (*(*xEC).function)(*xEpsilon,ECparameters) + *xdE;
		//evaluating test
		if (abs(*xdE)>1.0e-16) dE_test.push_back(abs((newdE-(*xdE))/(*xdE)));
		else 						dE_test.push_back(abs(newdE-(*xdE)));
		counter++;
		//test if too many runs
		if (counter>maxCounter)
			{
			cout << "epsilonFn error, more that " << maxCounter << " loops, consider reducing closenessdE" << endl;
			cout << "dE_test.back() = " << dE_test.back() << " , closenessdE = " << closenessdE << endl;
			cout << "dE = " << *xdE << " , minima[0] = " << (*xMinima)[0] << " , minima[1] = " << (*xMinima)[1];
			cout << " , epsilon = " << *xEpsilon << endl << endl;
			break;
			}
		}
	//*xdE = newdE;
	}


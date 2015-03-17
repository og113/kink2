/* -------------------------------------------------------------------------------------------------------------------------
	sphaleron
		- functions for program to solve boundary value ode to find sphaleron
------------------------------------------------------------------------------------------------------------------------- */

#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "sphaleron.h"

using namespace std;

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
		- func
		- jac
-------------------------------------------------------------------------------------------------------------------------*/

// func
int func (double t, const double y[], double f[], void *params) {
	f[0] = y[1];
	f[1] = -2.0*y[1]/t + y[0] - pow(y[0],3.0);
	f[2] = y[3];
	f[3] = -2.0*y[3]/t + y[2] - 3.0*y[2]*pow(y[0],2.0);
	return GSL_SUCCESS;
}

// jac
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params) {
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 4, 4);
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set_all (m, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
	gsl_matrix_set (m, 1, 0, 1.0 - 3.0*pow(y[0],2.0));
	gsl_matrix_set (m, 1, 1, -2.0/t);
	gsl_matrix_set (m, 2, 3, 1.0);
	gsl_matrix_set (m, 3, 0, -6.0*y[2]*y[0]);
	gsl_matrix_set (m, 3, 2, 1.0 - 3.0*pow(y[0],2.0));
	gsl_matrix_set (m, 3, 3, -2.0/t);
	dfdt[0] = 0.0;
	dfdt[1] = 2.0*y[1]/pow(t,2.0);
	dfdt[2] = 0.0;
	dfdt[3] = 2.0*y[3]/pow(t,2.0);
	return GSL_SUCCESS;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. pde functions
		- funcPDE
		- jacPDE
-------------------------------------------------------------------------------------------------------------------------*/

// a function to give 'force' for time evolution
int funcPDE (double t, const double y[], double f[], void *params) {
	struct force_params * parameters = (struct force_params *)params;
	double x0 = (parameters->x0);
	double x1 = (parameters->x1);
	uint Nx = (parameters->Nx);
	double sigma = (parameters->sigma);
	double dx = (x1-x0), x;
	dx /= (double)Nx;
	uint yLength = (uint)(sizeof(y)/sizeof(*y));
	if (yLength!=2*Nx){printf("funcPDE error: yLength = %i",yLength);}
	f[0] = y[1];
	f[1] = 2.0*(y[2]-y[0])/pow(dx,2.0) - y[0] + pow(y[0],3.0);
	f[1] *= sigma;
	for (uint j=1; j<(Nx-1); j++) {
		x = x0 + (double)j*dx;
		f[2*j] = y[2*j+1];
		f[2*j+1] = (y[2*(j+1)] + y[2*(j-1)]-2.0*y[2*j])/pow(dx,2.0) + (y[2*(j+1)] - y[2*(j-1)])/x/dx - y[2*j] + pow(y[2*j],3.0);
		f[2*j+1] *= sigma;
	}
	x = x0 + (double)(Nx-1.0)*dx;
	f[2*(Nx-1)] = y[2*(Nx-1)+1];
	f[2*(Nx-1)+1] = (y[2*(Nx-2)]-2.0*y[2*(Nx-1)])/pow(dx,2.0) + (- y[2*(Nx-2)])/x/dx - y[2*(Nx-1)] + pow(y[2*(Nx-1)],3.0);
	f[2*(Nx-1)+1] *= sigma;
	return GSL_SUCCESS;
}
	
// jacobian for PDE solver	
int jacPDE (double t, const double y[], double *dfdy, double dfdt[], void *params) {
	struct force_params * parameters = (struct force_params *)params;
	double x0 = (parameters->x0);
	double x1 = (parameters->x1);
	uint Nx = (parameters->Nx);
	double sigma = (parameters->sigma);
	double dx = (x1-x0), x;
	dx /= (double)Nx;
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2*Nx, 2*Nx);
	gsl_matrix * m = &dfdy_mat.matrix;
	gsl_matrix_set_all (m, 0.0);
	for (uint j=0; j<Nx; j++) {
		if (j==0) {
			gsl_matrix_set(m,2*j,2*j+1, 1.0 );
			gsl_matrix_set(m,2*j+1,2*j, sigma*( -2.0/pow(dx,2.0) - 1.0 + 3.0*pow(y[2*j],2.0) ) );
			gsl_matrix_set(m,2*j+1,2*(j+1), sigma*( 2.0/pow(dx,2.0) ) );
			dfdt[2*j] = 0.0;
			dfdt[2*j+1] = 0.0;
		}
		else if (j==(Nx-1)) {
			x = x0 + (double)(Nx-1.0)*dx;
			gsl_matrix_set(m,2*j,2*j+1, 1.0 );
			gsl_matrix_set(m,2*j+1,2*j, sigma*( -2.0/pow(dx,2.0) - 1.0 + 3.0*pow(y[2*j],2.0) ) );
			gsl_matrix_set(m,2*j+1,2*(j-1), sigma*( 1.0/pow(dx,2.0) - 1.0/x/dx ) );
			dfdt[2*j] = 0.0;
			dfdt[2*j+1] = 0.0;
		}
		else {
			x = x0 + (double)j*dx;
			gsl_matrix_set(m,2*j,2*j+1, 1.0 );
			gsl_matrix_set(m,2*j+1,2*j, sigma*( -2.0/pow(dx,2.0) - 1.0 + 3.0*pow(y[2*j],2.0) ) );
			gsl_matrix_set(m,2*j+1,2*(j+1), sigma*( 1.0/pow(dx,2.0) + 1.0/x/dx ) );
			gsl_matrix_set(m,2*j+1,2*(j-1), sigma*( 1.0/pow(dx,2.0) - 1.0/x/dx ) );
			dfdt[2*j] = 0.0;
			dfdt[2*j+1] = 0.0;
		}
	}
	return GSL_SUCCESS;
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. energy integrand
		- E_integrand
-------------------------------------------------------------------------------------------------------------------------*/

// energy integrand
double E_integrand (double x, void * parameters) {
	if (x<0.0) {
		printf ("error in EIntegrand, R<0\n");
		return 0.0;
	}
	else {
		struct paramsStruct * params = (struct paramsStruct *)parameters;
		double Y_0 = (params->Y_0);
		double y_R[4] = { Y_0, 0.0, 1.0, 0.0};
		int status;
		double t = 1.0e-16;
		gsl_odeiv2_system syse = {func, jac, 4, &paramsVoid};
		gsl_odeiv2_driver * de = gsl_odeiv2_driver_alloc_yp_new (&syse,  gsl_odeiv2_step_rk8pd, 1.0e-9, 1.0e-9, 0.0);
		status = gsl_odeiv2_driver_apply (de, &t, x, y_R);
		if (status != GSL_SUCCESS) {
			printf ("error in EIntegrand, return value=%d\n", status);
			gsl_odeiv2_driver_free (de);
			return (double)status;
		}
		else {
			double to_return = 4.0*pi*pow(x,2.0)*( 0.5*pow(y_R[1],2.0) + 0.5*pow(y_R[0],2.0) - 0.25*pow(y_R[0],4.0) );
			gsl_odeiv2_driver_free (de);
			return to_return;
		}
	}
}

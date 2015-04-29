/*----------------------------------------------------------------------------------------------------------------------------
	staticShooting
		program to solve boundary value problem for static soliton
		using the shooting method
----------------------------------------------------------------------------------------------------------------------------*/
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <limits>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "main.h"

//#define NDEBUG //NDEBUG is to remove error and bounds checking on vectors in SparseLU, for speed - only include once everything works

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		A - ode functions
		B - mass integrand
	
		1 - loading options, closenesses, argv inputs
		2 - getting inputs, checks etc.
		3 - assigning potential fucntions
		4 - defining quantites
		5 - shooting method
		6 - calculating mass
		7 - printing results
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	A. ode functions
		- Potential V, dV, ddV
		- void params
		- func
		- jac
		
	n.b. V, dV and ddV are assigned at runtime in the body of main
-------------------------------------------------------------------------------------------------------------------------*/

// Potential V, dV, ddV
Potential<double> V, dV, ddV;

// void params
struct void_params {};

// func
int func (double t, const double y[], double f[], void *params) {
	f[0] = y[1];
	f[1] = dV(y[0]);
	return GSL_SUCCESS;
}

// jac
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params) {
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set_all (m, 0.0);
	gsl_matrix_set (m, 0, 1, 1.0);
	gsl_matrix_set (m, 1, 0, ddV(y[0]));
	dfdt[0] = 0.0; // no explicit dependence on the variable t
	dfdt[1] = 0.0; 
	return GSL_SUCCESS;
}

/*-------------------------------------------------------------------------------------------------------------------------
	B. mass integrand
		- Mass_params
		- Mass_integrand
-------------------------------------------------------------------------------------------------------------------------*/

struct Mass_params {
	double Y_0;
	double L;
};
	
// Mass integrand
double Mass_integrand (double x, void * parameters) {
	struct Mass_params * params = (struct Mass_params *)parameters;
	double Y_0 = (params->Y_0);
	double L = (params->L);
	double y_R[2] = { Y_0, 0.0};
	int status;
	double t = -L/2.0;
	void_params paramsVoid;
	gsl_odeiv2_system syse = {func, jac, 2, &paramsVoid};
	gsl_odeiv2_driver * de = gsl_odeiv2_driver_alloc_yp_new (&syse,  gsl_odeiv2_step_rk8pd, 1.0e-9, 1.0e-9, 0.0);
	status = gsl_odeiv2_driver_apply (de, &t, x, y_R);
	if (status != GSL_SUCCESS) {
		printf ("error in MassIntegrand, return value=%d\n", status);
		gsl_odeiv2_driver_free (de);
		return (double)status;
	}
	else {
		double to_return = pow(y_R[1],2.0);
		gsl_odeiv2_driver_free (de);
		return to_return;
	}
}


int main(int argc, char** argv)
{
/*----------------------------------------------------------------------------------------------------------------------------
	1. loading options, argv inputs
		- loading options
		- argv inputs
		- defining timenumber
----------------------------------------------------------------------------------------------------------------------------*/

// loading opts
Options opts;
opts.load("optionsS");
//opts.print();

// defining timenumber
string timenumber = "";//currentDateTime();

// getting argv inputs
if (argc==2) timenumber = argv[1];
else if (argc % 2 && argc>1) {
	for (unsigned int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("loops")==0) opts.loops = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("printChoice")==0) opts.printChoice = argv[2*j+2];
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}
else if (argc != 1) {
	cerr << "must provide an even number of inputs in format '-name value':" << endl;
	for (int j=1; j<argc; j++) cerr << argv[j] << " ";
	cerr << endl;
	return 1;
}

/*----------------------------------------------------------------------------------------------------------------------------
	2. getting inputs, checks etc.
		- loading inputs
		- printing timenumber and parameters
		- clock
		- checks
----------------------------------------------------------------------------------------------------------------------------*/

// loading inputs
Parameters ps;
ps.load("inputsS");
if (ps.pot==3) {
	cerr << "no soliton for potential 3, retry with pot=1 or 2" << endl;
	return 1;
}
	
//printing timenumber and parameters
printf("%12s%12s\n","timenumber: ",timenumber.c_str());
printf("%14s%14s%14s%14s%14s%14s\n","N","L","dE","epsilon","minima[0]","minima[1]");
printf("%14i%14.4g%14.4g%14.4g%14.4g%14.4g\n",ps.N,ps.L,ps.dE,ps.epsilon,(ps.minima)[0],(ps.minima)[1]);
//ps.print();

//defining a time and starting the clock
clock_t time;
time = clock();

// loading closenesses
Closenesses closenesses;
closenesses.load("closenessesS");

// declaring Checks
Check checkMass("mass",closenesses.Action);
Check checkSoln("solution",closenesses.Soln);
Check checkProfile("phi input calculation",closenesses.Profile);

/*----------------------------------------------------------------------------------------------------------------------------
	3. assigning potential functions
		- assigning potential functions
		- assigning preliminary parameter structs
----------------------------------------------------------------------------------------------------------------------------*/

// assigning potential functions
if (ps.pot==1) {
	V((Potential<double>::PotentialType)&V1<double>,ps);
	dV((Potential<double>::PotentialType)&dV1<double>,ps);
	ddV((Potential<double>::PotentialType)&ddV1<double>,ps);
}
else if (ps.pot==2) {
	V((Potential<double>::PotentialType)&V2<double>,ps);
	dV((Potential<double>::PotentialType)&dV2<double>,ps);
	ddV((Potential<double>::PotentialType)&ddV2<double>,ps);
}
else {
	cerr << "pot option not available, pot = " << ps.pot << endl;
	return 1;
}

// assigning preliminary parameter structs
params_for_V paramsV0  = {ps.epsilon0, ps.A};
//params_for_V paramsV  = {ps.epsilon, ps.A};

/*----------------------------------------------------------------------------------------------------------------------------
	4. defining quantities
		- erg, linErg etc
		- p, minusDS, DDS
----------------------------------------------------------------------------------------------------------------------------*/

// Mass, posZero
double Mass;
double posZero;

//defining shooting test
double F, dF;
double aim = ps.minima[0];

//initializing phi (=p), dphi/dx (=dp)
vec p(ps.N+1);
p = Eigen::VectorXd::Zero(ps.N+1);
vec dp(ps.N+1);
dp = Eigen::VectorXd::Zero(ps.N+1);

/*----------------------------------------------------------------------------------------------------------------------------
	5. shooting method
		
----------------------------------------------------------------------------------------------------------------------------*/

// initial guess
double Y1 = -1.0e-2;
//cout << "initial y1: ";
//cin >> Y1;

// shooting loop
while (!checkSoln.good()) {
	runsCount++;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_yp_new (&sys,  gsl_odeiv2_step_rk8pd, 1.0e-6, 1.0e-6, 0.0);

	double y0[2] = { ps.minima[1], Y1};
	double y[2];
	memcpy(y, y0, sizeof(y0));
	double yMin = abs(y[0]-aim);
	p[0] = y[0];
	dp[0] = y[1];
	int status;
	unsigned int i, iMin = 0;
	
	for (i = 1; i < ps.N; i++) {
		status = gsl_odeiv2_driver_apply (d, &r, ri, y);
		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			printf ("i = %3i, r = %3g\n",i,r);
			break;
		}
		if ((y[0]-aim)<yMin && (y[0]-aim)>0.0) {
			yMin = y[0]-aim;
			iMin = i;
		}
		p[i] = y[0];
		dp[i] = y[1];
		//printf ("%.5e %.5e %.5e\n", r, y[0], y[1]);
		if ((y[0]-aim)<0.0) {
			iMin = i;
			if ((y[0]-aim)<-0.2) break;
		}
	}
	if (status != GSL_SUCCESS) break;
		
	F = y0Vec[iMin]-aim; //as final boundary condition is y(r1)=0.0;
	dF = y1Vec[iMin];
	printf("%16i%16i%16g%16g%16g%16.12g",runsCount,i,y[0],yMin,F-aim,Y0);
	if (abs(dF)>MIN_NUMBER) {
		Y0 += -F/dF;
		printf("%16.12g%16g\n",Y0,-F/dF);
	}
	gsl_odeiv2_driver_free (d);
	if (i==ps.N) F = y[0];
	checkSoln.add(absDiff(F,aim));
}

/*----------------------------------------------------------------------------------------------------------------------------
	6. calculating mass
		
----------------------------------------------------------------------------------------------------------------------------*/

//finding E
double E, Eerror;
E_params paramsE;
paramsE.Y_0 = Y0;
gsl_function E_integrand_gsl;
E_integrand_gsl.function = &E_integrand;
E_integrand_gsl.params = &paramsE;
gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
gsl_integration_qag(&E_integrand_gsl, -ps.L/2.0, ps.L/2.0, 1.0e-10, 1.0e-9, 1e4, 4, w, &E, &Eerror);
gsl_integration_workspace_free(w);
if (Eerror>1.0e-8) 	cout << "E error = " << Eerror << endl;
else 				cout << "E = " << E << endl << endl;


/*----------------------------------------------------------------------------------------------------------------------------
	7. printing results
		- stopping clock
		- printing results to terminal
		- printing results to file
		- printing (and plotting if vectors):
			- p
			- dp
		- return
	
----------------------------------------------------------------------------------------------------------------------------*/

//stopping clock
time = clock() - time;
double realtime = time/1000000.0;

//printing results to terminal
printf("\n");
printf("%8s%8s%8s%12s%14s%14s%14s\n","runs","time","N","L","dE","posZero","Mass");
printf("%8i%8.1g%8i%12g%14.4g%14.4g%14.4g\n",runs_count,realtime,ps.N,ps.L,ps.dE,posZero,Mass);
printf("\n");
printf("%60s\n","----------------------------------------------------------------------------------------------------");

//printing results to file
FILE * staticfile;
staticfile = fopen("data/staticShooting.dat","a");
fprintf(staticfile,"%16s%8i%12g%12g%14.4g%14.4g%14.4g\n",timenumber.c_str()\
			,ps.N,ps.L,ps.dE,posZero,Mass,checkSoln.back());
fclose(staticfile);

bool printEverything = false;

string prefix = "data/"+timenumber;
string suffix = ".dat";
so_simple.printMessage = true;

//printing output phi on Euclidean time part
Filename pFile = (string)(prefix+"staticShootingp"+suffix);
save(pFile,so_simple,p);
Filename plotFile = pFile;
plotFile.Directory = "pics";
plotFile.Suffix = ".png";
po_simple.output = plotFile;
plot(pFile,po_simple);

if (printEverything) {
	//printing output dp
	pFile.ID = "staticShootingdp";
	save(pFile,so_simple,dp);
}

if (!checkDelta.good()) {
	return 1;
}

return 0;
}

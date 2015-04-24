/*----------------------------------------------------------------------------------------------------------------------------
	static
		program to solve boundary value problem for static soliton
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
#include "main.h"

//#define NDEBUG //NDEBUG is to remove error and bounds checking on vectors in SparseLU, for speed - only include once everything works

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - loading options, closenesses, argv inputs
		2 - getting inputs, checks etc.
		3 - assigning potential fucntions
		4 - defining quantites
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/


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
ps.print();

//defining a time and starting the clock
clock_t time;
time = clock();

// loading closenesses
Closenesses closenesses;
closenesses.load("closenessesS");

// declaring Checks
Check checkMass("mass",closenesses.Action);
Check checkSoln("solution",closenesses.Soln);
Check checkSolnMax("solution max",closenesses.SolnMax);
Check checkDelta("delta",closenesses.Delta);
Check checkInv("matrix inversion",closenesses.Inv*ps.N*ps.NT);
Check checkProfile("phi input calculation",closenesses.Profile);

// do trivial or redundant checks?
bool trivialChecks = false;

/*----------------------------------------------------------------------------------------------------------------------------
	3. assigning potential functions
		- assigning potential functions
		- assigning preliminary parameter structs
----------------------------------------------------------------------------------------------------------------------------*/

// assigning potential functions
Potential<double> V, dV, ddV;
if (ps.pot==1) {
	V((Potential<double>::PotentialType)&V1<comp>,ps);
	dV((Potential<double>::PotentialType)&dV1<comp>,ps);
	ddV((Potential<double>::PotentialType)&ddV1<comp>,ps);
}
else if (ps.pot==2) {
	V((Potential<double>::PotentialType)&V2<comp>,ps);
	dV((Potential<double>::PotentialType)&dV2<comp>,ps);
	ddV((Potential<double>::PotentialType)&ddV2<comp>,ps);
}
else {
	cerr << "pot option not available, pot = " << ps.pot << endl;
	return 1;
}

// assigning preliminary parameter structs
params_for_V paramsV  = {ps.epsilon, ps.A}, paramsV0  = {ps.epsilon0, ps.A};

/*----------------------------------------------------------------------------------------------------------------------------
	4. defining quantities
		- erg, linErg etc
		- p, minusDS, DDS
----------------------------------------------------------------------------------------------------------------------------*/

// Mass
double Mass;	

//defining some quantities used to stop the Newton-Raphson loop when action stops varying
uint runs_count = 0;
uint min_runs = 3;

//initializing phi (=p), DDS and minusDS
vec p(ps.N+1);
p = Eigen::VectorXd::Zero(ps.N+1);
spMat DDS(ps.N+1,ps.N+1);
vec minusDS(ps.N+1);

/*----------------------------------------------------------------------------------------------------------------------------
	5. calculating input phi
		- finding phi profile between minima
		- assigning phi
		- fixing boundary conditions
		- printing input phi
		- Cp
----------------------------------------------------------------------------------------------------------------------------*/
	

//finding phi profile between minima
uint profileSize = ps.N; //more than the minimum
vector<double> phiProfile(profileSize);
vector<double> rhoProfile(profileSize);
double alphaL = -opts.alpha, alphaR = opts.alpha;
if (ps.R<opts.alpha) {
	cerr << "R is too small. Not possible to give thinwall input. It should be >> " << opts.alpha;
	return 1;
}
if (ps.pot==2) {
	double phiL = ps.minima0[1]-1.0e-2;
	double phiR = ps.minima0[0]+1.0e-2;
	for (uint j=0;j<profileSize;j++) {
		phiProfile[j] = phiL + (phiR-phiL)*j/(profileSize-1.0);
	}

	double profileError;
	gsl_function rho_integrand;
	rho_integrand.function = &rhoIntegrand;
	rho_integrand.params = &paramsV0;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
	w = gsl_integration_workspace_alloc(1e4);
	for (uint j=0;j<profileSize;j++) {
		gsl_integration_qags(&rho_integrand, phiProfile[j], 0, 1.0e-16, 1.0e-6, 1e4, w,\
												 &(rhoProfile[j]), &profileError);
		checkProfile.add(profileError);
		checkProfile.checkMessage();
	}
	gsl_integration_workspace_free(w);
	alphaL = rhoProfile[0];
	alphaR = rhoProfile.back();
}
for (uint j=0; j<ps.N; j++) {
	double x = -ps.L/2.0+ps.a*(double)j;
	if (x>alphaR) {
		p(j) = ps.minima[0];
	}
	else if (x<alphaL) {
		p(j) = ps.minima[1];
	}
	else if (ps.pot==2) {
		vector<double> rhoPos (profileSize,x);
		for (uint k=0; k<profileSize; k++) {
			rhoPos[k] -= rhoProfile[k];
		}
		uint minLoc = smallestLoc(rhoPos);
        p(j) = phiProfile[minLoc];
	}
	else {
		p(2*j) = (ps.minima[1]+ps.minima[0])/2.0\
					+ (ps.minima[0]-ps.minima[1])*tanh(x/2.0)/2.0;
	}
}
p(ps.N) = 0.5; // lagrange multiplier for zero mode

// printing input phi
SaveOptions so_simple;
so_simple.paramsIn = ps; so_simple.paramsOut = ps;
so_simple.vectorType = SaveOptions::simple;
so_simple.extras = SaveOptions::none;
so_simple.printMessage = true;
Filename earlyFile = (string)("data/"+timenumber+"staticE.dat");
save(earlyFile,so_simple,p);

// plotting input phi
PlotOptions po_simple;
po_simple.column = 1;
po_simple.style = "linespoints";
po_simple.printMessage = true;
Filename earlyPlotFile = (string)("data/"+timenumber+"staticE.png");
po_simple.output = earlyPlotFile;
plot(earlyFile,po_simple);

return 0;
}

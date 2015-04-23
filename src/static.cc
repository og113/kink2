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
string timenumber = currentDateTime();

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
ps.load("inputsP");
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
Check checkAction("mass",closenesses.Action);
Check checkSoln("solution",closenesses.Soln);
Check checkSolnMax("solution max",closenesses.SolnMax);
Check checkDelta("delta",closenesses.Delta);
Check checkInv("matrix inversion",closenesses.Inv*ps.N*ps.NT);

// do trivial or redundant checks?
bool trivialChecks = false;

/*----------------------------------------------------------------------------------------------------------------------------
	3. assigning potential functions
		- assigning potential functions
		- assigning preliminary parameter structs
----------------------------------------------------------------------------------------------------------------------------*/

// assigning potential functions
Potential V, dV, ddV;
if (ps.pot==1) {
	V((PotentialType)&V1<comp>,ps);
	dV((PotentialType)&dV1<comp>,ps);
	ddV((PotentialType)&ddV1<comp>,ps);
}
else if (ps.pot==2) {
	V((PotentialType)&V2<comp>,ps);
	dV((PotentialType)&dV2<comp>,ps);
	ddV((PotentialType)&ddV2<comp>,ps);
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

// erg, linErg etc
comp action = ps.action0;
double W;
double E;
double ergZero = 0.0;
cVec erg(ps.NT);	

//defining some quantities used to stop the Newton-Raphson loop when action stops varying
comp action_last = action;
uint runs_count = 0;
uint min_runs = 3;

//initializing phi (=p), DDS and minusDS
vec p(2*ps.N*ps.Nb+1);
p = Eigen::VectorXd::Zero(2*ps.N*ps.Nb+1);
spMat DDS(2*ps.N*ps.Nb+1,2*ps.N*ps.Nb+1);
vec minusDS(2*ps.N*ps.Nb+1);



return 0;
}

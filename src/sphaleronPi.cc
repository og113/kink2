/* ---------------------------------------------------------------------------------------------
	sphaleron_pi
		provides initial guess for pi consisting of
		sphaleron plus small amount of negative mode, constant in time
---------------------------------------------------------------------------------------------*/
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
#include <ctype.h>
#include <cstring> //for memcpy
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "simple.h"
#include "folder.h"
#include "parameters.h"
#include "print.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - defining main parameters
		2 - getting argv inputs
		3 - loading sphaleron and negEig
		4 - normalising
		5 - constructing intial guess for phi
		6 - printing result
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char ** argv) {
/* ---------------------------------------------------------------------------------------------
	1. defining main parameters
---------------------------------------------------------------------------------------------*/
unsigned int 	N = 1e3, Nt = 1e3;
double 			r0 = 1.0e-16, r1 = 5.0, t0 = 0.0, t1 = 0.80;
double			dr, dt;
double			amp = 0.5;
Parameters		ps_run;

/* ---------------------------------------------------------------------------------------------
	2. getting argv inputs
---------------------------------------------------------------------------------------------*/
if (argc>2 && argc%2) {
	for (int j=0; j<(int)(argc/2); j++) {
		string temp1 = argv[2*j+1];
		string temp2 = argv[2*j+2];
		if (temp1[0]=='-') temp1 = temp1.substr(1);
		if (temp1.compare("amp")==0) amp = stringToNumber<double>(temp2);
		else if (temp1.compare("t1")==0) t1 = stringToNumber<double>(temp2);
		else if (temp1.compare("r1")==0) r1 = stringToNumber<double>(temp2);
		else if (temp1.compare("r0")==0) r0 = stringToNumber<double>(temp2);
		else if (temp1.compare("N")==0) N = stringToNumber<unsigned int>(temp2);
		else if (temp1.compare("Nt")==0) Nt = stringToNumber<unsigned int>(temp2);
		else {
			cerr << "input " << temp1 << " not understood" << endl;
			return 1;
		}
	}
}
else if (argc==2) amp = atof(argv[1]);
else if (argc>1) cerr << "inputs not understood" << endl;

dr = (r1-r0)/(double)(N-1.0), dt = (t1-t0)/(double)(Nt-1.0);

printf("%8s%8s%8s%8s%8s%8s%8s\n","N","Nt","r1","r0","t1","t0","amp");
printf("%8i%8i%8g%8g%8g%8g%8g\n",N,Nt,r1,r0,t1,t0,amp);
printf("\n");

ps_run.pot = 3;
ps_run.N = N;
ps_run.Nb = Nt;
ps_run.r0 = r0;
ps_run.L = r1;
ps_run.Tb = t1;
ps_run.a = dr;
ps_run.b = dt;

/* ---------------------------------------------------------------------------------------------
	3. loading sphaleron and negEig
---------------------------------------------------------------------------------------------*/
vec sphaleron(N), negEig(N);

Filename sphaleronFile = (string)("data/stable/sphaleron_L_"+numberToString<double>(r1)+".dat");
Filename negEigFile = (string)("data/stable/eigVec_pot_3_L_"+numberToString<double>(r1)+".dat");
SaveOptions so_simple;
so_simple.paramsOut = ps_run;
so_simple.vectorType = SaveOptions::simple;
so_simple.extras = SaveOptions::none;
so_simple.printMessage = false;
load(sphaleronFile,so_simple,sphaleron);
load(negEigFile,so_simple,negEig);

/* ---------------------------------------------------------------------------------------------
	4. normalising
		- normalising sign
		- normalising magnitude
---------------------------------------------------------------------------------------------*/
if (negEig[0]<0) negEig *= -1.0;
double normSphaleron, normNegEig;
normSphaleron = sphaleron.norm();
normNegEig = negEig.norm();
negEig *= normSphaleron/normNegEig;

/* ---------------------------------------------------------------------------------------------
	5. constructing intial guess for phi
		- defining phi
		- assiging values
		
	n.b. time direction is reversed
---------------------------------------------------------------------------------------------*/

vec phi(N*Nt);

for (unsigned int k=0; k<Nt; k++) {
	for (unsigned int j=0; j<N; j++) {
		unsigned int 	m = (Nt-1-k) + j*Nt;
		double 			r = r0 + dr*j, t = t0 + dt*(Nt-1.0-k);
		phi[m] = r*(sphaleron[j] + amp*cos(3.91*t)*negEig[j]);
	}
}
	
/* ---------------------------------------------------------------------------------------------
	6. printing result
---------------------------------------------------------------------------------------------*/
Filename datOut = (string)("data/sphaleronPi_L_"+numberToString<double>(r1)+"_Tb_"+numberToString<double>(t1)+".dat");
Filename pngOut = datOut;
pngOut.Suffix = ".png";
pngOut.Directory = "pics";

Parameters ps_print;
ps_print = ps_run;
ps_print.N = 300; // must have output N=Nb as load into pi requires this
ps_print.Nb = 300;
ps_print.a = (r1-r0)/(ps_print.N-1.0);
ps_print.b = (t1-t0)/(ps_print.Nb-1.0);

SaveOptions so_p;
so_p.printMessage = true;
so_p.vectorType = SaveOptions::realB;
so_p.extras = SaveOptions::coords;
so_p.paramsIn = ps_run;
so_p.paramsOut = ps_print;

PlotOptions po_p;
po_p.printMessage = true;
po_p.gp = "gp/repi.gp";
po_p.style = "points";
po_p.output = pngOut;

save(datOut,so_p,phi);
plot(datOut,po_p);

return 0;
}

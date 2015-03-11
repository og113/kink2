	/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for functions to save, load and plot
 -------------------------------------------------------------------------------------------------------------------------*/
 
 #include <fstream>
 #include <iostream>
 #include <sstream>
 #include <iomanip>
 #include <cmath>
 #include <complex>
#include <Eigen/Dense>
#include "simple.h" // for countLines
#include "parameters.h"
#include "lattice.h"
#include "print.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. saveOptions
	2. save
	3. load
	4. plot
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	2. save
		- vec
			- static simpleSaveVec
			- static saveVecB
			- static saveVec
		- cVec
			- static simpleSaveVec
			- static saveVecB
			- static saveVec
		- mat
		- cMat
		- spMat
-------------------------------------------------------------------------------------------------------------------------*/

// save vec - simpleSaveVec
static void simpleSaveVector(const string& f, const saveOptions& opts, const vec& v) {
	fstream F;
	F.open((f).c_str(), ios::out);
	F.precision(16);
	F << left;
	uint length = vecToPrint.size();
	for (uint j=0; j<length; j++) {
		if (opts.extras==loc) 	F << setw(25) << j;
								F << setw(25) << v(j) << endl;
	}
	F.close();
}

// save vec - saveVecB
static void saveVecB (const string& f, const saveOptions& opts, const vec& v) {
	Parameters pi = opts.ParamsIn, po = opts.ParamsOut;
	vec vo;
	if ((pi.N!=po.N && po.N!=0) || (pi.Nb!=po.Nb && po.Nb!=0)) {
		vo = interpolate(v,pi,po);
	}
	else {
		vo = v;
	}
	fstream F;
	F.open((f).c_str(), ios::out);
	int x0 = intCoord(0,1,po.Nb);
	F.precision(16);
	for (lint j=0; j<po.Nb; j++) {
		uint x = intCoord(j,1,po.Nb);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		F << left;
		F << setw(24) << real(coordB(j,0,po)) << setw(25) << imag(coordB(j,0,po));
		F << setw(25) << real(coordB(j,1,po));
		if (vecToPrint.size()>po.N*po.Nb) {
			F << setw(25) << vo(2*j) << setw(25) << vo(2*j+1)  << endl;
		}
		else {
			F << setw(25) << vo(j) << endl;
		}
	}
	if (vo.size()>2*po.N*po.Nb) {
		F << endl;
		for (unsigned int k=0; k<(vo.size()-2*po.N*po.Nb);k++) {
			F << setw(25) << vo(2*po.N*po.Nb+k) << endl;
		}
	}
	F.close();
}

// save vec
void save(const string&, const saveOptions&, const vec&);

// save cVec - simpleSaveVec
void simpleSaveVector(const string& f, const saveOptions& opts, const cVec& v) {
	fstream F;
	F.open((f).c_str(), ios::out);
	F.precision(16);
	F << left;
	unsigned int length = vecToPrint.size();
	for (unsigned int j=0; j<length; j++) {
		if (opts.extras==loc) 	F << setw(25) << j;
								F << setw(25) << real(v(j)) << setw(25) << imag(v(j)) << endl;
	}
	F.close();
}

// save cVec
void save(const string&, const saveOptions&, const cVec&);

// save mat
void save(const string&, const saveOptions&, const mat&);

// save cMat
void save(const string&, const saveOptions&, const cMat&);

// save spMat
void save(const string&, const saveOptions&, const spMat&);

/*-------------------------------------------------------------------------------------------------------------------------
	3. load
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	4. plot
-------------------------------------------------------------------------------------------------------------------------*/

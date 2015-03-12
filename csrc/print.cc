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
#include <Eigen/Sparse>
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
			- static savecVecSimple
			- static savecVecB
			- static savecVec
		- mat
		- cMat
		- spMat
-------------------------------------------------------------------------------------------------------------------------*/

// save vec - simpleSaveVec
static void saveVecSimple(const string& f, const saveOptions& opts, const vec& v) {
	fstream F;
	F.open(f.c_str(), ios::out);
	F.precision(16);
	F << left;
	uint length = vecToPrint.size();
	if (opts.vectorType==complex && !length%2) length = (uint)(length/2);
	for (uint j=0; j<length; j++) {
		if (opts.extras==loc) 			F << setw(25) << j;
		if (opts.vectorType==complex)	F << setw(25) << v(2*j) << setw(25) << v(2*j+1);
		else							F << setw(25) << v(j);
										F << endl;
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
	F.open(f.c_str(), ios::out);
	int x0 = intCoord(0,1,po.Nb);
	F.precision(16);
	F << left;
	for (lint j=0; j<po.N*po.Nb; j++) {
		uint x = intCoord(j,1,po.Nb);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		switch(objs.extras) {
			case none:	break;
			case loc: 	F << setw(25) << j;
						break;
			case coord:	F << setw(25) << real(coordB(j,0,po)) << setw(25) << imag(coordB(j,0,po));
						F << setw(25) << real(coordB(j,1,po)); 
						break;
			default:	cerr << "save error: print extras option(" << objs.extras << ") not possible" << endl;
						break;
		}
		switch(objs.vectorType) {
			case realB:		F << setw(25) << vo(j) << endl;
							break;
			case complexB:	F << setw(25) << vo(2*j) << setw(25) << vo(2*j+1)  << endl;
							break;
			default:		cerr << "save error: print vectorType option(" << objs.vectorType << ") not possible" << endl;
							break;
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

// save vec - saveVec
static void saveVec(const string& f, const saveOptions&, const vec& v) {
	Parameters pi = opts.ParamsIn, po = opts.ParamsOut;
	vec vo;
	if ((pi.N!=po.N && po.N!=0) || (pi.NT!=po.NT && po.NT!=0)) {
		vo = interpolate(v,pi,po);
	}
	else {
		vo = v;
	}
	fstream F;
	F.open(f.c_str(), ios::out);
	unsigned int x0 = intCoord(0,1,po.NT);
	F.precision(16);
	F << left;
	for (unsigned long int j=0; j<po.N*po.NT; j++) {
		unsigned int x = intCoord(j,1,po.NT);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		switch(objs.extras) {
			case none:	break;
			case loc: 	F << setw(25) << j;
						break;
			case coord:	F << setw(25) << real(coord(j,0,po)) << setw(25) << imag(coord(j,0,po));
						F << setw(25) << real(coord(j,1,po)); 
						break;
			default:	cerr << "save error: print extras option(" << objs.extras << ") not possible" << endl;
						break;
		}
		switch(objs.vectorType) {
			case realB:		F << setw(25) << vo(j) << endl;
							break;
			case complexB:	F << setw(25) << vo(2*j) << setw(25) << vo(2*j+1)  << endl;
							break;
			default:		cerr << "save error: print vectorType option(" << objs.vectorType << ") not possible" << endl;
							break;
		}
	}
	if (vecToPrint.size()>2*N*NT) {
		F << endl;
		for (unsigned int j=0; j<(vecToPrint.size()-2*N*NT);j++) {
			F << setw(25) << vecToPrint(2*N*NT+j) << endl;
		}
	}
	F.close();
}

// save vec
void save(const string& f, const saveOptions& opts, const vec& v) {
	switch(vectorType) {
		case simple:	saveVecSimple(f, opts, v);
						break;
		case real:		saveVec(f,opts,v);
						break;
		case complex:	saveVec(f,opts,v);
						break;
		case realB:		saveVecB(f,opts,v);
						break;
		case complexB:	saveVecB(f,opts,v);
						break;
		default:		cerr << "save error: print vectorType option(" << objs.vectorType << ") not possible" << endl;
						break;
	}	
}

// save cVec - simpleSaveVec
static void savecVecSimple(const string& f, const saveOptions& opts, const cVec& v) {
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

// save cVec - saveVecB
static void saveVecB (const string& f, const saveOptions& opts, const cVec& v) {
	Parameters pi = opts.ParamsIn, po = opts.ParamsOut;
	cVec vo;
	if ((pi.N!=po.N && po.N!=0) || (pi.Nb!=po.Nb && po.Nb!=0)) {
		vo = interpolate(v,pi,po);
	}
	else {
		vo = v;
	}
	fstream F;
	F.open(f.c_str(), ios::out);
	int x0 = intCoord(0,1,po.Nb);
	F.precision(16);
	F << left;
	for (lint j=0; j<po.N*po.Nb; j++) {
		uint x = intCoord(j,1,po.Nb);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		switch(objs.extras) {
			case none:	break;
			case loc: 	F << setw(25) << j;
						break;
			case coord:	F << setw(25) << real(coordB(j,0,po)) << setw(25) << imag(coordB(j,0,po));
						F << setw(25) << real(coordB(j,1,po)); 
						break;
			default:	cerr << "save error: print extras option(" << objs.extras << ") not possible" << endl;
						break;
		}
		F << setw(25) << real(vo(j)) << setw(25) << imag(vo(j))  << endl;
	}
	if (vo.size()>po.N*po.Nb) {
		F << endl;
		for (unsigned int k=0; k<(vo.size()-po.N*po.Nb);k++) {
			F << setw(25) << vo(po.N*po.Nb+k) << endl;
		}
	}
	F.close();
}

// save cVec - saveVec
static void saveVec(const string& f, const saveOptions&, const cVec& v) {
	Parameters pi = opts.ParamsIn, po = opts.ParamsOut;
	cVec vo;
	if ((pi.N!=po.N && po.N!=0) || (pi.NT!=po.NT && po.NT!=0)) {
		vo = interpolate(v,pi,po);
	}
	else {
		vo = v;
	}
	fstream F;
	F.open(f.c_str(), ios::out);
	unsigned int x0 = intCoord(0,1,po.NT);
	F.precision(16);
	F << left;
	for (unsigned long int j=0; j<po.N*po.NT; j++) {
		unsigned int x = intCoord(j,1,po.NT);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		switch(objs.extras) {
			case none:	break;
			case loc: 	F << setw(25) << j;
						break;
			case coord:	F << setw(25) << real(coord(j,0,po)) << setw(25) << imag(coord(j,0,po));
						F << setw(25) << real(coord(j,1,po)); 
						break;
			default:	cerr << "save error: print extras option(" << objs.extras << ") not possible" << endl;
						break;
		}
		F << setw(25) << real(vo(j)) << setw(25) << imag(vo(j))  << endl;
	}
	if (vecToPrint.size()>N*NT) {
		F << endl;
		for (unsigned int j=0; j<(vecToPrint.size()-N*NT);j++) {
			F << setw(25) << vecToPrint(N*NT+j) << endl;
		}
	}
	F.close();
}

// save cVec
void save(const string&, const saveOptions&, const cVec&) {
	switch(vectorType) {
		case simple:	saveVecSimple(f, opts, v);
						break;
		case real:		saveVec(f,opts,v);
						break;
		case complex:	saveVec(f,opts,v);
						break;
		case realB:		saveVecB(f,opts,v);
						break;
		case complexB:	saveVecB(f,opts,v);
						break;
		default:		cerr << "save error: print vectorType option(" << objs.vectorType << ") not possible" << endl;
						break;
	}	
}

// save mat
void save(const string& f, const saveOptions& opts, const mat& m) {
	fstream F;
	F.open((f).c_str(), ios::out);
	F << left;
	F.precision(16);
	for (uint j=0; j<m.rows(); j++)
		{
		for (uint k=0; k<m.cols(); k++)
			{
			F << setw(25) << m(j,k) << endl;
			}
		}
	F.close();
}

// save cMat
void save(const string&, const saveOptions&, const cMat&) {
	fstream F;
	F.open((f).c_str(), ios::out);
	F << left;
	F.precision(16);
	for (uint j=0; j<m.rows(); j++)
		{
		for (uint k=0; k<m.cols(); k++)
			{
			F << setw(25) << real(m(j,k)) << setw(25) << imag(m(j,k)) << endl;
			}
		}
	F.close();
}

// save spMat
void save(const string& f, const saveOptions& opts, const spMat& m) {
	fstream F;
	F.open((f).c_str(), ios::out);
	F << left;
	F.precision(16);
	for (int l=0; l<m.outerSize(); ++l)
		{
		for (Eigen::SparseMatrix<double>::InnerIterator it(m,l); it; ++it)
			{
			F << setw(25) << it.row()+1 << setw(25) << it.col()+1 << setw(25) << it.value() << endl;
			}
		}
	F.close();
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. load
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	4. plot
-------------------------------------------------------------------------------------------------------------------------*/

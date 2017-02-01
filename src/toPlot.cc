/*----------------------------------------------------------------------------------------------------------------------------
	toPlot
		program to take in binary files and spit out ascii files ready to plot
----------------------------------------------------------------------------------------------------------------------------*/
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <string>
#include "simple.h"
#include "error.h"
#include "folder.h"
#include "parameters.h"
#include "print3.h"
#include "lattice.h"
#include "nr.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - defining main quantities
		2 - getting argv inputs
		3 - loading input
		4 - saving output
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/* -------------------------------------------------------------------------------------------------------------------------
	1. defining main quantities
-------------------------------------------------------------------------------------------------------------------------*/
Filename fi, fo;
Parameters p;

/* -------------------------------------------------------------------------------------------------------------------------
	2. getting argv inputs
-------------------------------------------------------------------------------------------------------------------------*/

// getting argv inputs
if (argc==2)
	fi = (string)argv[1];
else if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("fi")==0 || id.compare("fileIn")==0) 						fi = (string)argv[2*j+2];
		else if (id.compare("fo")==0 || id.compare("fileOut")==0) 					fo = (string)argv[2*j+2];
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

if (fi.empty()) {
	cerr << "Filenames empty:" << fi << ", " << fo << endl;
	return 1;
}

if (fo.empty()) {
	fo = fi;
	fo.Directory = "temp";
	fo.Suffix = ".dat";
}

cout << "input:  " << fi << endl;
cout << "output: " << fo << endl;

/* -------------------------------------------------------------------------------------------------------------------------
	3. loading input
-------------------------------------------------------------------------------------------------------------------------*/

// getting parameters from filename
p = filenameToParameters(fi);

// loading vector
vec v;
loadVectorBinary(fi,v);

/* -------------------------------------------------------------------------------------------------------------------------
	4. saving output
-------------------------------------------------------------------------------------------------------------------------*/

ofstream os;
os.open(((string)fo).c_str());
os << left << setprecision(16) << endl;

if (p.Na==0 && p.Nb==0 && p.Nc==0 && p.N>0 && p.LoR>0) {
	if (v.size()>=p.N) {
		// printing spatial vector
		double r;
		for (uint j=0; j<p.N; j++) {
			r = p.r0+j*p.a;
			os << setw(25) << r << setw(25) << v[j] << endl;
		}
	}
	else {
		for (uint j=0; j<v.size(); j++) {
			os << setw(25) << v[j] << endl;
		}
	}
}
else if (p.Na==0 && p.Nc==0 && p.Nb>0 && p.N>0&& p.LoR>0 && p.Tb>0) {
	if (v.size()>=p.N*p.Nb) {
		// printing euclidean time vector
		uint x;
		uint x0 = intCoord(0,1,p.Nb);
		for (uint j=0; j<p.N*p.Nb; j++) {
			x = intCoord(j,1,p.Nb);
			if (x!=x0) { //this is put in for gnuplot
				os << endl;
				x0 = x;
			}
			os << setw(25) << real(coordB(j,0,p));
			os << setw(25) << imag(coordB(j,0,p));
			os << setw(25) << real(coordB(j,1,p));
			os << setw(25) << v[2*j] << setw(25) << v[2*j+1]  << endl;
		}
	}
	else {
		for (uint j=0; j<v.size(); j++) {
			os << setw(25) << v[j] << endl;
		}
	}
}
else if (p.Na>0 && p.Nb>0 && p.Nc>0 && p.N>0 && p.LoR>0 && p.Tb>0) {
	if (v.size()>=p.N*p.NT) {
		// printing vector over full contour
		uint x;
		uint x0 = intCoord(0,1,p.NT);
		for (uint j=0; j<p.N*p.NT; j++) {
			x = intCoord(j,1,p.NT);
			if (x!=x0) { //this is put in for gnuplot
				os << endl;
				x0 = x;
			}
			os << setw(25) << real(coord(j,0,p));
			os << setw(25) << imag(coord(j,0,p));
			os << setw(25) << real(coord(j,1,p));
			os << setw(25) << v[2*j] << setw(25) << v[2*j+1]  << endl;
		}
	}
	else if (v.size()>=2*p.NT) {
		// energy plot
		for (uint j=0; j<p.NT; j++) {
			os << setw(25) << v[2*j];
			os << setw(25) << v[2*j+1] << endl;
		}
	}
	else {
		for (uint j=0; j<v.size(); j++) {
			os << setw(25) << v[j] << endl;
		}
	}
}
else {
	cerr << "toPlot error: not able to pull parameters from filename" << endl;
	cerr << fi << endl;
	cerr << p << endl;
}

os.close();

return 0;
}

/*-------------------------------------------------------------------------------------------------------------------------
 	program to analyse data
 -------------------------------------------------------------------------------------------------------------------------*/

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <utility> // for pair
#include <cstdlib> // for system
#include <algorithm> // for sort
#include <iterator>
#include <Eigen/Dense>
#include "simple.h"
#include "folder.h"
#include "print.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. defining key parameters
	2. argv inputs
	3. getting data from file
	4. forming derivatives
	5. printing results
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. defining key parameters
-------------------------------------------------------------------------------------------------------------------------*/

Filename f = (string)("results/main_pot_3.dat");
mat data;
uint colTb = 5, colTheta = 7, colE = 8, colN = 9, colS = 10, colW = 11;
uint colX, colY, colRhs, colRhs2 = 0, colConst;
string deriv = "dSdT|theta";

/*-------------------------------------------------------------------------------------------------------------------------
	2. getting argv inputs
-------------------------------------------------------------------------------------------------------------------------*/

// getting argv inputs
if (argc == 2) f = (string)(argv[1]);
else if (argc % 2 && argc>1) {
	for (unsigned int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("f")==0 || id.compare("file")==0) f = (string)(argv[2*j+2]);
		else if (id.compare("d")==0 || id.compare("deriv")==0) deriv = argv[2*j+2];
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}
else if (argc!=1) {
	cerr << "inputs not understood:" << endl;
	for (int j=1; j<argc; j++) cerr << argv[j] << " ";
	cerr << endl;
	return 1;
}

if (deriv.compare("dWdE|N")==0) {
	colX = colW;
	colY = colE;
	colRhs = colTb;
	colConst = colN;
}
else if (deriv.compare("dWdN|E")==0) {
	colX = colW;
	colY = colN;
	colRhs = coltheta;
	colConst = colE;
}
else if (deriv.compare("dNdE|W")==0) {
	colX = colN;
	colY = colE;
	colRhs = colTb;
	colRhs2 = colTheta;
	colConst = colW;
}
else if (deriv.compare("dSdT|theta")==0) {
	colX = colS;
	colY = colTb;
	colRhs = colE;
	colConst = colTheta;
}
else if (deriv.compare("dSdtheta|T")==0) {
	colX = colS;
	colY = coltheta;
	colRhs = colN;
	colConst = colTb;
}
else {
	cerr << "derivative option: " << deriv << " not understood" << endl;
	return 1;
}


/*-------------------------------------------------------------------------------------------------------------------------
	3. getting data from file
-------------------------------------------------------------------------------------------------------------------------*/

SaveOptions so;
so_in.extras = SaveOptions::none;
so_in.printMessage = true;

load(f,so,data);

/*-------------------------------------------------------------------------------------------------------------------------
	4. forming derivatives
		- defining data by column
		- possible derivatives:
			- dWdE|N = -T
			- dWdN|E = -theta
			- dN/dE|W = -T/theta
			- dS/dT|theta = E
			- dS/dtheta|T = N
	
	n.b. T = 2Tb
-------------------------------------------------------------------------------------------------------------------------*/

vec dXdY_Z(data.rows()-1), rhs(data.rows()-1), diff(data.rows()-1);

for (uint j=0; j<(data.rows()-1); j++) {
	dXdY_Z(j) = (data(j+1,colX)-data(j,colX))/(data(j+1,colY)-data(j,colY));
	rhs(j) = (data(j+1,colRhs)+data(j,colRhs))/2.0;
	if (colRhs2!=0)
		rhs(j) /= (data(j+1,colRhs2)+data(j,colRhs2))/2.0;
}

if (deriv.compare("dWdE|N")==0 || deriv.compare("dWdN|E")==0 || deriv.compare("dNdE|W")==0)
	rhs *= -1.0;
	

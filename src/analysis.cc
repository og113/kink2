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

int main(int argc, char** argv) {

/*-------------------------------------------------------------------------------------------------------------------------
	1. defining key parameters

n.b. the first column is column 0
-------------------------------------------------------------------------------------------------------------------------*/

Filename f = (string)("results/mainReduced_pot_3.dat");
mat data;
uint colTb = 4, colTheta = 6, colE = 7, colN = 8, colS = 9, colW = 10;
uint colX, colY, colRhs, colRhs2 = 0, colConst;
string deriv = "dSdT-theta";
double closenessConst;

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
		else if (id.compare("c")==0 || id.compare("closeness")==0) closenessConst = stringToNumber<double>(argv[2*j+2]);
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

if (deriv.compare("dWdE-N")==0) {
	colX = colW;
	colY = colE;
	colRhs = colTb;
	colConst = colN;
	closenessConst = 5.0e-2;
}
else if (deriv.compare("dWdN-E")==0) {
	colX = colW;
	colY = colN;
	colRhs = colTheta;
	colConst = colE;
	closenessConst = 5.0e-2;
}
else if (deriv.compare("dNdE-W")==0) {
	colX = colN;
	colY = colE;
	colRhs = colTb;
	colRhs2 = colTheta;
	colConst = colW;
	closenessConst = 5.0e-2;
}
else if (deriv.compare("dSdT-theta")==0) {
	colX = colS;
	colY = colTb;
	colRhs = colE;
	colConst = colTheta;
	closenessConst = 1.0e-5;
}
else if (deriv.compare("dSdtheta-T")==0) {
	colX = colS;
	colY = colTheta;
	colRhs = colN;
	colConst = colTb;
	closenessConst = 1.0e-5;
}
else {
	cerr << "derivative option: " << deriv << " not understood" << endl;
	return 1;
}


/*-------------------------------------------------------------------------------------------------------------------------
	3. getting data from file
-------------------------------------------------------------------------------------------------------------------------*/

SaveOptions so;
so.extras = SaveOptions::none;
so.printMessage = true;

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

vec Y(data.rows()-1), dXdY_Z(data.rows()-1), rhs(data.rows()-1), diff(data.rows()-1), constError(data.rows()-1);

for (uint j=0; j<(data.rows()-1); j++) {
	Y(j) = (data(j+1,colY)+data(j,colY))/2.0;
	dXdY_Z(j) = (data(j+1,colX)-data(j,colX))/(data(j+1,colY)-data(j,colY));
	rhs(j) = (data(j+1,colRhs)+data(j,colRhs))/2.0;
	if (colRhs2!=0)
		rhs(j) /= (data(j+1,colRhs2)+data(j,colRhs2))/2.0;
	if (deriv.compare("dWdE-N")==0 || deriv.compare("dWdN-E")==0 || deriv.compare("dNdE-W")==0)
		rhs(j) *= -1.0;
	if (colRhs==colTb)
		rhs(j) *= 2.0;
	if (colRhs2==colTb)
		rhs(j) /= 2.0;
	diff(j) = absDiff(rhs(j),dXdY_Z(j));
	constError(j) = absDiff(data(j+1,colConst),data(j,colConst));
	if (abs(constError(j))>closenessConst) {
		if (j==0) {
			cerr << "required parameter not held constant" << endl;
			return 1;
		}
		Y.conservativeResize(j);
		dXdY_Z.conservativeResize(j);
		rhs.conservativeResize(j);
		diff.conservativeResize(j);
		constError.conservativeResize(j);
		break;
	}
}

	
/*-------------------------------------------------------------------------------------------------------------------------
	5. printing results
		- saving vectors to file
		- plotting dXdY_Z vs rhs
		- printing maxCoeffs
-------------------------------------------------------------------------------------------------------------------------*/

// saving vectors to file
Filename fo = f;
fo.ID = f.ID + "_deriv_" + deriv;
fo.Directory = "data";
so.vectorType = SaveOptions::simple;
so.extras = SaveOptions::none;

save(fo,so,Y);
so.vectorType = SaveOptions::append;
so.printMessage = false;
save(fo,so,dXdY_Z);
save(fo,so,rhs);
save(fo,so,diff);
save(fo,so,constError);

//so.extras = SaveOptions::loc;
//Filename fm = (string)"data/dataMatrix.dat";
//save(fm,so,data);

// plotting vectors
PlotOptions po;
po.style = "linespoints";
Filename plotFile = fo;
plotFile.Directory = "pics";
plotFile.Suffix = ".png";
po.output = plotFile;
po.column = 1;
po.column2 = 2;
po.column3 = 3;
po.printMessage = true;

plot(fo,po);

// printing maxCoeffs
cout << "diff.maxCoeff()       = " << diff.maxCoeff() << endl;
cout << "constError.maxCoeff() = " << constError.maxCoeff() << endl;

return 0;
}

/*-------------------------------------------------------------------------------------------------------------------------
 	program to clean up data
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
	0. declaring some functions that is will use
	1. defining key parameters
	2. argv inputs
	3. getting data from file
	4. removing rouges
	5. replacing file
	6. functions outside main
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	0. function declarations
-------------------------------------------------------------------------------------------------------------------------*/

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeCol(Eigen::MatrixXd& matrix, unsigned int colToRemove);

int main(int argc, char** argv) {

/*-------------------------------------------------------------------------------------------------------------------------
	1. defining key parameters

n.b. the first column is column 0
-------------------------------------------------------------------------------------------------------------------------*/

Filename f = ("results/mainResults2.dat");
mat data;
uint colError = 13;
double closeness = 1.0e-5;

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
		else if (id.compare("col")==0 || id.compare("column")==0) colError = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("close")==0 || id.compare("closeness")==0) closeness = stringToNumber<double>(argv[2*j+2]);
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

/*-------------------------------------------------------------------------------------------------------------------------
	3. getting data from file
-------------------------------------------------------------------------------------------------------------------------*/

SaveOptions so;
so.extras = SaveOptions::none;
so.printType = SaveOptions::ascii;
so.printMessage = true;

load(f,so,data);

/*-------------------------------------------------------------------------------------------------------------------------
	4. removing rogues
		- defining data by column
		- possible derivatives:
			- dWdE|N = -T
			- dWdN|E = -theta
			- dN/dE|W = -T/theta
			- 2dS/dT|theta = E
			- 2dS/dtheta|T = N
	
	n.b. T = 2Tb
-------------------------------------------------------------------------------------------------------------------------*/

vec error(data.rows());
vector<uint> toRemove;

cout << "removing rows: " << endl;
for (uint j=0; j<data.rows(); j++) {
	error(j) = data(j,colError);
	if (abs(error(j))>closeness) {
		toRemove.push_back(j);
		cout << j << endl;
	}
}

for (uint j=0; j<toRemove.size(); j++) {
	removeRow(data,toRemove[toRemove.size()-1-j]);
}

	
/*-------------------------------------------------------------------------------------------------------------------------
	5. printing results
		- saving vectors to file
		- plotting dXdY_Z vs rhs
		- printing maxCoeffs
-------------------------------------------------------------------------------------------------------------------------*/

// saving vectors to file
Filename fo = f;
fo.ID = f.ID + "fix";
save(fo,so,data);

return 0;
}

// removeRow
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

// removeCol
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

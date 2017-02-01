/*-------------------------------------------------------------------------------------------------------------------------
 	program to simplify data from conventions of print.cc to those of print3.cc
 -------------------------------------------------------------------------------------------------------------------------*/

#include <string>
#include <fstream>
#include <iostream>
#include "simple.h"
#include "folder.h"
#include "parameters.h"
#include "print.h"
#include "print3.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	0. enums
	1. defining key parameters
	2. getting argv
	3. loading and saving
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	0. enums
-------------------------------------------------------------------------------------------------------------------------*/

struct DataType {
	enum Option { vector, matrix, sparse};
};

/*-------------------------------------------------------------------------------------------------------------------------
	1. defining key parameters
-------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {

string fi, fo;
string dataType;
bool ascii = false;

/*-------------------------------------------------------------------------------------------------------------------------
	2. getting argv
-------------------------------------------------------------------------------------------------------------------------*/

// getting argv inputs
if (argc % 2 && argc>1) {
	for (int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("f")==0 || id.compare("fi")==0 || id.compare("fileIn")==0) 	fi = (argv[2*j+2]);
		else if (id.compare("fo")==0 || id.compare("fileOut")==0) 					fo = (argv[2*j+2]);
		else if (id.compare("type")==0 || id.compare("dataType")==0) 				dataType = (argv[2*j+2]);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

if (fi.empty() || fo.empty()) {
	cerr << "must give a input file and output file to simplifyData" << endl;
	return 1;
}

DataType::Option dt = DataType::vector;
if (!dataType.empty()) {
	if (dataType.compare("vector")==0)
		dt = DataType::vector;
	else if (dataType.compare("mat")==0 || printOpts.compare("matrix")==0)
		dt = DataType::matrix
	else if (dataType.compare("spMat")==0 || printOpts.compare("sparse")==0)
		dt = DataType::sparse;
	else
		cerr << "data type not understood: " << dataType << endl;
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. loading and saving
-------------------------------------------------------------------------------------------------------------------------*/

if (ascii) {
	switch(dt) {
		case DataType::vector:		vec v;
									load(fi,
									
									
									
									break;
		case DataType::matrix:		col = 2;
									break;
		case DataType::sparse:		col = 4;
									break;
		default:					cerr << "data type not understood: " << dataType << endl;
									break;
	}
}
else {
	switch(dt) {
		case DataType::vector:		SaveOptions so_tp;
									so_tp.printType = SaveOptions::binary;
									so_tp.paramsIn = psu; // AH! WHAT IS ALL THIS ABOUT? WHAT IF I DON'T HAVE THE PARAMETERS?
									so_tp.paramsOut = ps;
									so_tp.vectorType = SaveOptions::complex;
									so_tp.extras = SaveOptions::coords;
									so_tp.zeroModes = 2;
									so_tp.printMessage = true;
									
									
									
									
									break;
		case DataType::matrix:		
									
									
									
									
									break;
		case DataType::sparse:		
									
									
									
									break;
		default:					cerr << "data type not understood: " << dataType << endl;
									break;
	}
}

		

return 0;
}

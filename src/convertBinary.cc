/*
	convertBinary
		program to convert ascii files to binary
*/

#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include "folder.h"
#include "print.h"
#include "simple.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - defining main quantities
		2 - getting argv inputs
		3 - deciding what to save and saving
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/* -------------------------------------------------------------------------------------------------------------------------
	1. defining main quantities
-------------------------------------------------------------------------------------------------------------------------*/
Filename fin, fout;
SaveOptions opts;
opts.printType = SaveOptions::ascii;

/* -------------------------------------------------------------------------------------------------------------------------
	2. getting argv inputs
-------------------------------------------------------------------------------------------------------------------------*/
if (argc==2) {
	fin = (string)argv[1];
}
else {
	cerr << "inputs not understood" << endl;
	return 1;
}

if ((fin.Suffix).compare(".dat")!=0) {
	cerr << "should input .dat file" << endl;
	return 1;
}
	
fout = fin;
fout.Suffix = ".data";

/* -------------------------------------------------------------------------------------------------------------------------
	3. deciding what to save and saving
-------------------------------------------------------------------------------------------------------------------------*/

string id = fin.ID;
bool isStatic = false;
bool isMain = false;
if (id.back()=='E') id = id.substr(0,id.size()-1);
if ((id.substr(0,4)).compare("main")==0) {
	isMain = true;
	id = id.substr(4);
}
if ((id.substr(0,6)).compare("static")==0) {
	id = id.substr(6);
	isStatic = true;
}

if (id.compare("DDS")==0) {
	cerr << "cannot yet print sparse matrices in binary" << endl;
	return 1;	
}
if ((id.substr(0,5)).compare("omega")==0 || (id.substr(0,5)).compare("freqs")==0 \
	|| (id.substr(0,5)).compare("modes")==0) {
	mat m;
	load(fin,opts,m);
	opts.printType = SaveOptions::binary;
	save(fout,opts,m);
}
else if (id.compare("p")==0 || id.compare("pi")==0 || id.compare("tp")==0 ||id.compare("minusDS")==0 || id.compare("delta")==0 \
 		|| id.compare("chiX")==0  || id.compare("chiT")==0) {
 		if (!isStatic) {
 			opts.vectorType = ((isMain || id.compare("tp")==0)? SaveOptions::real: SaveOptions::realB);
 			cerr << "cannot continue, cannot get paramsIn and paramsOut retrospectively" << endl;
 			return 1;
			vec v;
			load(fin,opts,v);
			opts.printType = SaveOptions::binary;
			save(fout,opts,v);
		}
}
else if (id.compare("Step")==0) {
	cerr << "cannot convert Step file " << fin << endl;
	return 1;
}
else {
	cerr << "cannot convert given file " << fin << endl;
	return 1;
}


return 0;
}

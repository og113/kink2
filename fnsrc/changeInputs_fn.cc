/*-------------------------------------------------------------------------------------------------------------------------
 	program to change inputs or options, as a function not a main
 -------------------------------------------------------------------------------------------------------------------------*/

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "changeInputs_fn.h"
#include "simple.h"
#include "folder.h"
#include "parameters.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. defining key parameters
	2. getting argv
	3. changing inputs or options
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

int changeInputs_fn(int argc, vector<string> argv) {

/*-------------------------------------------------------------------------------------------------------------------------
	1. defining key parameters
-------------------------------------------------------------------------------------------------------------------------*/

string name, value;
string fi = "inputsP", fo;

/*-------------------------------------------------------------------------------------------------------------------------
	2. getting argv
-------------------------------------------------------------------------------------------------------------------------*/

// getting argv inputs
if (argc % 2 && argc>1) {
	for (int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("f")==0 || id.compare("fi")==0 || id.compare("fileIn")==0) fi = (string)(argv[2*j+2]);
		else if (id.compare("fo")==0 || id.compare("fileOut")==0) fo = (string)(argv[2*j+2]);
		else if (id.compare("n")==0 || id.compare("name")==0) name = argv[2*j+2];
		else if (id.compare("v")==0 || id.compare("value")==0) value = argv[2*j+2];
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

if (name.empty() || value.empty()) {
	cerr << "must give a name and value to changeParams" << endl;
	return 1;
}
if (fo.empty()) fo = fi;

/*-------------------------------------------------------------------------------------------------------------------------
	3. changing inputs or options
-------------------------------------------------------------------------------------------------------------------------*/
if (fi.find("inputs")!=string::npos && fo.find("inputs")!=string::npos) {
	Parameters ps;
	ps.load(fi);
	if (name[0]=='N' || name.compare("pot")==0) {
		ps.changeParameters(name,stringToNumber<uint>(value));
	}
	else {
		ps.changeParameters(name,stringToNumber<double>(value));
	}
	ps.save(fo);
}
else if (fi.find("options")!=string::npos && fi.find("options")!=string::npos) {
	Options opts;
	opts.load(fi);
	if (name.compare("loops")==0) {
		opts.changeOptions(name,stringToNumber<uint>(value));
	}
	else if (name.compare("alpha")==0 || name.compare("open")==0 || name.compare("amp")==0 || name.compare("loopMin")==0 \
			|| name.compare("loopMax")==0 || name.compare("epsiTb")==0 || name.compare("epsiTheta")==0) {
		opts.changeOptions(name,stringToNumber<double>(value));
	}		
	else {
		opts.changeOptions(name,value);
	}
	opts.save(fo);
}

return 0;
}

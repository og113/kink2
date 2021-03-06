/*----------------------------------------------------------------------------------------------------------------------------
	plot
		program to plot things using gnuplot
----------------------------------------------------------------------------------------------------------------------------*/
#include <Eigen/Sparse>
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
#include "print.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - defining main quantities
		2 - getting argv inputs
		3 - deciding plot options
		4 - plotting
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/* -------------------------------------------------------------------------------------------------------------------------
	1. defining main quantities
-------------------------------------------------------------------------------------------------------------------------*/
Filename f;

/* -------------------------------------------------------------------------------------------------------------------------
	2. getting argv inputs
-------------------------------------------------------------------------------------------------------------------------*/
if (argc==2) {
	f = (string)argv[1];
}
else
	cerr << "inputs not understood" << endl;

/* -------------------------------------------------------------------------------------------------------------------------
	3. deciding plot options
-------------------------------------------------------------------------------------------------------------------------*/
if ((f.Suffix).compare(".data")==0) {
	SaveOptions so;
	so.printType = SaveOptions::binary;
	so.printMessage = false;
	vec v;
	load(f,so,v);
	so.printType = SaveOptions::ascii;
	f.Suffix = ".dat";
	save(f,so,v);
}

PlotOptions po;
po.printMessage = true;
po.output = "gui";
string id = f.ID;
po.column2 = 0;
po.column3 = 0;
bool isStatic = false;
if (id.back()=='E') id = id.substr(0,id.size()-1);
if ((id.substr(0,4)).compare("main")==0) {
	id = id.substr(4);
}
if ((id.substr(0,6)).compare("static")==0) {
	id = id.substr(6);
	isStatic = true;
}

if (id.compare("DDS")==0 || (id.substr(0,5)).compare("omega")==0 \
	|| (id.substr(0,5)).compare("modes")==0) {
	cerr << "cannot plot matrices" << endl;
	return 1;	
}
else if (id.compare("p")==0 || id.compare("pi")==0 || id.compare("minusDS")==0 || id.compare("delta")==0 \
 		|| id.compare("chiX")==0  || id.compare("chiT")==0) {
	if (!isStatic) {
		po.gp = "gp/repi.gp";
		po.style = "points";
	}
	else {
		po.column = 1;
	}
}
else if (id.compare("Step")==0) {
	po.column = 8;
	po.column2 = 9;
	po.style = "points";
}
else {
	po.column = 1;
}

/* -------------------------------------------------------------------------------------------------------------------------
	4. plotting
-------------------------------------------------------------------------------------------------------------------------*/
plot(f,po);

return 0;
}

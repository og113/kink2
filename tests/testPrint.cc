/*
	main for program to test print.cc
*/

#include <iostream>
#include <Eigen/Dense>
#include "folder.h"
#include "parameters.h"
#include "print.h"
#include "simple.h"

using namespace std;

int main() {
cout << "test print: " << endl;

Parameters p;
p.load("inputsM");

SaveOptions so;
so.vectorType = SaveOptions::complex;
so.paramsIn = p;
so.paramsOut = p;
so.printMessage = false;
so.zeroModes = 2;


Filename f, g;
f = "data/000mainp_fLoop_0_loop_0.dat";
g = "data/000mainp_fLoop_0_loop_0.data";

vec u, v, w;

uint N = 10;
for (uint j=0; j<N; j++) {
	so.printType = SaveOptions::ascii;
	load(f,so,u);

	so.printType = SaveOptions::binary;
	save(g,so,u);
	load(g,so,v);
}

cout << "u.size() = " << u.size() << endl;
cout << "v.size() = " << v.size() << endl;
cout << "u(0)     = " << u(0) << endl;
cout << "v(0)     = " << v(0) << endl;
cout << "u.norm() = " << u.norm() << endl;
cout << "v.norm() = " << v.norm() << endl;

w = u-v;

cout << "test: " << w.norm() << endl;

return 0;
}

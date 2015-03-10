/*
	program to test thetaT.cc and thetaT.h
*/

#include <iostream>
#include <complex>
#include "potentials.h"
#include "thetaT.h"

using namespace std;

typedef complex<double> comp;

int main() {

params_for_V paramsV;
paramsV.epsi = 0.05;
paramsV.aa = 0.4;

ec_params paramsEC;
paramsEC.aa = 0.4;
paramsEC.minima0 = -1.0;
paramsEC.minima1 = 1.0;
paramsEC.de = 0.05;

//struct void_params paramsVoid;

double d = 0.9;
comp c(-1.0,0.0);

PotentialType pot = &V2<comp>;

Potential classV(paramsV,pot);

cout << classV(c) << endl;

cout << ec(d,&paramsEC) << endl;

cout << rhoIntegrand(d,&paramsV) << endl;

return 0;
}

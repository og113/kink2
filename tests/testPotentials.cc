/*
	main for program to test potentials.cc, potentials.h and fnptrs.h
*/

#include <iostream>
#include "potentials.h"

using namespace std;

int main() {
cout << "test: " << endl;

struct params_for_V paramsV;
paramsV.epsi = 0.05;
paramsV.aa = 0.4;

//struct void_params paramsVoid;

double d = 0.9;
comp c(-1.0,0.0);

//void* v = &paramsV;

//cout << V1(c,v) << " " << endl;

//PotentialType pot = &V2<comp>;

Potential classV;

classV((PotentialType)&V2<comp>,paramsV);

cout << classV(c) << endl;
cout << classV(d) << endl;

return 0;
}

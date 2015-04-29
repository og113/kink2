/*
	main for program to test potentials.cc, potentials.h and fnptrs.h
*/

#include <iostream>
#include <Eigen/Dense>
#include "folder.h"
#include "parameters.h"
#include "potentials.h"
#include "print.h"

using namespace std;

int main() {
cout << "test: " << endl;

struct params_for_V paramsV;
paramsV.epsi = 0.05;
paramsV.aa = 0.4;

//struct void_params paramsVoid;

double d = 0.9;
comp C(-1.0,0.0);

//void* v = &paramsV;

//cout << V1(C,v) << " " << endl;

//PotentialType pot = &V2<comp>;

Potential<comp> classV;

classV((Potential<comp>::PotentialType)&V2<comp>,paramsV);

cout << classV(C) << endl;
cout << classV(d) << endl;

uint N = 1e5;
double x0 = -5.0, x1 = 5.0;
double h = (x1-x0)/(double)(N-1.0);

Parameters ps;
ps.load("inputsP");

Potential<double> V, dV, ddV;
V((Potential<double>::PotentialType)&V2<double>,ps);
dV((Potential<double>::PotentialType)&dV2<double>,ps);
ddV((Potential<double>::PotentialType)&ddV2<double>,ps);

vec a(N), b(N), c(N);

double locMax = x0;
double maxDiff = 0.0;

for (uint j=0; j<N; j++) {
	double x = x0 + (double)j*h;
	a(j) = (V(x+h)+V(x-h)-2.0*V(x))/pow(h,2.0);
	b(j) = (dV(x+h)-dV(x-h))/2.0/h;
	c(j) = ddV(x);
	if (absDiff(a(j),b(j))>maxDiff) {
		maxDiff = absDiff(a(j),b(j));
		locMax = x;
	}
	if (absDiff(a(j),c(j))>maxDiff) {
		maxDiff = absDiff(a(j),c(j));
		locMax = x;
	}
	if (absDiff(b(j),c(j))>maxDiff) {
		maxDiff = absDiff(b(j),c(j));
		locMax = x;
	}	
}

SaveOptions so;
so.paramsIn = ps; so.paramsOut = ps;
so.vectorType = SaveOptions::simple;
so.extras = SaveOptions::loc;
so.printMessage = true;

Filename fileABC = (string)"tests/data/abc.dat";

save(fileABC,so,a);
so.vectorType = SaveOptions::append;
so.extras = SaveOptions::none;
so.printMessage = false;
save(fileABC,so,b);
save(fileABC,so,c);

PlotOptions po;
po.column = 1;
po.column2 = 2;
po.column3 = 3;
po.style = "linespoints";
po.printMessage = true;

Filename plotABC = (string)"tests/data/abc.png";
po.output = plotABC;

plot(fileABC,po);

cout << "stepSize = " << h << endl;
cout << "maxDiff = " << maxDiff << endl;
cout << "locMax = " << locMax << endl;

return 0;
}

/*
	main for program to test lattice.cc
*/

#include <iostream>
#include "lattice.h"

using namespace std;

int main() {
cout << "test lattice: " << endl;

PrimaryParameters p1;

p1.pot = 1;
p1.N = 130;
p1.Na = 300;
p1.Nb = 80;
p1.Nc = 2;
p1.LoR = 3.2;
p1.dE = 0.05;
p1.Tb = 18.0;
p1.theta = 0.0;
p1.reg = 0.0;

cout << p1;

Parameters p(p1);

uint arb = 1030;
cout << "test1: " << c(intCoord(arb,0,p),intCoord(arb,1,p),p)-arb << endl;

cout << "coord(0,0,p) = " << coord(0,0,p) << endl;
cout << "coord(0,1,p) = " << coord(0,1,p) << endl;
cout << "coord(p.N*p.NT-1,0,p) = " << coord(p.N*p.NT-1,0,p) << endl;
cout << "coord(p.N*p.NT-1,1,p) = " << coord(p.N*p.NT-1,1,p) << endl;

cout << "neigh(0,0,1,p) = " << neigh(0,0,1,p) << endl;
cout << "neigh(0,1,1,p) = " << neigh(0,1,1,p) << endl;
cout << "neigh(p.N*p.NT-1,0,1,p) = " << neigh(p.N*p.NT-1,0,1,p) << endl;
cout << "neigh(p.N*p.NT-1,1,1,p) = " << neigh(p.N*p.NT-1,1,1,p) << endl;

cout << "p.a = " << p.a << ", p.b = " << p.b << endl;
cout << "dtFn(0,p) = " << dtFn(0,p) << endl;
cout << "DtFn(0,p) = " << DtFn(0,p) << endl;
cout << "dxFn(0,p) = " << dxFn(0,p) << endl;
cout << "DxFn(0,p) = " << DxFn(0,p) << endl;

cout << "interpolate, vecComplex and vecReal untested" << endl;

return 0;
}

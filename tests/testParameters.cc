/*
	program to test parameters.cc and parameters.h
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "parameters.h"

using namespace std;

int main() {

cout << "nothing tested yet" << endl;

PrimaryParameters p1, P1;

p1.pot = 1;
p1.N = 130;
p1.Na = 300;
p1.Nb = 80;
p1.Nc = 2;
p1.LoR = 3.2;
p1.dE = 0.005;
p1.Tb = 18.0;
p1.theta = 0.0;
p1.reg = 0.0;

string filename = "data/tests/testParameters_p1";
p1.save(filename);
cout << filename << " printed" << endl;

P1.load(filename);
cout << "P1: " << endl << P1;

SecondaryParameters p2;
p2.setSecondaryParameters(P1);

cout << "p2: " << endl << p2;

Parameters p(p1,p2);
cout << "p: " << endl;
p.print();
return 0;
}

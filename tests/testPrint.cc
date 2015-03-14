/*
	main for program to test print.cc
*/

#include <iostream>
#include "print.h"

using namespace std;

int main() {
cout << "test print: " << endl;

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

cout << "save, load and plot untested" << endl;

return 0;
}

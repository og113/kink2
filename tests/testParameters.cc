/*
	program to test parameters.cc and parameters.h
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "parameters.h"

using namespace std;

int main() {

cout << "nothing tested yet" << endl;

PrimaryParameters p1;

p1.pot = 1;
p1.N = 1;
p1.Na = 1;
p1.Nb = 1;
p1.Nc = 1;
p1.LoR = 1.0;
p1.dE = 1.0;
p1.Tb = 1.0;
p1.theta = 1.0;
p1.reg = 1.0;

p1.save("data/tests/testParameters_p1");

return 0;
}

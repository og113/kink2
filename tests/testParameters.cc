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
p1.dE = 0.05;
p1.Tb = 18.0;
p1.theta = 0.0;
p1.reg = 0.0;

string filename = "tests/data/testParameters_p1";
p1.save(filename);
cout << filename << " printed" << endl;

P1.load(filename);
cout << "P1: " << endl << P1;

SecondaryParameters p2;
p2.setSecondaryParameters(P1);

cout << "p2: " << endl << p2;

Parameters p(p1);
//p.load(filename);
p.changeParameters("theta",0.1);
p.changeParameters("Tb",17.9);

cout << "p: " << endl;
p.print();

Options o;
o.alpha = 0.5;
o.open = 1.0;
o.amp = 0.5;
o.zmx = "nC2";
o.zmt = "nC2";
o.bds = "all";
o.inF = "p";
o.minTimenumberLoad = "0";
o.maxTimenumberLoad = "0";
o.minLoopLoad = "0";
o.maxLoopLoad = "0";
o.loopChoice = "n";
o.loopMin = 0.0;
o.loopMax = 0.0;
o.loops = 1;
o.printChoice = "n";

string str = "tests/data/testParameters_o";
o.save(str);
o.print();
cout << "o: " << endl << o;

return 0;
}

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
cout << "o: " << endl << o << endl;

Closenesses c;
c.Step = 1.0e-2;
c.Action = 1.0e-2;
c.Soln = 1.0e-6;
c.SolnMax = 1.0e-5;
c.Delta = 1.0;
c.Inv = 1.0e-16;
c.Con = 1.0e-2;
c.Lin = 5.0e-2;
c.True = 5.0e-2;
c.Latt = 2.0e-1;
c.Reg = 1.0e-2;
c.IE = 1.0e-5;
c.Contm = 5.0e-2;
c.OS = 5.0e-2;
c.AB = 5.0e-2;
c.ABNE = 5.0e-2;
c.LR = 1.0e-12;
c.DT = 1.0e3;
c.Profile = 1.0e-5;

str = "tests/data/testParameters_c";
c.save(str);
cout << "c: " << endl << c;

return 0;
}

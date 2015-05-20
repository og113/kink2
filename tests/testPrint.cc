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

//Parameters p;
//p.load("inputsM");

SaveOptions so;
so.printType = SaveOptions::binary;
so.vectorType = SaveOptions::simple;
//so.paramsIn = p;
//so.paramsOut = p;
so.printMessage = false;

clock_t time;
time = clock();

cout << "testing vector save and load" << endl;

Filename f, g;
f = "tests/data/u.data";
g = "tests/data/v.data";

vec u, v, w;
uint size = 1000;
u = Eigen::VectorXd::Random(size);

uint N = 10;
for (uint j=0; j<N; j++) {
	save(f,so,u);
	load(f,so,v);
	save(g,so,v);
	load(g,so,u);
}

cout << "u.size()          = " << u.size() << endl;
cout << "v.size()          = " << v.size() << endl;
cout << "u(u.size()-1)     = " << u(u.size()-1) << endl;
cout << "v(v.size()-1)     = " << v(v.size()-1) << endl;
cout << "u.norm()          = " << u.norm() << endl;
cout << "v.norm()          = " << v.norm() << endl;

w = u-v;

cout << "test vector: " << w.norm() << endl;

cout << "testing matrix save and load" << endl;

f = "tests/data/m.data";
g = "tests/data/n.data";

mat m, n, o;
m = Eigen::MatrixXd::Random(size,size);

for (uint j=0; j<N; j++) {
	save(f,so,m);
	load(f,so,n);
	save(g,so,n);
	load(g,so,m);
}

cout << "m.size() = " << m.size() << endl;
cout << "n.size() = " << n.size() << endl;
cout << "m(0,0)   = " << m(0) << endl;
cout << "n(0,0)   = " << n(0) << endl;
cout << "m.norm() = " << m.norm() << endl;
cout << "n.norm() = " << n.norm() << endl;

o = m-n;

cout << "test matrix: " << o.norm() << endl;

time = clock() - time;
double realtime = time/1000000.0;

cout << endl << N << " load and save loops took " << realtime << " seconds" << endl << endl;

return 0;
}

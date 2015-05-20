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

Parameters p;
p.load("inputsM");

SaveOptions so;
so.printType = SaveOptions::ascii;
so.vectorType = SaveOptions::complex;
so.paramsIn = p;
so.paramsOut = p;
so.printMessage = false;
so.zeroModes = 2;

string d("data/001mainp_fLoop_0_loop_0.data");
string e((string)"temp/001mainp_fLoop_0_loop_0.data");
string f("data/000tp_fLoop_0_loop_0.data");
string g("temp/000tp_fLoop_0_loop_0.data");
vec r, s, t;
vec u, v, w;

clock_t time;
time = clock();

uint loops = 1;
for (uint j=0; j<loops; j++) {
	
	so.printType = SaveOptions::binary;
	load(d,so,r);
	load(f,so,u);

	/*so.printType = SaveOptions::binary;
	save(e,so,r);
	save(g,so,u);

	load(e,so,s);
	load(g,so,v);*/
}

time = clock() - time;
double realtime = time/1000000.0;

cout << endl << loops << " load and save loops took " << realtime << " seconds" << endl << endl;

/*cout << "r.size()          = " << r.size() << endl;
cout << "s.size()          = " << s.size() << endl;
cout << "r(r.size()-1)     = " << r(r.size()-1) << endl;
cout << "s(s.size()-1)     = " << s(s.size()-1) << endl;
cout << "r.norm()          = " << r.norm() << endl;
cout << "s.norm()          = " << s.norm() << endl;
cout << "u.size()          = " << u.size() << endl;
cout << "v.size()          = " << v.size() << endl;
cout << "u(u.size()-1)     = " << u(u.size()-1) << endl;
cout << "v(v.size()-1)     = " << v(v.size()-1) << endl;
cout << "u.norm()          = " << u.norm() << endl;
cout << "v.norm()          = " << v.norm() << endl;
t = r-s;
w = u-v;
cout << "test 1            = " << t.norm() << endl;
cout << "test 2            = " << w.norm() << endl;*/

return 0;
}

/*
	main for program to test stepper.cc
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include "stepper.h"
#include "simple.h"
#include "print.h"
#include "folder.h"

using namespace std;

double f(const Point2d& p) {
	return pow(p.X,2.0)+pow(p.Y,2.0);
}


int main() {
cout << "test stepper: " << endl;
cout << "pi = " << pi << endl;
cout << "MIN_NUMBER = " << MIN_NUMBER << endl;

StepperOptions sto;
sto.stepType = StepperOptions::constSimple;
sto.directed = StepperOptions::local;
sto.closeness = 1.0e-2;
sto.angle0 = 0.0;
sto.epsi_x = 0.05;
sto.epsi_y = 0.1;

Point2d P(5.0,0.0);
Point2d Q(5.0,1.0);
Point2d R(6.0,1.0);
Point2d S(7.0,1.0);
Point2d T(7.0,2.0);

vector<FxyPair> fxy;
fxy.push_back(FxyPair(Q,f(Q)));
fxy.push_back(FxyPair(R,f(R)));
fxy.push_back(FxyPair(S,f(S)));
fxy.push_back(FxyPair(T,f(T)));
fxy.push_back(FxyPair(P,f(P)));

cout << "1st closest = " << find_nth_closest(fxy, f(P), 1) << endl;
cout << "2nd closest = " << find_nth_closest(fxy, f(P), 2) << endl;
cout << "3rd closest = " << find_nth_closest(fxy, f(P), 3) << endl;
cout << "4th closest = " << find_nth_closest(fxy, f(P), 4) << endl;

cout << "initializing stepper: ";
Stepper st(sto,P);
st.addResult(f(P));
cout << "done" << endl;

string file = "tests/data/testStepper.dat";
ofstream os;
os.open(file.c_str());
cout << "stepping variables:" << endl;
cout << P << setw(15) << f(P) << endl;
for (unsigned int j=0; j<1000; j++) {
	st.step();
	P = st.point();
	//cout << P << setw(15) << f(P) << endl;
	os << P << setw(15) << f(P) << endl;
	st.addResult(f(P));
	if (j%100==0)
		cout << st.stepAngle() << endl;
}
os.close();

cout << "plotting steps" << endl;
PlotOptions po;
po.column = 1;
po.column2 = 2;
po.style = "points";
po.output = "tests/data/testStepper.png";
po.printMessage = true;

plot(file,po);

return 0;
}

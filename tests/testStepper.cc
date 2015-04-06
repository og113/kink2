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
sto.angle0 = pi/2.0;
sto.epsi_x = 0.01;
sto.epsi_y = 0.02;

Point2d P(5.0,0.0);
double maxError = 0.0;
double f0 = f(P);

cout << "initializing stepper: ";
Stepper st(sto,P);
st.addResult(f(P));
cout << "done" << endl;

string file = "tests/data/testStepper.dat";
ofstream os;
os.open(file.c_str());
cout << "stepping variables:" << endl;
for (unsigned int j=0; j<1e4; j++) {
	if (st.keep()) {
		os << P << setw(15) << f(P) << endl;
		//cout << st.stepAngle() << endl;
		maxError = (absDiff(f(P),f0)>maxError? absDiff(f(P),f0): maxError);
	}
	st.step();
	P = st.point();
	st.addResult(f(P));
}
os.close();

cout << "plotting steps" << endl;
PlotOptions po;
po.column = 1;
po.column2 = 2;
po.style = "linespoints";
po.output = "tests/data/testStepper1.png";
po.printMessage = true;

plot(file,po);

cout << "maxError = " << maxError << endl;
cout << "closeness = " << sto.closeness << endl;

return 0;
}

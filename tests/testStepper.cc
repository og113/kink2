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
sto.constant = true;
sto.angle = 0.0;
sto.epsi_x = 0.05;
sto.epsi_y = 0.1;

Point2d P(5.0,0.0);

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
	st.stepVariables();
	P = st.point();
	//cout << P << setw(15) << f(P) << endl;
	os << P << setw(15) << f(P) << endl;
	st.addResult(f(P));
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

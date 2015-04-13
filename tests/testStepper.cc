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
	//return pow(p.X,2.0)+pow(p.Y,2.0);
	//return pow(p.X,-2.0)+pow(p.Y,-2.0);
	//return p.X*exp(-pow(p.X,2.0)-pow(p.Y,2.0));
	return exp(+pow(p.X,2.0)+pow(p.Y,2.0))/p.X/p.Y;
}


int main() {
cout << "test stepper: " << endl;
string stepType;
cout << "input either t (constTaylor) or p (constPlane): ";
cin >> stepType;

uint 			loops 				= 1e1;
uint			avgLoops			= 1e0;
uint 			parameterLoops 		= 1;
double			closenessMin 		= 0.005;
double 			closenessMax 		= 0.025;

cout << "loops          = " << loops << endl;
cout << "avgLoops       = " << avgLoops << endl;
cout << "parameterLoops = " << parameterLoops << endl;
cout << "closenessMin   = " << closenessMin << endl;
cout << "closenessMax   = " << closenessMax << endl;
cout << endl;

StepperOptions closenessStepperOpts;
closenessStepperOpts.stepType = StepperOptions::straight;
closenessStepperOpts.directed = StepperOptions::undirected;
closenessStepperOpts.closeness = 1.0;
closenessStepperOpts.angle0 = 0.0;
closenessStepperOpts.epsi_x = (closenessMax-closenessMin)/(parameterLoops-1.0);
closenessStepperOpts.epsi_y = 0.0;
Stepper closenessStepper(closenessStepperOpts);

printf("%20s%20s%20s%20s\n","closeness","avgLocal","minLocal","maxLocal");

for (uint l=0; l<parameterLoops; l++) {
	if (l==0) {
		Point2d closenessStart(closenessMin,0.0);
		closenessStepper.setStart(closenessStart);
	}
	else {
		closenessStepper.step();
	}
	
	double avgLocal = 0.0;
	double minLocal = 1.0e16;
	double maxLocal = 0.0;
	for (uint m=0; m<avgLoops; m++) {
		StepperOptions sto;
		if (stepType.compare("t")==0)
			sto.stepType = StepperOptions::constTaylor;
		else if (stepType.compare("p")==0)
			sto.stepType = StepperOptions::constPlane;
		else {
			cerr << "stepType, " << stepType << ", not understood" << endl;
			return 1;
		}
		sto.directed = StepperOptions::undirected;
		sto.closeness = closenessStepper.x();
		sto.angle0 = pi/2.0;
		sto.epsi_x = 0.01;
		sto.epsi_y = 0.02;

		Point2d P(1.0,1.0);
		double maxError = 0.0;
		double f0 = f(P);

		//cout << "initializing stepper: ";
		Stepper st(sto,P);
		st.addResult(f(P));
		//cout << "done" << endl;

		string file = "tests/data/testStepper.dat";
		ofstream os;
		os.open(file.c_str());
		//cout << "stepping variables:" << endl;
		//cout << "             loops: " << loops << endl;
		bool coutEveryLoop = false;
		bool coutAngle = false;
		if (coutEveryLoop) {
			cout << P << setw(15) << 0.0 << setw(15) << f(P) << setw(15) << "y" << endl;
		}
		for (unsigned int j=0; j<loops; j++) {
			if (st.keep()) {
				os << P << setw(15) << f(P) << endl;
				if (!coutEveryLoop && coutAngle) cout << st.stepAngle() << endl;
				maxError = (absDiff(f(P),f0)>maxError? absDiff(f(P),f0): maxError);
			}
			st.step();
			if (coutEveryLoop) {
			string keep = (st.keep()? "y": "n");
			cout << P << setw(15) << st.stepAngle() << setw(15) << f(P);
			}
			P = st.point();
			st.addResult(f(P));
			if (coutEveryLoop) {
			string keep = (st.keep()? "y": "n");
			cout << setw(15) << keep << endl;
			}
		}
		os.close();

		/*
		cout << "plotting steps" << endl;
		PlotOptions po;
		po.column = 1;
		po.column2 = 2;
		po.style = "linespoints";
		po.output = "tests/data/testStepper.png";
		po.printMessage = true;

		plot(file,po);

		cout << "steps = " << st.steps() << endl;
		cout << "avgLocal = " << (double)loops/(double)st.steps() << endl;
		cout << "maxError = " << maxError << endl;
		cout << "closeness = " << sto.closeness << endl;
		*/
		if (st.steps()>0)
			avgLocal += (double)loops/(double)st.steps();
		else
			avgLocal += (double)loops;
		if ((double)loops/(double)st.steps()<minLocal)
			minLocal = (double)loops/(double)st.steps();
		if ((double)loops/(double)st.steps()>maxLocal)
			maxLocal = (double)loops/(double)st.steps();
	}
	avgLocal /= (double)avgLoops;
	closenessStepper.addResult(1.0);
	printf("%20g%20g%20g%20g\n",closenessStepper.x(),avgLocal,minLocal,maxLocal);
}

return 0;
}

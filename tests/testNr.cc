/*
	main for program to test nr.cc
*/

#include <iostream>
#include "nr.h"

using namespace std;

int main() {
cout << "test nr: " << endl;

// parameters
PrimaryParameters pr1;
pr1.pot=3;
pr1.N = 2;
pr1.Na = 0;
pr1.Nb = 2;
pr1.Nc = 0;
pr1.LoR = 1.0;
pr1.dE = 0.01;
pr1.Tb = 1.0;
pr1.theta = 0.0;
pr1.reg = 0.01;
Parameters pr(pr1);

// scalars
uint length = pr.N*pr.NT;
comp kinetic = 0.0;

// vectors
vec p, mds;
cVec f;
p = Eigen::VectorXd::Random(2*length);
mds = Eigen::VectorXd::Zero(2*length);
f = Eigen::VectorXcd::Random(length);

// matrix
spMat dds(2*length,2*length);

for (uint i=0; i<length; i++) {

	Kinetic_nr(i,0,p,pr,f,kinetic);
	Kinetic_nr(i,0,p,pr,f,kinetic);
	
	mdKinetic_nr(i,0,p,pr,f,mds);
	mdKinetic_nr(i,1,p,pr,f,mds);
	
	for (uint j=0; j<length; j++) {
		ddKinetic_nr(i,j,0,p,pr,f,dds);
		ddKinetic_nr(i,j,1,p,pr,f,dds);
	}
}

cout << "nothing much done" << endl;

return 0;
}

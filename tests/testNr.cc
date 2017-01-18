/*
	main for program to test nr.cc
*/

#include <iostream>
#include "lattice.h"
#include "nr.h"

using namespace std;

int main() {
cout << "test nr: " << endl;

// parameters
PrimaryParameters pr1;
pr1.pot=3;
pr1.N = 4;
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
comp action = 0.0;

// vectors
vec p, mds;
cVec f;
p = Eigen::VectorXd::Constant(2*length,0.0);
p(4) += 1.0e-4;
mds = Eigen::VectorXd::Zero(2*length);
f = Eigen::VectorXcd::Random(length);

// matrix
spMat dds(2*length,2*length);

// potential
Potential<comp> V, dV, ddV;
V((Potential<comp>::PotentialType)&V3<comp>,pr);
dV((Potential<comp>::PotentialType)&dV3<comp>,pr);
ddV((Potential<comp>::PotentialType)&ddV3<comp>,pr);

// assigning preliminary parameter structs
params_for_V paramsV  = {pr.epsilon, pr.A};

for (uint i=0; i<length; i++) {
	//uint t= intCoord(i,0,pr); //coordinates
	uint x = intCoord(i,1,pr);
			
	paramsV.epsi = pr.r0+x*pr.a;
	V.setParams(paramsV);
	dV.setParams(paramsV);
	ddV.setParams(paramsV);

	Kinetic_nr(i,0,p,pr,f,action);
	Kinetic_nr(i,0,p,pr,f,action);
	Potential_nr(i,p,pr,V,f,action);
	
	mdKinetic_nr(i,0,p,pr,f,mds);
	mdKinetic_nr(i,1,p,pr,f,mds);
	mdPotential_nr(i,p,pr,dV,f,mds);
	
	for (uint j=0; j<length; j++) {
		ddKinetic_nr(i,j,0,p,pr,f,dds);
		ddKinetic_nr(i,j,1,p,pr,f,dds);
		ddPotential_nr(i,j,p,pr,ddV,f,dds);
	}
}

cout << "p = ";
for (uint j=0; j<length; j++) {
	cout << p(j);
	if (j<(length-1))
		cout << ", ";
}
cout << endl;

cout << "mds = ";
for (uint j=0; j<length; j++) {
	cout << mds(j);
	if (j<(length-1))
		cout << ", ";
}
cout << endl;

cout << "action = " << action << endl;

cout << "nothing much done" << endl;

return 0;
}

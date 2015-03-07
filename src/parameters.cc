/*-------------------------------------------------------------------------------------------------------------------------
	definitions for functions and classes for dealing with parameters
-------------------------------------------------------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <cmath>
#include "simple.h"
#include "potentials.h"
#include "parameters.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - errors
	2 - PrimaryParameters
	3 - SecondaryParameters
	4 - Parameters
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1.parameter related errors
		- unset
		- load
-------------------------------------------------------------------------------------------------------------------------*/

// unset
string ParameterError::Unset::message() const{
	return "Parameter " + pName + " not set.";
}

// load
string ParameterError::Load::message() const{
	return "Error: parameters not loaded from " + filename + ".";
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. PrimaryParameters member functions
		- save
		- load
-------------------------------------------------------------------------------------------------------------------------*/

//save
void PrimaryParameters::save(const string& filename) const {
	ofstream os;
	os.open(filename.c_str());
	if (!os.good()) {
		FileError::StreamNotGood e(filename);
		throw e;
	}
	os << left;
	os << setw(20) << "pot" << setw(20) << pot << endl;
	os << setw(20) << "N" << setw(20) << N << endl;
	os << setw(20) << "Na" << setw(20) << Na << endl;
	os << setw(20) << "Nb" << setw(20) << Nb << endl;
	os << setw(20) << "Nc" << setw(20) << Nc << endl;
	os << setw(20) << "LoR" << setw(20) << LoR << endl;
	os << setw(20) << "dE" << setw(20) << dE << endl;
	os << setw(20) << "Tb" << setw(20) << Tb << endl;
	os << setw(20) << "theta" << setw(20) << theta << endl;
	os << setw(20) << "reg" << setw(20) << reg << endl;
	os << endl;
	os.close();
}

//load
void PrimaryParameters::load(const string& filename) {
	ifstream is;
	is.open(filename.c_str());
	if (!is.good()) {
		FileError::StreamNotGood e(filename);
		throw e;
	}
	string line;
	bool set = false;
	while(getline(is,line)) {
		if (line[0]=="#") continue;
		if (line.empty()) break;
		set = true;
		ss >> "pot" >> pot;
		ss >> "N" >> N;
		ss >> "Na" >> Na;
		ss >> "Nb" >> Nb;
		ss >> "Nc" >> Nc;
		ss >> "LoR" >> LoR;
		ss >> "dE" >> dE;
		ss >> "Tb" >> Tb;
		ss >> "theta" >> theta;
		ss >> "reg" >> reg;
		ss.close();
	}
	if (!set) {
		ParameterError::Load e(filename);
		throw e;
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	3.SecondaryParameters member functions
		- set SecondaryParameters from PrimaryParameters
-------------------------------------------------------------------------------------------------------------------------*/

// set secondary parameters
void SecondaryParameters::setSecondaryParameters (const struct PrimaryParameters& pp) {
	NT = Na + Nb + Nc;
	A = 0.4;
	Gamma = exp(-pp.theta);
	
	#define WRAP(FN, TYPE) double FN##_wrapped(double x, void* parameters) { return FN(x, (TYPE*)parameters); }
	WRAP(float_function, float)
	WRAP(double_function, double)
	
	double epsilon0;
	//potential functions
	if (pp.pot==1) {
		Vd = &V1;
		dVd = &dV1;
		ddVd = &ddV1;
		epsilon0 = 0.0;
		epsilon = pp.dE;
		r0 = 0.0;
	}
	else if (pp.pot==2) {
		Vd = &V2;
		dVd = &dV2;
		ddVd = &ddV2;
		epsilon0 = 0.74507774287199924;
		epsilon = 0.75;
		r0 = 0.0;
	}
	else if (pp.pot==3) {
		Vd = &V3;
		dVd = &dV3;
		ddVd = &ddV3;
		epsilon0 = 0.0;
		epsilon = 0.0;
		r0 = MIN_NUMBER;
	}
	else
		cerr << "pot option not available, pot = " << pot << endl;
	paramsV  = {epsilon, A};
	paramsV0 = {epsilon0, A};
	paramsVoid = {};
	
	// epsilon, minima, mass2, action_0
	if (pp.pot!=3) {
		//gsl function for dV(phi)
		gsl_function F;
		F.function = Vd;
		F.params = &paramsV;	

		//finding preliminary roots of dV(phi)=0
		minima[0] = brentMinimum(&F, -1.0, -3.0, 0.0);
		minima[1] = brentMinimum(&F, 1.2, 0.5, 3.0);

		//gsl function for V(root2)-V(root1)-dE
		struct ec_params ec_params = { A, minima[0], minima[1], dE};
		gsl_function EC;
		EC.function = &ec;
		EC.params = &ec_params;

		//evaluating epsilon, new root and dE may change slightly
		epsilonFn(&F,&EC,&dE,&epsilon,&minima);

		//evaluating mass about false vac
		mass2 = ddVd(minima[0],&paramsV);

		//finding root0 of dV0(phi)=0;
		vector<double> minima0(3);
		if (pot[0]=='1')
			{
			minima0[0] = -1.0; minima0[1] = 1.0;
			}
		else if (pot[0]=='2')
			{
			gsl_function V0;
			V0.function = Vd;
			V0.params = &paramsV0;	
			minima0[0] = brentMinimum(&V0, -1.0, -3.0, 0.0);
			minima0[0] = brentMinimum(&V0, 1.2, 0.5, 3.0);
			struct ec_params ec0_params = { A, minima0[0], minima0[1], 0.0};
			gsl_function EC0;
			EC0.function = &ec;
			EC0.params = &ec0_params;
			double dE0 = 0.0;
			epsilonFn(&V0,&EC0,&dE0,&epsilon0,&minima0);
			}

		//finding S1
		double S1error;
		gsl_function S1_integrand;
		S1_integrand.function = &s1Integrand;
		S1_integrand.params = &paramsV0;
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
		gsl_integration_qag(&S1_integrand, minima0[0], minima0[1], MIN_NUMBER, 1.0e-8, 1e4, 4, w, &S1, &S1error);
		gsl_integration_workspace_free(w);
		if (S1error>1.0e-8) { cout << "S1 error = " << S1error << endl;}
	}
	else {
		mass2 = 1.0;
		minima[0] = 0.0; //only one minimum
		R = 10.0; alpha = 0.0; //not used
	}
	
	// R, L, a, b, Ta, Tc

}

/*-------------------------------------------------------------------------------------------------------------------------
	4. Parameters member functions
		- empty constructor
		- constructor from PrimaryParameters and SecondaryParameters
		- print to shell
		- update secondary parameters
		- change parameters based on change in one
-------------------------------------------------------------------------------------------------------------------------*/

// empty constructor
Parameters::Parameters(): PrimaryParmaeters(), SecondaryParameters()	{}			

// constructor using primary and secondary parameters
Parameters::Parameters(const PrimaryParameters& p1, const SecondaryParmeters& p2): \
					PrimaryParmaeters(p1), SecondaryParameters(p2)	{}				

// print to shell
void Parameters::printParameters() const {
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","N","Na","Nb","Nc","L","Ta","Tb","R","dE","theta","reg");
	printf("%8i%8i%8i%8i%8g%8g%8g%8g%8g%8g%8g\n",\
			N,Na,Nb,Nc,\
			L,Ta,Tb,R,\
			dE,theta,reg);
	printf("\n");
}

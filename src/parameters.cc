/*-------------------------------------------------------------------------------------------------------------------------
	definitions for functions and classes for dealing with parameters
-------------------------------------------------------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "simple.h"
#include "potentials.h"
#include "fnptrs.h"
#include "gsl_extras.h"
#include "thetaT.h"
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
		if (line[0]=='#') continue;
		if (line.empty()) break;
		set = true;
		istringstream ss(line);
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
		- wrapped functions to take void* instead of a parameter struct, declared and defined here
		- set SecondaryParameters from PrimaryParameters
-------------------------------------------------------------------------------------------------------------------------*/

// wrapped functions
#define WRAP(FN) double FN##_wrapped(double x, void* parameters) { return FN(x, (params_for_V*)parameters); }
WRAP(Vd)

// set secondary parameters
void SecondaryParameters::setSecondaryParameters (const struct PrimaryParameters& pp) {
	NT = Na + Nb + Nc; 										////////// NT
	A = 0.4;												////////// A
	Gamma = exp(-pp.theta);									////////// Gamma
	
	double epsilon0;
	params_for_V paramsV, paramsV0;
	minima.reserve(2);
	//potential functions
	if (pp.pot==1) {
		Vd = &V1;
		dVd = &dV1;
		ddVd = &ddV1;
		epsilon0 = 0.0;
		epsilon = pp.dE;
		r0 = 0.0;											////////// r0
	}
	else if (pp.pot==2) {
		Vd = &V2;
		dVd = &dV2;
		ddVd = &ddV2;
		epsilon0 = 0.74507774287199924;
		epsilon = 0.75;
		r0 = 0.0;											////////// r0
	}
	else if (pp.pot==3) {
		Vd = &V3;
		dVd = &dV3;
		ddVd = &ddV3;
		epsilon0 = 0.0;
		epsilon = 0.0;
		r0 = MIN_NUMBER;									////////// r0
	}
	else
		cerr << "pot option not available, pot = " << pp.pot << endl;
	paramsV  = {epsilon, A};
	paramsV0 = {epsilon0, A};
	
	// epsilon, minima, mass2, action_0
	if (pp.pot!=3) {
		//gsl function for V(phi)
		gsl_function F;
		F.function = Vd_wrapped;
		F.params = &paramsV;	

		//finding preliminary minima of V(phi)
		minima[0] = brentMinimum(&F, -1.0, -3.0, 0.0);		////////// minima
		minima[1] = brentMinimum(&F, 1.2, 0.5, 3.0);		////////// minima

		//gsl function for V(root2)-V(root1)-dE
		struct ec_params ec_params = { A, minima[0], minima[1], pp.dE};
		gsl_function EC;
		EC.function = &ec;
		EC.params = &ec_params;

		//evaluating epsilon, new root and dE may change slightly
		epsilonFn(&F,&EC,&pp.dE,&epsilon,&minima);

		//evaluating mass about false vac
		mass2 = ddVd(minima[0],&paramsV);					////////// mass2

		//finding root0 of dV0(phi)=0;
		vector<double> minima0(3);
		if (pp.pot==1) {
			minima0[0] = -1.0; minima0[1] = 1.0;
		}
		else if (pp.pot==2) {
			gsl_function V0;
			V0.function = Vd_wrapped;
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
		double S1, S1error;
		gsl_function S1_integrand;
		S1_integrand.function = &s1Integrand;
		S1_integrand.params = &paramsV0;
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
		gsl_integration_qag(&S1_integrand, minima0[0], minima0[1], MIN_NUMBER, 1.0e-8, 1e4, 4, w, &S1, &S1error);
		gsl_integration_workspace_free(w);
		if (S1error>1.0e-8) cerr << "S1 error = " << S1error << endl;
		
		R = dE/S1;											////////// R
		action0 = -pi*epsilon*pow(R,2)/2.0 + pi*R*S1;		////////// action0
		L = pp.LoR*R;										////////// L
		if (pp.Tb<R) {
			double angle = asin(Tb/R);
			double Ltemp = 1.5*(1.5*Tb*tan(angle));
			if (Ltemp<L) L=Ltemp; //making sure to use the smaller of the two possible Ls
		}
	}
	else if (pot==3) {
		mass2 = 1.0;										////////// mass2
		minima[0] = 0.0, minima[1] = 0.0; 					////////// minima
		R = 10.0;											////////// R
		action0 = 8.0*pow(pi,2.0)/3.0;						////////// action0
		L = pp.LoR*R;										////////// L
	}
	
	a = L/(N-1.0);											////////// a
	b = Tb/(Nb-1.0);										////////// b
	Ta = b*(double)Na;										////////// Ta
	Tc = b*(double)Nc;										////////// Tc
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
Parameters::Parameters(): PrimaryParameters(), SecondaryParameters()	{}			

// constructor using primary and secondary parameters
Parameters::Parameters(const PrimaryParameters& p1, const SecondaryParmeters& p2): \
					PrimaryParmaeters(p1), SecondaryParameters(p2)	{}				

// print to shell
void Parameters::print() const {
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s8s%8s\n","N","Na","Nb","Nc","L","Ta","Tb","Tc","R","dE","theta","reg");
	printf("%8i%8i%8i%8i%8.4g%8.4g%8.4g%8.4g%8.4g%8.4g%8.4g%8.4g\n",\
			N,Na,Nb,Nc,\
			L,Ta,Tb,Tc,R,\
			dE,theta,reg);
	printf("\n");
}

// update parameters, uses setSecondaryParameters
void updateSecondaryParameters () {
	SecondaryParameters::setSecondaryParameters(PrimaryParameters);
}

// change all parameters due to change in one, uint
void changeParameters (const string& pName, const uint& pValue, struct Parameters&) {

}

// change all parameters due to change in one, double
void changeParameters (const string& pName, const double& pValue, struct Parameters&) {

}

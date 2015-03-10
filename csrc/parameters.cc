/*-------------------------------------------------------------------------------------------------------------------------
	definitions for functions and classes for dealing with parameters
-------------------------------------------------------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "simple.h"
#include "potentials.h"
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
		- operator<<
		- save
		- load
-------------------------------------------------------------------------------------------------------------------------*/


// operator<<
ostream& operator<<(ostream& os, const PrimaryParameters& p1) {
	os << left;
	os << setw(20) << "pot" << setw(20) << p1.pot << endl;
	os << setw(20) << "N" << setw(20) << p1.N << endl;
	os << setw(20) << "Na" << setw(20) << p1.Na << endl;
	os << setw(20) << "Nb" << setw(20) << p1.Nb << endl;
	os << setw(20) << "Nc" << setw(20) << p1.Nc << endl;
	os << setw(20) << "LoR" << setw(20) << p1.LoR << endl;
	os << setw(20) << "dE" << setw(20) << p1.dE << endl;
	os << setw(20) << "Tb" << setw(20) << p1.Tb << endl;
	os << setw(20) << "theta" << setw(20) << p1.theta << endl;
	os << setw(20) << "reg" << setw(20) << p1.reg << endl;
	os << endl;
	return os;
}

//save
void PrimaryParameters::save(const string& filename) const {
	ofstream os;
	os.open(filename.c_str());
	if (!os.good()) {
		FileError::StreamNotGood e(filename);
		throw e;
	}
	os << *this;
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
	string dross;
	is >> dross >> pot;
	is >> dross >> N;
	is >> dross >> Na;
	is >> dross >> Nb;
	is >> dross >> Nc;
	is >> dross >> LoR;
	is >> dross >> dE;
	is >> dross >> Tb;
	is >> dross >> theta;
	is >> dross >> reg;
	is.close();
}

/*-------------------------------------------------------------------------------------------------------------------------
	3.SecondaryParameters member functions
		- wrapped functions to take void* instead of a parameter struct, declared and defined here
		- set SecondaryParameters from PrimaryParameters
		- operator<<
-------------------------------------------------------------------------------------------------------------------------*/

//function pointer Vd, for gsl
static double (*Vd_local) (const double& phi, const params_for_V& parameters);
static double (*dVd_local) (const double& phi, const params_for_V& parameters);
static double (*ddVd_local) (const double& phi, const params_for_V& parameters);

// wrapped functions
#define WRAP2(FN) double FN##_wrapped2(double x, void* parameters) { return FN(x, *((params_for_V*)parameters)); }
WRAP2(Vd_local)
#undef WRAP2

// set secondary parameters
void SecondaryParameters::setSecondaryParameters (const struct PrimaryParameters& pp) {
	NT = pp.Na + pp.Nb + pp.Nc; 							////////// NT
	A = 0.4;												////////// A
	Gamma = exp(-pp.theta);									////////// Gamma
	
	double epsilon0;
	params_for_V paramsV, paramsV0;
	minima.reserve(2);
	//potential functions
	if (pp.pot==1) {
		Vd_local = &V1;
		dVd_local = &dV1;
		ddVd_local = &ddV1;
		epsilon0 = 0.0;
		epsilon = pp.dE;
		r0 = 0.0;											////////// r0
	}
	else if (pp.pot==2) {
		Vd_local = &V2;
		dVd_local = &dV2;
		ddVd_local = &ddV2;
		epsilon0 = 0.74507774287199924;
		epsilon = 0.75;
		r0 = 0.0;											////////// r0
	}
	else if (pp.pot==3) {
		Vd_local = &V3;
		dVd_local = &dV3;
		ddVd_local = &ddV3;
		epsilon0 = 0.0;
		epsilon = 0.0;
		r0 = MIN_NUMBER;									////////// r0
	}
	else
		cerr << "pot option not available, pot = " << pp.pot << endl;
	paramsV.epsi = epsilon;
	paramsV.aa = A;
	paramsV0.epsi = epsilon0;
	paramsV0.aa = A;
	
	// epsilon, minima, mass2, action0
	if (pp.pot!=3) {
		//gsl function for V(phi)
		gsl_function F;
		F.function = Vd_local_wrapped2;
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
		cout << "ec = " << endl;
		//cout << ec(epsilon,&ec_params) << endl;
		epsilonFn(&F,&EC,&pp.dE,&epsilon,&minima);
		cout << "test5" << endl;

		//evaluating mass about false vac
		mass2 = ddVd_local(minima[0],paramsV);					////////// mass2

		//finding root0 of dV0(phi)=0;
		vector<double> minima0(3);
		if (pp.pot==1) {
			minima0[0] = -1.0; minima0[1] = 1.0;
		}
		else if (pp.pot==2) {
			gsl_function V0;
			V0.function = Vd_local_wrapped2;
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
		
		R = pp.dE/S1;											////////// R
		action0 = -pi*epsilon*pow(R,2)/2.0 + pi*R*S1;		////////// action0
		L = pp.LoR*R;										////////// L
		if (pp.Tb<R) {
			double angle = asin(pp.Tb/R);
			double Ltemp = 1.5*(1.5*pp.Tb*tan(angle));
			if (Ltemp<L) L=Ltemp; //making sure to use the smaller of the two possible Ls
		}
	}
	else if (pp.pot==3) {
		mass2 = 1.0;										////////// mass2
		minima[0] = 0.0, minima[1] = 0.0; 					////////// minima
		R = 10.0;											////////// R
		action0 = 8.0*pow(pi,2.0)/3.0;						////////// action0
		L = pp.LoR*R;										////////// L
	}
	
	a = L/(pp.N-1.0);										////////// a
	b = pp.Tb/(pp.Nb-1.0);									////////// b
	Ta = b*(double)pp.Na;									////////// Ta
	Tc = b*(double)pp.Nc;									////////// Tc
}

// operator<<
ostream& operator<<(ostream& os, const SecondaryParameters& p2) {
	os << left;
	os << setw(20) << "NT" << setw(20) << p2.NT << endl;
	os << setw(20) << "epsilon" << setw(20) << p2.epsilon << endl;
	os << setw(20) << "R" << setw(20) << p2.R << endl;
	os << setw(20) << "Gamma" << setw(20) << p2.Gamma << endl;
	os << setw(20) << "r0" << setw(20) << p2.r0 << endl;
	os << setw(20) << "L" << setw(20) << p2.L << endl;
	os << setw(20) << "a" << setw(20) << p2.a << endl;
	os << setw(20) << "b" << setw(20) << p2.b << endl;
	os << setw(20) << "Ta" << setw(20) << p2.Ta << endl;
	os << setw(20) << "Tc" << setw(20) << p2.Tc << endl;
	os << setw(20) << "A" << setw(20) << p2.A << endl;
	os << setw(20) << "minima[0]" << setw(20) << p2.minima[0] << endl;
	os << setw(20) << "minima[1]" << setw(20) << p2.minima[1] << endl;
	os << setw(20) << "mass2" << setw(20) << p2.mass2 << endl;
	os << setw(20) << "action0" << setw(20) << p2.action0 << endl;
	os << endl;
	return os;
}

/*-------------------------------------------------------------------------------------------------------------------------
	4. Parameters member functions
		- empty constructor
		- constructor from PrimaryParameters and SecondaryParameters
		- print to shell
		- set secondary parameters
		- change parameters based on change in one
-------------------------------------------------------------------------------------------------------------------------*/

// empty constructor
Parameters::Parameters(): PrimaryParameters(), SecondaryParameters()	{}			

// constructor using primary and secondary parameters
Parameters::Parameters(const PrimaryParameters& p1, const SecondaryParameters& p2): \
					PrimaryParameters(p1), SecondaryParameters(p2)	{}				

// print to shell
void Parameters::print() const {
	printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n","N","Na","Nb","Nc","L","Ta","Tb","Tc","R","dE","theta","reg");
	printf("%8i%8i%8i%8i%8.4g%8.4g%8.4g%8.4g%8.4g%8.4g%8.4g%8.4g\n",\
			N,Na,Nb,Nc,\
			L,Ta,Tb,Tc,R,\
			dE,theta,reg);
	printf("\n");
}

// set parameters from PrimaryParameters
void Parameters::setSecondaryParameters () {
	PrimaryParameters p1;
	p1.pot = pot;
	p1.N = N;
	p1.Na = Na;
	p1.Nb = Nb;
	p1.Nc = Nc;
	p1.LoR = LoR;
	p1.dE = dE;
	p1.Tb = Tb;
	p1.theta = theta;
	p1.reg = reg;
	SecondaryParameters::setSecondaryParameters(p1);
}

// change all parameters due to change in one, uint
void Parameters::changeParameters (const string& pName, const uint& pValue, struct Parameters&) {

}

// change all parameters due to change in one, double
void Parameters::changeParameters (const string& pName, const double& pValue, struct Parameters&) {

}

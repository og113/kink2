/*-------------------------------------------------------------------------------------------------------------------------
	definitions for functions and classes for dealing with parameters
-------------------------------------------------------------------------------------------------------------------------*/

#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include "simple.h"
#include "potentials.h"
#include "gsl_extras.h"
#include "parameters.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - errors
	2 - PrimaryParameters
	3 - SecondaryParameters
	4 - Parameters
	5 - options
	6 - closenesses
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. parameter related errors
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
		- empty
		- operator==
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
	return os;
}

//save
void PrimaryParameters::save(const string& filename) const {
	try {
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
	catch (FileError::StreamNotGood & e) {
		cerr << "PrimaryParameters::save error:" << endl;
		cerr << e;
		return;
	}
}

//load
void PrimaryParameters::load(const string& filename) {
	try {
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
	catch (FileError::StreamNotGood& e) {
		cerr << "PrimaryParameters::load error:" << endl;
		cerr << e;
		return;
	}
}

// empty
bool PrimaryParameters::empty() const {
	return (pot==0 && N==0 && Na==0 && Nb==0 && Nc==0 && abs(LoR)<MIN_NUMBER && abs(dE)<MIN_NUMBER \
				&& abs(Tb)<MIN_NUMBER && abs(theta)<MIN_NUMBER && abs(reg)<MIN_NUMBER);
}

// operator==
bool operator==(const PrimaryParameters& l, const PrimaryParameters& r){
	return (l.pot==r.pot && l.N==r.N && l.Na==r.Na && l.Nb==r.Nb && l.Nc==r.Nc && abs(l.LoR-r.LoR)<MIN_NUMBER && \
				abs(l.dE-r.dE)<MIN_NUMBER && abs(l.Tb-r.Tb)<MIN_NUMBER && abs(l.theta-r.theta)<MIN_NUMBER && abs(l.reg-r.reg)<MIN_NUMBER);
}

// writeBinary
ostream& PrimaryParameters::writeBinary(ostream& os) const {
	os.write(reinterpret_cast<const char*>(&pot),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&N),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Na),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Nb),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Nc),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&LoR),sizeof(double));
	os.write(reinterpret_cast<const char*>(&dE),sizeof(double));
	os.write(reinterpret_cast<const char*>(&Tb),sizeof(double));
	os.write(reinterpret_cast<const char*>(&theta),sizeof(double));
	os.write(reinterpret_cast<const char*>(&reg),sizeof(double));
	return os;
}

// readBinary
istream& PrimaryParameters::readBinary(istream& is) {
	is.read(reinterpret_cast<char*>(&pot),sizeof(uint));
	is.read(reinterpret_cast<char*>(&N),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Na),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Nb),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Nc),sizeof(uint));
	is.read(reinterpret_cast<char*>(&LoR),sizeof(double));
	is.read(reinterpret_cast<char*>(&dE),sizeof(double));
	is.read(reinterpret_cast<char*>(&Tb),sizeof(double));
	is.read(reinterpret_cast<char*>(&theta),sizeof(double));
	is.read(reinterpret_cast<char*>(&reg),sizeof(double));
	return is;
}

/*-------------------------------------------------------------------------------------------------------------------------
	3.SecondaryParameters member functions
		- wrapped functions
		- VdV
		- dVddV
		- struct ec_params
		- ec (energy change: V(minima[1])-V(minima[0])-dE)
		- S1 integrand
		- rho integrand
		- set SecondaryParameters from PrimaryParameters
		- operator<<
		
	N.B. "static" before a function makes it local to the .cc file
-------------------------------------------------------------------------------------------------------------------------*/

//function pointer Vd, dVd, ddVd for gsl
static double (*Vd_local) (const double& phi, const params_for_V& parameters);
static double (*dVd_local) (const double& phi, const params_for_V& parameters);
static double (*ddVd_local) (const double& phi, const params_for_V& parameters);

// wrapped functions
#define WRAP(FN) double FN##_wrapped(double x, void* parameters) { return FN(x, *((params_for_V*)parameters)); }
WRAP(Vd_local)
WRAP(dVd_local)
WRAP(ddVd_local)
#undef WRAP

/*
//V FDF gsl function
static void VdV (double x, void * parameters, double * f, double* df) 
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	*f =  Vd_local_wrapped(x,params);
	*df = dVd_local_wrapped(x,params);
	}
	
//dV FDF gsl functions
static void dVddV (double x, void * parameters, double * f, double* df) 
	{
	struct params_for_V * params = (struct params_for_V *)parameters;
	*f =  dVd_local_wrapped(x,params);
	*df = ddVd_local_wrapped(x,params);
	}
*/
	
//energy change parameter struct
struct ec_params {double aa; double minima0; double minima1; double de; };

//energy change gsl function : V(minima[1])-V(minima[0])-dE
static double ec (double epsi, void * parameters) {
	struct ec_params * paramsIn = (struct ec_params *)parameters;
	struct params_for_V paramsOut;
	paramsOut.epsi = epsi;
	paramsOut.aa = (paramsIn->aa);
	double minima1 = (paramsIn->minima1);
	double minima0 = (paramsIn->minima0);
	double de = (paramsIn->de);
	return Vd_local(minima0,paramsOut) - Vd_local(minima1,paramsOut) - de;
}
	
//S1 integrand
static double s1Integrand (double x, void * parameters) {
	struct params_for_V * params = (struct params_for_V *)parameters;
	return pow(2.0*Vd_local(x,*params),0.5);
}

//rho integrand
double rhoIntegrand (double x, void * parameters) {
	struct params_for_V * params = (struct params_for_V *)parameters;
	return pow(2.0*Vd_local(x,*params),-0.5);
}

//program to find epsilon given gsl functions df and dE
static void epsilonFn (gsl_function * xF, gsl_function * xEC, const double * xdE, double * xEpsilon, vector<double>* xMinima)
	{
	double closenessdE = 1.0e-14;
	vector<double> dE_test(1);	dE_test[0] = 1.0;
	double newdE = *xdE;
	struct params_for_V * Fparameters = (struct params_for_V *) (*xF).params;
	struct ec_params * ECparameters = (struct ec_params *) (*xEC).params;
	unsigned int counter = 0;
	unsigned int maxCounter = 1e4;
	while (dE_test.back()>closenessdE)
		{
		//find roots of ec(epsilon)=0
		*xEpsilon = brentRootFinder(xEC,*xEpsilon,*xEpsilon/2.0,*xEpsilon*2.0);
		//assign new value of epsilon to xF
		(*Fparameters).epsi = *xEpsilon;
		(*xF).params = Fparameters;
		//finding new roots of dV(phi)=0
		(*xMinima)[0] = brentMinimum(xF,-1.0,-3.0,0.0);
		(*xMinima)[1] = brentMinimum(xF,1.2,0.5,3.0);
		//assign new roots to xECDF
		(*ECparameters).minima0 = (*xMinima)[0];
		(*ECparameters).minima1 = (*xMinima)[1];
		(*xEC).params = ECparameters;
		//evaluating new dE
		newdE = (*(*xEC).function)(*xEpsilon,ECparameters) + *xdE;
		//evaluating test
		if (abs(*xdE)>MIN_NUMBER) dE_test.push_back(abs((newdE-(*xdE))/(*xdE)));
		else 						dE_test.push_back(abs(newdE-(*xdE)));
		counter++;
		//test if too many runs
		if (counter>maxCounter)
			{
			cout << "epsilonFn error, more that " << maxCounter << " loops, consider reducing closenessdE" << endl;
			cout << "dE_test.back() = " << dE_test.back() << " , closenessdE = " << closenessdE << endl;
			cout << "dE = " << *xdE << " , minima[0] = " << (*xMinima)[0] << " , minima[1] = " << (*xMinima)[1];
			cout << " , epsilon = " << *xEpsilon << endl << endl;
			break;
			}
		}
	delete Fparameters;
	delete ECparameters;
	//*xdE = newdE;
	}

// set secondary parameters
void SecondaryParameters::setSecondaryParameters (const struct PrimaryParameters& pp) {
	NT = pp.Na + pp.Nb + pp.Nc; 							////////// NT
	A = 0.4;												////////// A
	Gamma = exp(-pp.theta);									////////// Gamma
	
	params_for_V paramsV, paramsV0;
	minima = vector<double>(2,0.0);
	minima0 = vector<double>(2,0.0);
	//potential functions
	if (pp.pot==1) {
		Vd_local = &V1;
		dVd_local = &dV1;
		ddVd_local = &ddV1;
		epsilon0 = 0.0;										////////// epsilon0
		epsilon = pp.dE;
		r0 = 0.0;											////////// r0
	}
	else if (pp.pot==2) {
		Vd_local = &V2;
		dVd_local = &dV2;
		ddVd_local = &ddV2;
		epsilon0 = 0.74507774287199924;						////////// epsilon0, guess
		epsilon = 0.75;
		r0 = 0.0;											////////// r0
	}
	else if (pp.pot==3) {
		Vd_local = &V3;
		dVd_local = &dV3;
		ddVd_local = &ddV3;
		epsilon0 = 0.0;										////////// epsilon0
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
	
		//finding root0 of dV0(phi)=0;
		if (pp.pot==1) {
			minima0[0] = -1.0; minima0[1] = 1.0;				////////// minima0
		}
		else if (pp.pot==2) {
			gsl_function V0;
			V0.function = Vd_local_wrapped;
			V0.params = &paramsV0;	
			minima0[0] = brentMinimum(&V0, -1.0, -3.0, 0.0);
			minima0[0] = brentMinimum(&V0, 1.2, 0.5, 3.0);
			struct ec_params ec0_params = { A, minima0[0], minima0[1], 0.0};
			gsl_function EC0;
			EC0.function = &ec;
			EC0.params = &ec0_params;
			double dE0 = 0.0;
			epsilonFn(&V0,&EC0,&dE0,&epsilon0,&minima0);		////////// epsilon0, minima0
		}
	
		if (abs(pp.dE)>MIN_NUMBER) {
			//gsl function for V(phi)
			gsl_function F;
			F.function = Vd_local_wrapped;
			F.params = &paramsV;	

			//finding preliminary minima of V(phi)
			minima[0] = brentMinimum(&F, -1.0, -3.0, 0.0);		////////// minima
			minima[1] = brentMinimum(&F, 1.2, 0.5, 3.0);		////////// minima

			//gsl function for V(root2)-V(root1)-dE
			struct ec_params ec_params = { A, minima[0], minima[1], pp.dE};
			gsl_function EC;
			EC.function = &ec;
			EC.params = &ec_params;

			//evaluating epsilon, new root
			epsilonFn(&F,&EC,&pp.dE,&epsilon,&minima);
		}
		else {
			minima = minima0;
			epsilon = epsilon0;
		}

		//evaluating mass about false vac
		mass2 = ddVd_local(minima[0],paramsV);					////////// mass2	

		//finding S1
		double S1, S1error;
		gsl_function S1_integrand;
		S1_integrand.function = &s1Integrand;
		S1_integrand.params = &paramsV0;
		gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
		gsl_integration_qag(&S1_integrand, minima0[0], minima0[1], MIN_NUMBER, 1.0e-8, 1e4, 4, w, &S1, &S1error);
		gsl_integration_workspace_free(w);
		if (S1error>1.0e-8) cerr << "S1 error = " << S1error << endl;
		
		if (abs(pp.dE)>MIN_NUMBER)
			R = S1/pp.dE;									////////// R
		else
			R = 10.0;										////////// R
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
		minima[0] = 0.0; minima[1] = 0.0; 					////////// minima
		minima0 = minima;									////////// minima0
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
	os << setw(20) << "epsilon0" << setw(20) << p2.epsilon0 << endl;
	os << setw(20) << "minima0[0]" << setw(20) << p2.minima0[0] << endl;
	os << setw(20) << "minima0[1]" << setw(20) << p2.minima0[1] << endl;
	return os;
}

/*-------------------------------------------------------------------------------------------------------------------------
	4. Parameters member functions
		- empty constructor
		- constructor from PrimaryParameters
		- constructor from PrimaryParameters and SecondaryParameters
		- copy
		- constructor from Parameters
		- operator=
		- print to shell
		- set secondary parameters
		- load
		- empty
		- change parameters based on change in one
		- <<
		- readBinary
		- writeBinary
-------------------------------------------------------------------------------------------------------------------------*/

// empty constructor
Parameters::Parameters(): PrimaryParameters(), SecondaryParameters()	{}			

// constructor from PrimaryParameters
Parameters::Parameters(const PrimaryParameters& p1): PrimaryParameters(p1) {
	setSecondaryParameters();
}

// constructor using primary and secondary parameters
Parameters::Parameters(const PrimaryParameters& p1, const SecondaryParameters& p2): \
					PrimaryParameters(p1), SecondaryParameters(p2)	{}

// print to shell
void Parameters::print() const {
	printf("%8s%8s%8s%8s%12s%12s%12s%12s%12s%12s%12s%12s\n","N","Na","Nb","Nc","L","Ta","Tb","Tc","R","dE","theta","reg");
	printf("%8i%8i%8i%8i%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g\n",\
			N,Na,Nb,Nc,\
			L,Ta,Tb,Tc,R,\
			dE,theta,reg);
	printf("\n");
}

// print File*
void Parameters::print(FILE* stream) const {
	fprintf(stream,"%8s%8s%8s%8s%12s%12s%12s%12s%12s%12s%12s%12s\n","N","Na","Nb","Nc","L","Ta","Tb","Tc","R","dE","theta","reg");
	fprintf(stream,"%8i%8i%8i%8i%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g\n",\
			N,Na,Nb,Nc,\
			L,Ta,Tb,Tc,R,\
			dE,theta,reg);
	fprintf(stream,"\n");
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

// load
void Parameters::load(const string& f) {
	PrimaryParameters::load(f);
	setSecondaryParameters();
}	

// empty
bool Parameters::empty() const {
	return PrimaryParameters::empty();
}

// change all parameters due to change in one, uint
bool Parameters::changeParameters (const string& pName, const uint& pValue) {
	bool anythingChanged = false;
	if ( pName.compare("N")==0) {
			if (N!=pValue) anythingChanged = true;
			N = pValue;
			a = L/(N-1);
		}
		else if ( pName.compare("Na")==0) {
			if (Na!=pValue) anythingChanged = true;
			Na = pValue;
			NT = Na + Nb + Nc;
			Ta = b*(double)Na;
			}
		else if ( pName.compare("Nb")==0) {
			if (Nb!=pValue) anythingChanged = true;
			Nb = pValue;
			NT = Na + Nb + Nc;
			b = Tb/(Nb-1.0);
			Ta = b*(double)Na;
			Tc = b*(double)Nc;
		}
		else if ( pName.compare("Nc")==0) {
			if (Nc!=pValue) anythingChanged = true;
			Nc = pValue;
			NT = Nc + Nb + Nc;
			Tc = b*(double)Nc;
		}
		else if ( pName.compare("pot")==0) { // would not recommend change pot this way
			if (pot!=pValue) anythingChanged = true;
			pot = pValue;
		}
		else {
			cerr << "Parameters::changeParameters error: " << pName << " not changed" << endl;
		}
	return anythingChanged;	
}

// change all parameters due to change in one, double
bool Parameters::changeParameters (const string& pName, const double& pValue) {
	bool anythingChanged = false;
	if ( pName.compare("L")==0) { //this does not changes the physics but simply the size of the box in space
		if (abs(L-pValue)>MIN_NUMBER) anythingChanged = true;
		L = pValue;
		LoR = L/R;
		a = L/(N-1);
	}
	if ( pName.compare("LoR")==0) { //this changes L, as above
		if (abs(LoR-pValue)>MIN_NUMBER) anythingChanged = true;
		LoR = pValue;
		L = LoR*R;
		a = L/(N-1);
	}
	else if ( pName.compare("Tb")==0) { //this paramter changes the physics for the periodic instanton,											//as Tb/R changes where R = R(epsilon)
		if (abs(Tb-pValue)>MIN_NUMBER) anythingChanged = true;
		b = b*pValue/Tb;
		Ta = Ta*pValue/Tb;
		Tc = Tc*pValue/Tb;
		Tb = pValue;
		if (Tb<R && pot!=3){
			double angle = asin(Tb/R);
			if (2.0*(1.5*Tb*tan(angle))<L) L=2.0*(1.5*Tb*tan(angle));
			a = L/(N-1.0);
			}
	}
	else if ( pName.compare("R")==0) { //this parameter changes the initial guess
		L = L*pValue/R; //all length scales scale with R
		a = a*pValue/R;
		b = b*pValue/R;
		Ta = Ta*pValue/R;
		Tb = Tb*pValue/R;
		Tc = Tc*pValue/R;
		R = pValue;
	}
	else if ( pName.compare("dE")==0) { //this parameter changes the physics of the potential
													//but it does not change Tb/R, where R(epsilon)
		if (abs(dE-pValue)>MIN_NUMBER) anythingChanged = true;
		R = R*dE/pValue; //R scales with 1/dE and the other length scales scale with R
		L = L*dE/pValue;
		a = a*dE/pValue;
		b = b*dE/pValue;
		Ta = Ta*dE/pValue;
		Tb = Tb*dE/pValue;
		Tc = Tc*dE/pValue;
		epsilon = epsilon*pValue/dE;
		dE = pValue;
	}
	else if ( pName.compare("theta")==0) {
		if (abs(theta-pValue)>MIN_NUMBER) anythingChanged = true;
		theta = pValue;
		Gamma = exp(-theta);
	}
	else {
			cerr << "Parameters::changeParameters error: " << pName << " not changed" << endl;
		}
	return anythingChanged;
}

// operator<< - just prints primary parameters
ostream& operator<<(ostream& os, const Parameters& p1) {
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
	return os;
}

// writeBinary
ostream& Parameters::writeBinary(ostream& os) const {
	return PrimaryParameters::writeBinary(os);
}

// readBinary
istream& Parameters::readBinary(istream& is) {
	PrimaryParameters::readBinary(is);
	setSecondaryParameters();
	return is;
}

/*-------------------------------------------------------------------------------------------------------------------------
	5. Options
		- changeOptions
		- operator<<
		- save
		- load
		- print
-------------------------------------------------------------------------------------------------------------------------*/

// changeOptions double
bool Options::changeOptions (const string& pName, const double& pValue) {
	bool anythingChanged = false;
	if (pName.compare("alpha")==0) {
		if (absDiff(alpha,pValue)>MIN_NUMBER) anythingChanged = true;
		alpha = pValue;
	}
	else if (pName.compare("open")==0) {
		if (absDiff(open,pValue)>MIN_NUMBER) anythingChanged = true;
		open = pValue;
	}
	else if (pName.compare("amp")==0) {
		if (absDiff(amp,pValue)>MIN_NUMBER) anythingChanged = true;
		amp = pValue;
	}
	else if (pName.compare("loopMin")==0) {
		if (absDiff(loopMin,pValue)>MIN_NUMBER) anythingChanged = true;
		loopMin = pValue;
	}
	else if (pName.compare("loopMax")==0) {
		if (absDiff(loopMax,pValue)>MIN_NUMBER) anythingChanged = true;
		loopMax = pValue;
	}
	else if (pName.compare("epsiTb")==0) {
		if (absDiff(epsiTb,pValue)>MIN_NUMBER) anythingChanged = true;
		epsiTb = pValue;
	}
	else if (pName.compare("epsiTheta")==0) {
		if (absDiff(epsiTheta,pValue)>MIN_NUMBER) anythingChanged = true;
		epsiTheta = pValue;
	}
	else
		cerr << "changeOptions error: " << pName << " not understood" << endl;
	return anythingChanged;
}

// changeOptions uint
bool Options::changeOptions (const string& pName, const uint& pValue) {
	bool anythingChanged = false;
	if (pName.compare("loops")==0) {
		if ((int)(loops-pValue)!=0) anythingChanged = true;
		loops = pValue;
	}
	else
		cerr << "changeOptions error: " << pName << " not understood" << endl;
	return anythingChanged;
}

// changeOptions string
bool Options::changeOptions (const string& pName, const string& pValue) {
	bool anythingChanged = false;
	if (pName.compare("zmx")==0) {
		if (zmx.compare(pValue)!=0) anythingChanged = true;
		zmx = pValue;
	}
	else if (pName.compare("zmt")==0) {
		if (zmx.compare(pValue)!=0) anythingChanged = true;
		zmt = pValue;
	}
	else if (pName.compare("bds")==0) {
		if (bds.compare(pValue)!=0) anythingChanged = true;
		bds = pValue;
	}
	else if (pName.compare("inF")==0) {
		if (inF.compare(pValue)!=0) anythingChanged = true;
		inF = pValue;
	}
	else if (pName.compare("minTimenumberLoad")==0) {
		if (minTimenumberLoad.compare(pValue)!=0) anythingChanged = true;
		minTimenumberLoad = pValue;
	}
	else if (pName.compare("maxTimenumberLoad")==0) {
		if (maxTimenumberLoad.compare(pValue)!=0) anythingChanged = true;
		maxTimenumberLoad = pValue;
	}
	else if (pName.compare("minfLoopLoad")==0) {
		if (minfLoopLoad.compare(pValue)!=0) anythingChanged = true;
		minfLoopLoad = pValue;
	}
	else if (pName.compare("maxfLoopLoad")==0) {
		if (maxfLoopLoad.compare(pValue)!=0) anythingChanged = true;
		maxfLoopLoad = pValue;
	}
	else if (pName.compare("minLoopLoad")==0) {
		if (minLoopLoad.compare(pValue)!=0) anythingChanged = true;
		minLoopLoad = pValue;
	}
	else if (pName.compare("maxLoopLoad")==0) {
		if (maxLoopLoad.compare(pValue)!=0) anythingChanged = true;
		maxLoopLoad = pValue;
	}
	else if (pName.compare("loopChoice")==0) {
		if (loopChoice.compare(pValue)!=0) anythingChanged = true;
		loopChoice = pValue;
	}
	else if (pName.compare("printChoice")==0) {
		if (printChoice.compare(pValue)!=0) anythingChanged = true;
		printChoice = pValue;
	}
	else
		cerr << "changeOptions error: " << pName << " not understood" << endl;
	return anythingChanged;
}

// operator<<
ostream& operator<<(ostream& os, const Options& o) {
	os << left;
	os << setw(20) << "alpha" << setw(20) << o.alpha << endl;
	os << setw(20) << "open" << setw(20) << o.open << endl;
	os << setw(20) << "amp" << setw(20) << o.amp << endl;
	os << setw(20) << "zmx" << setw(20) << o.zmx << endl;
	os << setw(20) << "zmt" << setw(20) << o.zmt << endl;
	os << setw(20) << "bds" << setw(20) << o.bds << endl;
	os << setw(20) << "inF" << setw(20) << o.inF << endl;
	os << setw(20) << "minTimenumberLoad" << setw(20) << o.minTimenumberLoad << endl;
	os << setw(20) << "maxTimenumberLoad" << setw(20) << o.maxTimenumberLoad << endl;
	os << setw(20) << "minfLoopLoad" << setw(20) << o.minfLoopLoad << endl;
	os << setw(20) << "maxfLoopLoad" << setw(20) << o.maxfLoopLoad << endl;
	os << setw(20) << "minLoopLoad" << setw(20) << o.minLoopLoad << endl;
	os << setw(20) << "maxLoopLoad" << setw(20) << o.maxLoopLoad << endl;
	os << setw(20) << "loopChoice" << setw(20) << o.loopChoice << endl;
	os << setw(20) << "loopMin" << setw(20) << o.loopMin << endl;
	os << setw(20) << "loopMax" << setw(20) << o.loopMax << endl;
	os << setw(20) << "epsiTb" << setw(20) << o.epsiTb << endl;
	os << setw(20) << "epsiTheta" << setw(20) << o.epsiTheta << endl;
	os << setw(20) << "loops" << setw(20) << o.loops << endl;
	os << setw(20) << "printChoice" << setw(20) << o.printChoice << endl;
	return os;
}

//save
void Options::save(const string& filename) const {
	ofstream os;
	os.open(filename.c_str());
	if (!os.good()) {
		FileError::StreamNotGood e(filename);
		cerr << e;
		return;
	}
	os << *this;
	os << endl;
	os.close();
}

//load
void Options::load(const string& filename) {
	ifstream is;
	is.open(filename.c_str());
	if (!is.good()) {
		FileError::StreamNotGood e(filename);
		cerr << e;
		return;
	}
	string dross;
	is >> dross >> alpha;
	is >> dross >> open;
	is >> dross >> amp;
	is >> dross >> zmx;
	is >> dross >> zmt;
	is >> dross >> bds;
	is >> dross >> inF;
	is >> dross >> minTimenumberLoad;
	is >> dross >> maxTimenumberLoad;
	is >> dross >> minfLoopLoad;
	is >> dross >> maxfLoopLoad;
	is >> dross >> minLoopLoad;
	is >> dross >> maxLoopLoad;
	is >> dross >> loopChoice;
	is >> dross >> loopMin;
	is >> dross >> loopMax;
	is >> dross >> epsiTb;
	is >> dross >> epsiTheta;
	is >> dross >> loops;
	is >> dross >> printChoice;
	is.close();
}

// print
void Options::print() const {
	printf("%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n","alpha","open","zmx","zmt","bds","inF","loopChoice","loopMin","loopMax","loops");
	printf("%12.1g%12.1g%12s%12s%12s%12s%12s%12.4g%12.4g%12i\n",\
			alpha,open,zmx.c_str(),zmt.c_str(),\
			bds.c_str(),inF.c_str(),loopChoice.c_str(),loopMin,loopMax,\
			loops);
	printf("\n");
}

// print File*
void Options::print(FILE* stream) const {
	fprintf(stream,"%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n","alpha","open","zmx","zmt","bds","inF","loopChoice","loopMin","loopMax","loops");
	fprintf(stream,"%12.1g%12.1g%12s%12s%12s%12s%12s%12.4g%12.4g%12i\n",\
			alpha,open,zmx.c_str(),zmt.c_str(),\
			bds.c_str(),inF.c_str(),loopChoice.c_str(),loopMin,loopMax,\
			loops);
	fprintf(stream,"\n");
}

/*-------------------------------------------------------------------------------------------------------------------------
	5. Closenesses
		- operator<<
		- save
		- load
-------------------------------------------------------------------------------------------------------------------------*/
	
// operator<<
ostream& operator<<(ostream& os, const Closenesses& c) {
	os << left;
	os << setw(20) << "Step" << setw(20) << c.Step << endl;
	os << setw(20) << "Action" << setw(20) << c.Action << endl;
	os << setw(20) << "Soln" << setw(20) << c.Soln << endl;
	os << setw(20) << "SolnMax" << setw(20) << c.SolnMax << endl;
	os << setw(20) << "Delta" << setw(20) << c.Delta << endl;
	os << setw(20) << "Inv" << setw(20) << c.Inv << endl;
	os << setw(20) << "Con" << setw(20) << c.Con << endl;
	os << setw(20) << "Lin" << setw(20) << c.Lin << endl;
	os << setw(20) << "True" << setw(20) << c.True << endl;
	os << setw(20) << "Latt" << setw(20) << c.Latt << endl;
	os << setw(20) << "Reg" << setw(20) << c.Reg << endl;
	os << setw(20) << "IE" << setw(20) << c.IE << endl;
	os << setw(20) << "Contm" << setw(20) << c.Contm << endl;
	os << setw(20) << "OS" << setw(20) << c.OS << endl;
	os << setw(20) << "AB" << setw(20) << c.AB << endl;
	os << setw(20) << "ABNE" << setw(20) << c.ABNE << endl;
	os << setw(20) << "LR" << setw(20) << c.LR << endl;
	os << setw(20) << "DT" << setw(20) << c.DT << endl;
	os << setw(20) << "Profile" << setw(20) << c.Profile << endl;
	return os;
}

// save
void Closenesses::save(const string& filename) const {
	ofstream os;
	os.open(filename.c_str());
	if (!os.good()) {
		FileError::StreamNotGood e(filename);
		cerr << e;
	}
	os << *this;
	os << endl;
	os.close();
}

//load
void Closenesses::load(const string& filename) {
	ifstream is;
	is.open(filename.c_str());
	if (!is.good()) {
		FileError::StreamNotGood e(filename);
		cerr << e;
	}
	string dross;
	is >> dross >> Step;
	is >> dross >> Action;
	is >> dross >> Soln;
	is >> dross >> SolnMax;
	is >> dross >> Delta;
	is >> dross >> Inv;
	is >> dross >> Con;
	is >> dross >> Lin;
	is >> dross >> True;
	is >> dross >> Latt;
	is >> dross >> Reg;
	is >> dross >> IE;
	is >> dross >> Contm;
	is >> dross >> OS;
	is >> dross >> AB;
	is >> dross >> ABNE;
	is >> dross >> LR;
	is >> dross >> DT;
	is >> dross >> Profile;
	is.close();
}

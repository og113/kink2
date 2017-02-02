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
	5 - ParametersRange
	6 - options
	7 - closenesses
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

// Size
const uint PrimaryParameters::Size = 10;
const uint ParametersRange::Size = 10;

// nameVector()
vector<string> PrimaryParameters::nameVector() const {
	vector<string> v(Size);
	v[0] = "Pot";
	v[1] =  "N";
	v[2] =  "Na";
	v[3] =  "Nb";
	v[4] =  "Nc";
	v[5] =  "LoR";
	v[6] =  "DE";
	v[7] =  "Tb";
	v[8] =  "Theta";
	v[9] =  "Reg";
	return v;
}

// getParameter
void PrimaryParameters::getParameter(const PrimaryParameters::Label& pLabel, uint& pValue) const {
	switch (pLabel) {
		case pot:
			pValue = Pot;
			break;
		case n:
			pValue = N;
			break;
		case na:
			pValue = Na;
			break;
		case nb:
			pValue = Nb;
			break;
		case nc:
			pValue = Nc;
			break;
		default:
			cerr << "Parameters::getParameter (uint) error: " << pLabel << " not found" << endl;
			break;
		}
}

// getParameter
void PrimaryParameters::getParameter(const PrimaryParameters::Label& pLabel, double& pValue) const {
	switch (pLabel) {
		case lor:
			pValue = LoR;
			break;
		case de:
			pValue = DE;
			break;
		case tb:
			pValue = Tb;
			break;
		case theta:
			pValue = Theta;
			break;
		case reg:
			pValue = Reg;
			break;
		default:
			cerr << "Parameters::getParameter (double) error: " << pLabel << " not found" << endl;
			break;
		}
}

// valueVector()
vector<string> PrimaryParameters::valueVector() const {
	vector<string> v(Size);
	v[0] = nts(Pot);
	v[1] =  nts(N);
	v[2] =  nts(Na);
	v[3] =  nts(Nb);
	v[4] =  nts(Nc);
	v[5] =  nts(LoR,16);
	v[6] =  nts(DE,16);
	v[7] =  nts(Tb,16);
	v[8] =  nts(Theta,16);
	v[9] =  nts(Reg,16);
	return v;
}

// operator<<
ostream& operator<<(ostream& os, const PrimaryParameters& p1) {
	os << left;
	vector<string> nv = p1.nameVector();
	vector<string> vv = p1.valueVector();
	for (uint j=0; j<p1.Size; j++) {
		os << setw(20) << nv[j] << setw(20) << vv[j] << endl;
	}
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

// load
void PrimaryParameters::load(const vector<string>& v) {
	if (v.size()!=Size) {
		cerr << "PrimaryParameters::load error: vector of size " << v.size() << "!=" << Size << endl;
	}
	Pot = stn<uint>(v[0]);
	N = stn<uint>(v[1]);
	Na = stn<uint>(v[2]);
	Nb = stn<uint>(v[3]);
	Nc = stn<uint>(v[4]);
	LoR = stn<number>(v[5]);
	DE = stn<number>(v[6]);
	Tb = stn<number>(v[7]);
	Theta = stn<number>(v[8]);
	Reg = stn<number>(v[9]);
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
	vector<string> v(Size);
	for (uint j=0; j<Size; j++)
		is >> dross >> v[j];
	is.close();
	load(v);
	}
	catch (FileError::StreamNotGood& e) {
		cerr << "PrimaryParameters::load error:" << endl;
		cerr << e;
		return;
	}
}

// empty
bool PrimaryParameters::empty() const {
	return (Pot==0 && N==0 && Na==0 && Nb==0 && Nc==0 && abs(LoR)<MIN_NUMBER && abs(DE)<MIN_NUMBER \
				&& abs(Tb)<MIN_NUMBER && abs(Theta)<MIN_NUMBER && abs(Reg)<MIN_NUMBER);
}

// operator==
bool operator==(const PrimaryParameters& l, const PrimaryParameters& r){
	return (l.Pot==r.Pot && l.N==r.N && l.Na==r.Na && l.Nb==r.Nb && l.Nc==r.Nc && abs(l.LoR-r.LoR)<MIN_NUMBER && \
				abs(l.DE-r.DE)<MIN_NUMBER && abs(l.Tb-r.Tb)<MIN_NUMBER && abs(l.Theta-r.Theta)<MIN_NUMBER && abs(l.Reg-r.Reg)<MIN_NUMBER);
}

// writeBinary
ostream& PrimaryParameters::writeBinary(ostream& os) const {
	os.write(reinterpret_cast<const char*>(&Pot),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&N),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Na),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Nb),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&Nc),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&LoR),sizeof(double));
	os.write(reinterpret_cast<const char*>(&DE),sizeof(double));
	os.write(reinterpret_cast<const char*>(&Tb),sizeof(double));
	os.write(reinterpret_cast<const char*>(&Theta),sizeof(double));
	os.write(reinterpret_cast<const char*>(&Reg),sizeof(double));
	return os;
}

// readBinary
istream& PrimaryParameters::readBinary(istream& is) {
	is.read(reinterpret_cast<char*>(&Pot),sizeof(uint));
	is.read(reinterpret_cast<char*>(&N),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Na),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Nb),sizeof(uint));
	is.read(reinterpret_cast<char*>(&Nc),sizeof(uint));
	is.read(reinterpret_cast<char*>(&LoR),sizeof(double));
	is.read(reinterpret_cast<char*>(&DE),sizeof(double));
	is.read(reinterpret_cast<char*>(&Tb),sizeof(double));
	is.read(reinterpret_cast<char*>(&Theta),sizeof(double));
	is.read(reinterpret_cast<char*>(&Reg),sizeof(double));
	return is;
}

// copying
void PrimaryParameters::copy(const PrimaryParameters& rhs) {
	Pot = rhs.Pot;
	N = rhs.N;
	Na = rhs.Na;
	Nb = rhs.Nb;
	Nc = rhs.Nc;
	Na = rhs.Na;
	LoR = rhs.LoR;
	DE = rhs.DE;
	Tb = rhs.Tb;
	Theta = rhs.Theta;
	Reg = rhs.Reg;
}

// copying
PrimaryParameters& PrimaryParameters::operator=(const PrimaryParameters& rhs) {
	copy(rhs);
	return *this;
}

/*-------------------------------------------------------------------------------------------------------------------------
	3.SecondaryParameters member functions
		- wrapped functions
		- VdV
		- dVddV
		- struct ec_params
		- ec (energy change: V(minima[1])-V(minima[0])-DE)
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

//energy change gsl function : V(minima[1])-V(minima[0])-DE
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

//program to find epsilon given gsl functions df and DE
static void epsilonFn (gsl_function * xF, gsl_function * xEC, const double * xDE, double * xEpsilon, vector<double>* xMinima)
	{
	double closenessDE = 1.0e-14;
	vector<double> DE_test(1);	DE_test[0] = 1.0;
	double newDE = *xDE;
	struct params_for_V * Fparameters = (struct params_for_V *) (*xF).params;
	struct ec_params * ECparameters = (struct ec_params *) (*xEC).params;
	unsigned int counter = 0;
	unsigned int maxCounter = 1e4;
	while (DE_test.back()>closenessDE)
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
		//evaluating new DE
		newDE = (*(*xEC).function)(*xEpsilon,ECparameters) + *xDE;
		//evaluating test
		if (abs(*xDE)>MIN_NUMBER) DE_test.push_back(abs((newDE-(*xDE))/(*xDE)));
		else 						DE_test.push_back(abs(newDE-(*xDE)));
		counter++;
		//test if too many runs
		if (counter>maxCounter)
			{
			cout << "epsilonFn error, more that " << maxCounter << " loops, consider reducing closenessDE" << endl;
			cout << "DE_test.back() = " << DE_test.back() << " , closenessDE = " << closenessDE << endl;
			cout << "DE = " << *xDE << " , minima[0] = " << (*xMinima)[0] << " , minima[1] = " << (*xMinima)[1];
			cout << " , epsilon = " << *xEpsilon << endl << endl;
			break;
			}
		}
	delete Fparameters;
	delete ECparameters;
	//*xDE = newDE;
	}

// set secondary parameters
void SecondaryParameters::setSecondaryParameters (const struct PrimaryParameters& pp) {
	NT = pp.Na + pp.Nb + pp.Nc; 							////////// NT
	A = 0.4;												////////// A
	Gamma = exp(-pp.Theta);									////////// Gamma
	
	params_for_V paramsV, paramsV0;
	minima = vector<double>(2,0.0);
	minima0 = vector<double>(2,0.0);
	//Potential functions
	if (pp.Pot==1) {
		Vd_local = &V1;
		dVd_local = &dV1;
		ddVd_local = &ddV1;
		epsilon0 = 0.0;										////////// epsilon0
		epsilon = pp.DE;
		r0 = 0.0;											////////// r0
	}
	else if (pp.Pot==2) {
		Vd_local = &V2;
		dVd_local = &dV2;
		ddVd_local = &ddV2;
		epsilon0 = 0.74507774287199924;						////////// epsilon0, guess
		epsilon = 0.75;
		r0 = 0.0;											////////// r0
	}
	else if (pp.Pot==3) {
		Vd_local = &V3;
		dVd_local = &dV3;
		ddVd_local = &ddV3;
		epsilon0 = 0.0;										////////// epsilon0
		epsilon = 0.0;
		r0 = MIN_NUMBER;									////////// r0
	}
	else {
		cerr << "Pot option not available, Pot = " << pp.Pot << endl;
		cerr << pp.Pot << endl;
	}
	paramsV.epsi = epsilon;
	paramsV.aa = A;
	paramsV0.epsi = epsilon0;
	paramsV0.aa = A;
	
	// epsilon, minima, mass2, action0
	if (pp.Pot!=3) {
	
		//finding root0 of dV0(phi)=0;
		if (pp.Pot==1) {
			minima0[0] = -1.0; minima0[1] = 1.0;				////////// minima0
		}
		else if (pp.Pot==2) {
			gsl_function V0;
			V0.function = Vd_local_wrapped;
			V0.params = &paramsV0;	
			minima0[0] = brentMinimum(&V0, -1.0, -3.0, 0.0);
			minima0[0] = brentMinimum(&V0, 1.2, 0.5, 3.0);
			struct ec_params ec0_params = { A, minima0[0], minima0[1], 0.0};
			gsl_function EC0;
			EC0.function = &ec;
			EC0.params = &ec0_params;
			double DE0 = 0.0;
			epsilonFn(&V0,&EC0,&DE0,&epsilon0,&minima0);		////////// epsilon0, minima0
		}
	
		if (abs(pp.DE)>MIN_NUMBER) {
			//gsl function for V(phi)
			gsl_function F;
			F.function = Vd_local_wrapped;
			F.params = &paramsV;	

			//finding preliminary minima of V(phi)
			minima[0] = brentMinimum(&F, -1.0, -3.0, 0.0);		////////// minima
			minima[1] = brentMinimum(&F, 1.2, 0.5, 3.0);		////////// minima

			//gsl function for V(root2)-V(root1)-DE
			struct ec_params ec_params = { A, minima[0], minima[1], pp.DE};
			gsl_function EC;
			EC.function = &ec;
			EC.params = &ec_params;

			//evaluating epsilon, new root
			epsilonFn(&F,&EC,&pp.DE,&epsilon,&minima);
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
		
		if (abs(pp.DE)>MIN_NUMBER)
			R = S1/pp.DE;									////////// R
		else
			R = 10.0;										////////// R
		action0 = -PI*epsilon*pow(R,2)/2.0 + PI*R*S1;		////////// action0
		L = pp.LoR*R;										////////// L
		if (pp.Tb<R) {
			double angle = asin(pp.Tb/R);
			double Ltemp = 1.5*(1.5*pp.Tb*tan(angle));
			if (Ltemp<L) L=Ltemp; //making sure to use the smaller of the two possible Ls
		}
	}
	else if (pp.Pot==3) {
		mass2 = 1.0;										////////// mass2
		minima[0] = 0.0; minima[1] = 0.0; 					////////// minima
		minima0 = minima;									////////// minima0
		R = 1.0;											////////// R
		action0 = 8.0*pow(PI,2.0)/3.0;						////////// action0
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
	printf("%8s%8s%8s%8s%8s%12s%12s%12s%12s%12s%12s%12s%12s\n","Pot","N","Na","Nb","Nc","L","Ta","Tb","Tc","R","DE","Theta","Reg");
	printf("%8i%8i%8i%8i%8i%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g\n",\
			Pot,N,Na,Nb,Nc,\
			L,Ta,Tb,Tc,R,\
			DE,Theta,Reg);
	printf("\n");
}

// print File*
void Parameters::print(FILE* stream) const {
	fprintf(stream,"%8s%8s%8s%8s%8s%12s%12s%12s%12s%12s%12s%12s%12s\n","Pot","N","Na","Nb","Nc","L","Ta","Tb","Tc","R","DE","Theta","Reg");
	fprintf(stream,"%8i%8i%8i%8i%8i%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g%12.5g\n",\
			Pot,N,Na,Nb,Nc,\
			L,Ta,Tb,Tc,R,\
			DE,Theta,Reg);
	fprintf(stream,"\n");
}

// getParameter, inherited
void Parameters::getParameter(const PrimaryParameters::Label& pLabel, uint& pValue) const {
	PrimaryParameters::getParameter(pLabel,pValue);
}

// getParameter, inherited
void Parameters::getParameter(const PrimaryParameters::Label& pLabel, double& pValue) const {
	PrimaryParameters::getParameter(pLabel,pValue);
}

// nameVector, inherited
vector<string> Parameters::nameVector() const {
	return PrimaryParameters::nameVector();
}

// valueVector, inherited
vector<string> Parameters::valueVector() const {
	return PrimaryParameters::valueVector();
}

// set parameters from PrimaryParameters
void Parameters::setSecondaryParameters () {
	PrimaryParameters p1;
	p1.Pot = Pot;
	p1.N = N;
	p1.Na = Na;
	p1.Nb = Nb;
	p1.Nc = Nc;
	p1.LoR = LoR;
	p1.DE = DE;
	p1.Tb = Tb;
	p1.Theta = Theta;
	p1.Reg = Reg;
	SecondaryParameters::setSecondaryParameters(p1);
}

// load
void Parameters::load(const vector<string>& vv) {
	PrimaryParameters::load(vv);
	setSecondaryParameters();
}	

// load
void Parameters::load(const string& f) {
	PrimaryParameters::load(f);
	setSecondaryParameters();
}	

// empty, inherited
bool Parameters::empty() const {
	return PrimaryParameters::empty();
}

// getLabel
PrimaryParameters::Label Parameters::getLabel(const string& pName) const {
	PrimaryParameters::Label label = static_cast<Parameters::Label>(0);
	if (pName.compare("Pot")==0 || pName.compare("pot")==0)
		label = PrimaryParameters::pot;
	else if (pName.compare("N")==0)
		label = PrimaryParameters::n;
	else if (pName.compare("Na")==0)
		label = PrimaryParameters::na;
	else if (pName.compare("Nb")==0)
		label = PrimaryParameters::nb;
	else if (pName.compare("Nc")==0)
		label = PrimaryParameters::nc;
	else if (pName.compare("LoR")==0)
		label = PrimaryParameters::lor;
	else if (pName.compare("DE")==0 || pName.compare("dE")==0)
		label = PrimaryParameters::de;
	else if (pName.compare("Tb")==0 )
		label = PrimaryParameters::tb;
	else if (pName.compare("Theta")==0 || pName.compare("theta")==0)
		label = PrimaryParameters::theta;
	else if (pName.compare("Reg")==0 || pName.compare("reg")==0)
		label = PrimaryParameters::reg;
	else {
		cerr << "PrimaryParameters error: label " << pName << " not understood" << endl;
	}
	return label;
}

// change all parameters due to change in one, uint
void Parameters::changeParameters (const PrimaryParameters::Label& label, const uint& pValue) {
	switch (label) {
		case pot:
			Pot = pValue;
			break;
		case n:
			N = pValue;
			a = L/(N-1);
			break;
		case na:
			Na = pValue;
			NT = Na + Nb + Nc;
			Ta = b*(double)Na;
			break;
		case nb:
			Nb = pValue;
			NT = Na + Nb + Nc;
			b = Tb/(Nb-1.0);
			Ta = b*(double)Na;
			Tc = b*(double)Nc;
			break;
		case nc:
			Nc = pValue;
			NT = Na + Nb + Nc;
			Tc = b*(double)Nc;
			break;
		default:
			cerr << "Parameters::changeParameters error: " << label << " not changed" << endl;
			break;
	}
}

// change all parameters due to change in one, uint
void Parameters::changeParameters (const string& pName, const uint& pValue) {
	PrimaryParameters::Label label = getLabel(pName);
	return changeParameters(label,pValue);
}

// change all parameters due to change in one, double
void Parameters::changeParameters (const PrimaryParameters::Label& label, const double& pValue) {
	switch (label) {
		case lor:
			LoR = pValue;
			L = LoR*R;
			a = L/(N-1);
			break;
		case de:
			//this parameter changes the physics of the Potential
			//but it does not change Tb/R, where R(epsilon)
			R = R*DE/pValue; //R scales with 1/DE and the other length scales scale with R
			L = L*DE/pValue;
			a = a*DE/pValue;
			b = b*DE/pValue;
			Ta = Ta*DE/pValue;
			Tb = Tb*DE/pValue;
			Tc = Tc*DE/pValue;
			epsilon = epsilon*pValue/DE;
			DE = pValue;
			break;
		case tb:
			//this paramter changes the physics for the periodic instanton,	
			//as Tb/R changes where R = R(epsilon)
			b = b*pValue/Tb;
			Ta = Ta*pValue/Tb;
			Tc = Tc*pValue/Tb;
			Tb = pValue;
			if (Tb<R && Pot!=3){
				double angle = asin(Tb/R);
				if (2.0*(1.5*Tb*tan(angle))<L) L=2.0*(1.5*Tb*tan(angle));
				a = L/(N-1.0);
			}
			break;
		case theta:
			Theta = pValue;
			Gamma = exp(-Theta);
			break;
		case reg:
			Reg = pValue;
			break;
		default:
			cerr << "Parameters::changeParameters error: " << label << " not changed" << endl;
			break;
	}
}


// change all parameters due to change in one, double
// n.b. i have removed the option to change secondary parameters such as L and R
void Parameters::changeParameters (const string& pName, const double& pValue) {
	PrimaryParameters::Label label = getLabel(pName);
	return changeParameters(label,pValue);
}

// change all parameters due to change in one, uint
void Parameters::step (const ParametersRange& pr, const PrimaryParameters::Label& label, const uint& j) {
	if ((uint)label<6) { // THIS IS HIGHLY IMPLEMENTATION DEPENDENT!
		uint current = 0, max = 0, min = 0;
		getParameter(label,current);
		(pr.Max).getParameter(label,max);
		(pr.Min).getParameter(label,min);
		uint pValue = current+j*(max-min)/((pr.Steps)[label-1]-1.0);
		changeParameters(label,pValue);
	}
	else {
		double current = 0.0, max = 0.0, min = 0.0;
		getParameter(label,current);
		(pr.Max).getParameter(label,max);
		(pr.Min).getParameter(label,min);
		double pValue = current+j*(max-min)/((pr.Steps)[label-1]-1.0);
		changeParameters(label,pValue);
	}
}

// step
void Parameters::step(const ParametersRange& pr, const Parameters::Label& label) {
	step(pr,label,1);
}

// operator<< - just prints primary parameters
ostream& operator<<(ostream& os, const Parameters& p1) {
	os << static_cast<const PrimaryParameters&>(p1);
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

// copying
void Parameters::copy(const Parameters& rhs) {
	PrimaryParameters::copy(static_cast<const PrimaryParameters&>(rhs));
	setSecondaryParameters();
}

// copying
Parameters& Parameters::operator=(const Parameters& rhs) {
	copy(rhs);
	return *this;
}

/*-------------------------------------------------------------------------------------------------------------------------
	5. ParametersRange
-------------------------------------------------------------------------------------------------------------------------*/

// ParametersRange empty constructor
ParametersRange::ParametersRange(): Min(), Max() {
	Steps = vector<uint>(Size,0);
}

// ParametersRange constructor
ParametersRange::ParametersRange(const Parameters& min, const Parameters& max, const vector<uint>& steps): \
									Min(min), Max(max) {
	if (steps.size()==Size)
		Steps = steps;
	else {
		Steps = vector<uint>(Size,0);
		cerr << "ParametersRange Error: Initialization not possible as input steps.size() = " << steps.size() << " not " << Size << endl;
	} 							
}

// toStep
bool ParametersRange::toStep(const Parameters::Label& stepNum) const {
	return (Steps[stepNum-1]>0);
}

// totalSteps
uint ParametersRange::totalSteps() const {
	uint N = 1;
	for (uint j=0; j<ParametersRange::Size; j++)
		N *= (Steps[j]>0? Steps[j]: 1);
	return N;
}

// position
Parameters ParametersRange::position(const uint& pos) const {
	Parameters p =  Min;
	uint n, Nratio, local;
	Nratio = totalSteps();
	local = pos;
	for (uint j=0; j<ParametersRange::Size; j++) {
		if (Steps[ParametersRange::Size-j-1]>0) {
			Nratio /= Steps[ParametersRange::Size-j-1];
			n = local/Nratio;
			local -= n*Nratio;
			p.step(*this,static_cast<Parameters::Label>(ParametersRange::Size-j),n);
		}
	}
	p.setSecondaryParameters();
	return p;
}

// neigh
Parameters ParametersRange::neigh(const uint& pos) const {
	if (pos==0) {
		cerr << "ParametersRange::neigh error: position 0 has no neighbour" << endl;
		return Min;
	}
	uint N = 1, k=0, npos = pos;
	while (k<ParametersRange::Size && npos==pos) {
		if (Steps[k]>0) {
			if (pos%(N*Steps[k])!=0)
				npos -= N;
			N *= Steps[k];
		}
		k++;
	}
	return position(npos);
}

// save
void ParametersRange::save(const string& filename) const {
	ofstream os;
	os.open(filename.c_str());
	if (!os.good()) {
		cerr << "ParametersRange::save error:" << endl;
		cerr << "stream for " << filename << " not good." << endl;
		return;
	}
	os << *this;
	os << endl;
	os.close();
}

// load
void ParametersRange::load(const string& filename) {
	ifstream is;
	is.open(filename.c_str());
	if (!is.good()) {
		cerr << "ParametersRange::load error:" << endl;
		cerr << "stream for " << filename << " not good" << endl;
		return;
	}
	string dross;
	// THIS NEEDS REDOING SO THAT THERE AREN'T MISTAKES WHEN SPACES DISAPPEAR
	is >> dross >> dross >> Min.Pot >> dross >> Max.Pot >> dross >> Steps[0] >> dross;
	is >> dross >> dross >> Min.N >> dross >> Max.N >> dross >> Steps[1] >> dross;
	is >> dross >> dross >> Min.Na >> dross >> Max.Na >> dross >> Steps[2] >> dross;
	is >> dross >> dross >> Min.Nb >> dross >> Max.Nb >> dross >> Steps[3] >> dross;
	is >> dross >> dross >> Min.Nc >> dross >> Max.Nc >> dross >> Steps[4] >> dross;
	is >> dross >> dross >> Min.LoR >> dross >> Max.LoR >> dross >> Steps[5] >> dross;
	is >> dross >> dross >> Min.DE >> dross >> Max.DE >> dross >> Steps[6] >> dross;
	is >> dross >> dross >> Min.Tb >> dross >> Max.Tb >> dross >> Steps[7] >> dross;
	is >> dross >> dross >> Min.Theta >> dross >> Max.Theta >> dross >> Steps[8] >> dross;
	is >> dross >> dross >> Min.Reg >> dross >> Max.Reg >> dross >> Steps[9] >> dross;
	is.close();
	Min.setSecondaryParameters();
	Max.setSecondaryParameters();
}

// empty
bool ParametersRange::empty() const {
	return (Min.empty() && Max.empty() && Steps.empty());
}

// writeBinary
ostream& ParametersRange::writeBinary(ostream& os) const {
	Min.writeBinary(os);
	Max.writeBinary(os);
	os.write(reinterpret_cast<const char*>(&Steps[0]),Size*sizeof(uint));
	return os;
}

// readBinary
istream& ParametersRange::readBinary(istream& is) {
	Min.readBinary(is);
	Max.readBinary(is);
	Min.setSecondaryParameters();
	Max.setSecondaryParameters();
	is.read(reinterpret_cast<char*>(&Steps[0]),Size*sizeof(uint));
	return is;
}

// operator<<
ostream& operator<<(ostream& os, const ParametersRange& pr) {
	os << left;
	vector<string> nv = (pr.Min).nameVector();
	vector<string> vvMin = (pr.Min).valueVector();
	vector<string> vvMax = (pr.Max).valueVector();
	for (uint j=0; j<pr.Size; j++)
		os << setw(20) << nv[j] << "[ " << setw(12) << vvMin[j] << " , " \
						<< setw(12) << vvMax[j] << " , " << setw(12) << (pr.Steps)[j] << " ]" << endl;
	return os;
}

// operator==
bool operator==(const ParametersRange& l, const ParametersRange& r) {
	return (l.Min==r.Min && l.Max==r.Max && l.Steps==r.Steps);
}

/*-------------------------------------------------------------------------------------------------------------------------
	6. Options
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
	7. Closenesses
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

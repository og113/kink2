/*-------------------------------------------------------------------------------------------------------------------------
	definitions for functions and classes for dealing with parameters
-------------------------------------------------------------------------------------------------------------------------*/

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
static double rhoIntegrand (double x, void * parameters) {
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
		if (abs(*xdE)>1.0e-16) dE_test.push_back(abs((newdE-(*xdE))/(*xdE)));
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
	//*xdE = newdE;
	}

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

		//evaluating epsilon, new root and dE may change slightly
		epsilonFn(&F,&EC,&pp.dE,&epsilon,&minima);

		//evaluating mass about false vac
		mass2 = ddVd_local(minima[0],paramsV);					////////// mass2

		//finding root0 of dV0(phi)=0;
		vector<double> minima0(3);
		if (pp.pot==1) {
			minima0[0] = -1.0; minima0[1] = 1.0;
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
		
		R = S1/pp.dE;										////////// R
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
		- load
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

// load
void Parameters::load(const string& f) {
	PrimaryParameters::load(f);
	setSecondaryParameters();
}	

// change all parameters due to change in one, uint
void Parameters::changeParameters (const string& pName, const uint& pValue) {
	if ( pName.compare("N")==0) {
			N = pValue;
			a = L/(N-1);
		}
		else if ( pName.compare("Na")==0) {
			Na = pValue;
			NT = Na + Nb + Nc;
			Ta = b*(double)Na;
			}
		else if ( pName.compare("Nb")==0) {
			Nb = pValue;
			NT = Na + Nb + Nc;
			b = Tb/(Nb-1.0);
			Ta = b*(double)Na;
			Tc = b*(double)Nc;
		}
		else if ( pName.compare("Nc")==0) {
			Nc = pValue;
			NT = Nc + Nb + Nc;
			Tc = b*(double)Nc;
		}
		else if ( pName.compare("pot")==0) {
			pot = pValue;
		}
		else {
			cerr << "Parameters::changeParameters error: " << pName << " not changed" << endl;
		}
	}

// change all parameters due to change in one, double
void Parameters::changeParameters (const string& pName, const double& pValue) {
	if ( pName.compare("L")==0) { //this does not changes the physics but simply the size of the box in space
		L = pValue;
		LoR = L/R;
		a = L/(N-1);
	}
	if ( pName.compare("LoR")==0) { //this changes L, as above
		LoR = pValue;
		L = LoR*R;
		a = L/(N-1);
	}
	else if ( pName.compare("Tb")==0) { //this paramter changes the physics for the periodic instanton,											//as Tb/R changes where R = R(epsilon)
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
		theta = pValue;
		Gamma = exp(-theta);
	}
	else {
			cerr << "Parameters::changeParameters error: " << pName << " not changed" << endl;
		}
}

/*-------------------------------------------------------------------------------------------------------------------------
	5. Options
		- operator<<
		- save
		- load
-------------------------------------------------------------------------------------------------------------------------*/

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
	os << setw(20) << "minLoopLoad" << setw(20) << o.minLoopLoad << endl;
	os << setw(20) << "maxLoopLoad" << setw(20) << o.maxLoopLoad << endl;
	os << setw(20) << "loopMin" << setw(20) << o.loopMin << endl;
	os << setw(20) << "loopMax" << setw(20) << o.loopMax << endl;
	os << setw(20) << "printChoice" << setw(20) << o.printChoice << endl;
	os << endl;
	return os;
}

//save
void Options::save(const string& filename) const {
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
void Options::load(const string& filename) {
	ifstream is;
	is.open(filename.c_str());
	if (!is.good()) {
		FileError::StreamNotGood e(filename);
		throw e;
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
	is >> dross >> minLoopLoad;
	is >> dross >> maxLoopLoad;
	is >> dross >> loopChoice;
	is >> dross >> loopMin;
	is >> dross >> loopMax;
	is >> dross >> printChoice;
	is.close();
}

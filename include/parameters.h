/*-------------------------------------------------------------------------------------------------------------------------
	declarations for functions and classes for dealing with parameters
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __PARAMETERS_H_INCLUDED__
#define __PARAMETERS_H_INCLUDED__

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include "error.h"

typedef unsigned int uint;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - errors
	2 - PrimaryParameters
	3 - SecondaryParameters
	4 - Parameters
	5 - ParametersRange
	6 - Options
	7 - Closenesses
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. parameter related errors
-------------------------------------------------------------------------------------------------------------------------*/

class ParameterError {
public:

	class Unset: public SimpleError{
	public:
		Unset(const string& s) : pName(s) {}					// constructor
		virtual string		message() const;					// message to be passed for printing
	private:
		string	pName;
	};
	
	class Load: public SimpleError{
	public:
		Load(const string& s) : filename(s) {}					// constructor
		virtual string		message() const;					// message to be passed for printing
	private:
		string	filename;
	};
};

/*-------------------------------------------------------------------------------------------------------------------------
	2. PrimaryParameters
		- PrimaryParameters
		- operator<<
		- operator==
	
	N.B. for structs the compiler automatically writes the following member functions, unless user provided:
	 empty constructor; copy constructor; assignment operator=; destructor
-------------------------------------------------------------------------------------------------------------------------*/

// PrimaryParameters
struct PrimaryParameters {
	enum Label { pot=1, n=2, na=3, nb=4, nc=5, lor=6, de=7, tb=8, theta=9, reg=10};
	static const uint Size;
	uint Pot;							// Potential
	uint N;
	uint Na;
	uint Nb;
	uint Nc;
	double LoR; 						// L/R
	double DE;
	double Tb;  						// BC section includes both corner points
	double Theta;
	double Reg;
	void getParameter(const PrimaryParameters::Label& pLabel, uint& pValue) const;
	void getParameter(const PrimaryParameters::Label& pLabel, double& pValue) const;
	vector<string> nameVector() const;
	vector<string> valueVector() const;
	void save(const string& filename) const;
	void load(const vector<string>& vv);
	void load(const string& filename);
	bool empty() const;
	ostream& writeBinary(ostream&) const;
	istream& readBinary(istream&);
	// copying
	void copy(const PrimaryParameters&);
	PrimaryParameters& operator=(const PrimaryParameters&);
};

// operator<<
ostream& operator<<(ostream&, const PrimaryParameters&);

// operator==
bool operator==(const PrimaryParameters& lhs, const PrimaryParameters& rhs);

/*-------------------------------------------------------------------------------------------------------------------------
	3. SecondaryParameters
		- SecondaryParameters
		- operator<<
		- VdV
		- dVddV
		- struct ec_params
		- ec (energy change: V(minima[1])-V(minima[0])-DE)
		- S1 integrand
		- rho integrand
		- epsilonFn
-------------------------------------------------------------------------------------------------------------------------*/

// SecondaryParameters
struct SecondaryParameters {
	uint NT;
	double epsilon;						// determined by DE
	double R; 							// size of bubble
	double Gamma; 						// equals exp(-Theta)
	double r0;							// minumum r (or x)
	double L;
	double a; 							// step sizes in each spatial dimension
	double b; 							// step sizes in time
	double Ta;
	double Tc;
	double A;							// mostly relevant for Pot=3
	vector<double> minima;				// minima of V
	double mass2; 						// as derived from V''
	double action0;						// action normalisation
	double epsilon0;					// value of epsilon at DE=0
	vector<double> minima0;				// positions of minima of V at DE=0
	void setSecondaryParameters (const struct PrimaryParameters&);				// sets secondary parameters using primary ones
};

// operator<<
ostream& operator<<(ostream&, const SecondaryParameters&);

//V FDF gsl function	
//void VdV (double x, void * parameters, double * f, double* df);

//dV FDF gsl function	
//void dVddV (double x, void * parameters, double * f, double* df);
	
//energy change parameter struct
//struct ec_params {double aa; double minima0; double minima1; double de; };

//energy change F gsl function : V(minima[1])-V(minima[0])-DE
//double ec (double epsi, void * parameters);
	
//S1 integrand
//double s1Integrand (double x, void * parameters);

//rho integrand
double rhoIntegrand (double x, void * parameters);

//program to find epsilon given gsls function df and DE
//void epsilonFn (gsl_function * xF, gsl_function * xEC, const double * xDE, double * xEpsilon, vector<double>* xMinima);

/*-------------------------------------------------------------------------------------------------------------------------
	4. Parameters
-------------------------------------------------------------------------------------------------------------------------*/

// ParametersRange declaration
struct ParametersRange;

struct Parameters: PrimaryParameters, SecondaryParameters {
	Parameters();																// empty constructor
	Parameters(const PrimaryParameters& p1);									// constructor using primary parameters
	Parameters(const PrimaryParameters& p1, const SecondaryParameters& p2);		// constructor using primary and secondary parameters
	~Parameters() {}
	void load(const string&);													// uses PrimaryParameters::load
	void load(const vector<string>& vv);										// uses PrimaryParameters::load
	void setSecondaryParameters();												// uses setSecondaryParameters
	PrimaryParameters::Label getLabel(const string& pName) const;
	void getParameter(const PrimaryParameters::Label& pLabel, uint& pValue) const;
	void getParameter(const PrimaryParameters::Label& pLabel, double& pValue) const;
	vector<string> nameVector() const;
	vector<string> valueVector() const;
	void changeParameters (const PrimaryParameters::Label& pLabel, const double& pValue); 			// change all due to change in one
	void changeParameters (const PrimaryParameters::Label& pLabel, const uint& pValue); 			// change all due to change in one
	void changeParameters (const string& pName, const double& pValue); 			// change all due to change in one
	void changeParameters (const string& pName, const uint& pValue); 	
	void step(const ParametersRange&, const PrimaryParameters::Label&);
	void step(const ParametersRange&, const PrimaryParameters::Label&, const uint&);
	void print() const;
	void print(FILE * stream) const;
	bool empty() const;
	ostream& writeBinary(ostream&) const;
	istream& readBinary(istream&);
	// copying
	void copy(const Parameters&);
	Parameters& operator=(const Parameters&);
};

// operator<<
ostream& operator<<(ostream&, const Parameters&);

/*-------------------------------------------------------------------------------------------------------------------------
	5. ParametersRange
-------------------------------------------------------------------------------------------------------------------------*/

// ParametersRange declaration
struct ParametersRange {
	static const uint Size;
	ParametersRange();
	ParametersRange(const Parameters& min, const Parameters& max, const vector<uint>& steps);
	Parameters 		Min;
	Parameters 		Max;
	vector<uint> 	Steps;
	uint			totalSteps() const;
	bool			toStep(const Parameters::Label&) const;
	Parameters		position(const uint&) const;
	Parameters		neigh(const uint&) const;
	void save(const string& filename) const;
	void load(const string& filename);
	bool empty() const;
	ostream& writeBinary(ostream&) const;
	istream& readBinary(istream&);
};

// operator<<
ostream& operator<<(ostream&, const ParametersRange&);

// operator==
bool operator==(const Parameters& lhs, const ParametersRange& rhs);

/*-------------------------------------------------------------------------------------------------------------------------
	6. Options
		- Options
		- operator<<
-------------------------------------------------------------------------------------------------------------------------*/

// Options
struct Options {
	double alpha;					// range over which thin wall spans
	double open;					// openness of input bubble between 0 and 1
	double amp;						// amplitude of negative mode oscillation
	string zmx;						// x translation zero mode fixing
	string zmt;						// t translation zero mode fixing
	string bds;						// initial boundary condition
	string inF;						// p, f, m, b - to load
	string minTimenumberLoad;
	string maxTimenumberLoad;
	string minfLoopLoad;
	string maxfLoopLoad;
	string minLoopLoad;
	string maxLoopLoad;
	string loopChoice;
	double loopMin;
	double loopMax;
	double epsiTb;
	double epsiTheta;
	uint loops;
	string printChoice;
	bool changeOptions (const string& pName, const double& pValue);
	bool changeOptions (const string& pName, const uint& pValue);
	bool changeOptions (const string& pName, const string& pValue);
	void save(const string& filename) const;
	void load(const string& filename);
	void print(FILE * stream) const;
	void print() const;
};

// operator<<
ostream& operator<<(ostream&, const Options&);

/*-------------------------------------------------------------------------------------------------------------------------
	7. Closenesses
		- Closenesses
		- operator<<
-------------------------------------------------------------------------------------------------------------------------*/

// Closenesses
struct Closenesses {
	double Step;
	double Action;
	double Soln;
	double SolnMax;
	double Delta;
	double Inv;
	double Con;
	double Lin;
	double True;
	double Latt;
	double Reg;
	double IE;
	double Contm;
	double OS;
	double AB;
	double ABNE;
	double LR;	
	double DT;
	double Profile;
	void save(const string& filename) const;
	void load(const string& filename);
};

// operator<<
ostream& operator<<(ostream&, const Closenesses&);

#endif // __PARAMETERS_H_INCLUDED__

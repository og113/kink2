/*-------------------------------------------------------------------------------------------------------------------------
	declarations for functions and classes for dealing with parameters
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __PARAMETERS_H_INCLUDED__
#define __PARAMETERS_H_INCLUDED__

#include <iostream>
#include <iomanip>
#include <map>
#include "map_extras.h"
#include "error.h"

typedef unsigned int uint;

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
	
	class Save: public SimpleError{
	public:
		Load(const string& s) : filename(s) {}					// constructor
		virtual string		message() const;					// message to be passed for printing
	private:
		string	filename;
	};
};

/*-------------------------------------------------------------------------------------------------------------------------
	2. PrimaryParameters
	
	N.B. for structs the compiler automatically writes the following member functions, unless user provided:
	 empty constructor; copy constructor; assignment operator=; destructor
-------------------------------------------------------------------------------------------------------------------------*/

struct PrimaryParameters {
	uint pot;							// potential
	uint N;
	uint Na;
	uint Nb;
	uint Nc;
	double LoR; 						// L/R
	double dE;
	double Tb;  						// BC section includes both corner points
	double theta;
	double reg;
	void save(const string& filename);
	void load(const string& filename);
};

/*-------------------------------------------------------------------------------------------------------------------------
	3. SecondaryParameters
-------------------------------------------------------------------------------------------------------------------------*/

struct SecondaryParameters {
	uint NT;
	double epsilon;						// determined by dE
	double R; 							// size of bubble
	double Gamma; 						// equals exp(-theta)
	double r0;							// minumum r (or x)
	double L;
	double a; 							// step sizes in each spatial dimension
	double b; 							// step sizes in time
	double Ta;
	double Tc;
	double A;							// mostly relevant for pot=3
	vector<double> minima(2);			// minima of V
	double mass2; 						// as derived from V''
	double action_0;					// action normalisation
	void setSecondaryParameters (const struct PrimaryParameters&);				// sets secondary parameters using primary ones
};

/*-------------------------------------------------------------------------------------------------------------------------
	4. Parameters
-------------------------------------------------------------------------------------------------------------------------*/

struct Parameters: PrimaryParameters, SecondaryParameters {
	Parameters();																// empty constructor
	Parameters(const PrimaryParameters& p1, const SecondaryParmeters& p2);		// constructor using primary and secondary parameters
	void updateSecondaryParameters ();											// uses setSecondaryParameters
	void changeParameters (const string& pName, const T& pValue, struct Parameters&); // change all due to change in one
	void print(const struct Parameters&);
};

#endif // __PARAMETERS_H_INCLUDED__

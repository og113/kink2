/*-------------------------------------------------------------------------------------------------------------------------
	declarations for functions giving newton-raphson quantities
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __NR_H_INCLUDED__
#define __NR_H_INCLUDED__

#include <complex>
#include "simple.h"
#include "parameters.h"
#include "potentials.h"
#include "folder.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	0 - typedefs
	1 - accessory functions
	2 - nr scalar functions
	3 - nr vector functions
	4 - nr matrix functions
	5 - filename functions
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	0 - typedefs
----------------------------------------------------------------------------------------------------------------------------*/

// most typedefs in simple.h

struct Complex_nr {
	enum Option { real, imaginary, both};
};

/*-------------------------------------------------------------------------------------------------------------------------
	1 - accessory functions
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	2 - nr scalar functions
-------------------------------------------------------------------------------------------------------------------------*/

// Kinetic_nr
void Kinetic_nr (const uint& j, const uint& dir, const vec& p, const Parameters& pr, const cVec& f, comp& result);

// Potential_nr
void Potential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& v, const cVec& f, comp& result);

// Boundary_nr
void Boundary_nr (const uint& j, const vec& p, const Parameters& pr, const Eigen::MatrixXd& omega_1, comp& result);

// Kinetic_nr
void Kinetic_nr (const uint& j, const uint& dir, const vec& p, const Parameters& pr, const cVec& f, cVec& result, const uint& pos);

// Potential_nr
void Potential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& v, const cVec& f, cVec& result, const uint& pos);

// Boundary_nr
void Boundary_nr (const uint& j, const vec& p, const Parameters& pr, const Eigen::MatrixXd& omega_1, cVec& result, const uint& pos);

/*-------------------------------------------------------------------------------------------------------------------------
	3 - nr vector functions
-------------------------------------------------------------------------------------------------------------------------*/

// mdKinetic_nr
void mdKinetic_nr (const uint& j, const uint& dir, const vec& p, const Parameters& pr, const cVec& f, vec& mds\
					, const Complex_nr::Option& opt=Complex_nr::both);

// mdPotential_nr
void mdPotential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& dv, const cVec& f, vec& mds\
					, const Complex_nr::Option& opt=Complex_nr::both);
					
// mdBoundary_nr
void mdBoundary_nr (const uint& j, const vec& p, const Parameters& pr, const Eigen::MatrixXd& omega_1\
					, vec& mds, const Complex_nr::Option& opt=Complex_nr::both);

/*-------------------------------------------------------------------------------------------------------------------------
	4 - nr matrix functions
-------------------------------------------------------------------------------------------------------------------------*/

// ddKinetic_nr
void ddKinetic_nr (const uint& j, const uint& dir, const vec& p, const Parameters& pr, const cVec& f, spMat& dds\
					, const Complex_nr::Option& opt=Complex_nr::both);

// ddPotential_nr
void ddPotential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& ddv, const cVec& f, spMat& dds\
					, const Complex_nr::Option& opt=Complex_nr::both);
					
// ddBoundary_nr
void ddBoundary_nr (const uint& j, const vec& p, const Parameters& pr, const Eigen::MatrixXd& omega_1\
					, spMat& dds, const Complex_nr::Option& opt=Complex_nr::both);

/*-------------------------------------------------------------------------------------------------------------------------
	5 - filename functions
-------------------------------------------------------------------------------------------------------------------------*/

// filenameMain
Filename filenameMain(const Parameters&, const string& base, const string& subfolder, const string& id, const string& suffix);

// filenameSpatial
Filename filenameSpatial(const Parameters&, const string& base, const string& subfolder, const string& id, const string& suffix);

// filenamePeriodic
Filename filenamePeriodic(const Parameters&, const string& base, const string& subfolder, const string& id, const string& suffix);

#endif // __NR_H_INCLUDED__

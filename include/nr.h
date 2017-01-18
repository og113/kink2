/*-------------------------------------------------------------------------------------------------------------------------
	declarations for functions giving newton-raphson quantities
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __NR_H_INCLUDED__
#define __NR_H_INCLUDED__

#include <complex>
#include "simple.h"
#include "parameters.h"
#include "potentials.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	0 - typedefs
	1 - accessory functions
	2 - nr scalar functions
	3 - nr vector functions
	4 - nr matrix functions
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	0 - typedefs
----------------------------------------------------------------------------------------------------------------------------*/

// most typedefs in simple.h

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

/*-------------------------------------------------------------------------------------------------------------------------
	3 - nr vector functions
-------------------------------------------------------------------------------------------------------------------------*/

// mdKinetic_nr
void mdKinetic_nr (const uint& j, const uint& dir, const vec& p, const Parameters& pr, const cVec& f, vec& mds);

// mdPotential_nr
void mdPotential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& v, const cVec& f, vec& mds);

/*-------------------------------------------------------------------------------------------------------------------------
	4 - nr matrix functions
-------------------------------------------------------------------------------------------------------------------------*/

// ddKinetic_nr
void ddKinetic_nr (const uint& j, const uint& dir, const vec& p, const Parameters& pr, const cVec& f, spMat& dds);

// ddPotential_nr
void ddPotential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& v, const cVec& f, spMat& dds);

#endif // __NR_H_INCLUDED__

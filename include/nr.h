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

// KineticS_nr
void KineticS_nr (const uint& j, const vec& p, const Parameters& pr, const cVec& f, comp& result);

// KineticT_nr
void KineticT_nr (const uint& j, const vec& p, const Parameters& pr, const cVec& f, comp& result);

// Potential_nr
void Potential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& v, const cVec& f, comp& result);

/*-------------------------------------------------------------------------------------------------------------------------
	3 - nr vector functions
-------------------------------------------------------------------------------------------------------------------------*/

// mdKineticS_nr
void mdKineticS_nr (const uint& j, const vec& p, const Parameters& pr, const cVec& f, vec& mds);


// mdKineticT_nr
void mdKineticT_nr (const uint& j, const vec& p, const Parameters& pr, const cVec& f, vec& mds);

/*-------------------------------------------------------------------------------------------------------------------------
	4 - nr matrix functions
-------------------------------------------------------------------------------------------------------------------------*/

// ddKineticS_nr
void ddKineticS_nr (const uint& j, const vec& p, const Parameters& pr, const cVec& f, spMat& dds);

// ddKineticT_nr
void ddKineticT_nr (const uint& j, const vec& p, const Parameters& pr, const cVec& f, spMat& dds);

#endif // __NR_H_INCLUDED__

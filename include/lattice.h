/*-------------------------------------------------------------------------------------------------------------------------
	declarations for functions and classes for dealing with the representation of the lattice
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __LATTICE_H_INCLUDED__
#define __LATTICE_H_INCLUDED__

#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "simple.h"
#include "error.h"
#include "parameters.h"

typedef unsigned int uint;
typedef long unsigned int lint;
typedef complex<double> comp;
typedef Eigen::VectorXd vec;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - integer lattice coord fns
	2 - real (or complex) lattice coord fns
	3 - neigh fns
	4 - dx, dt etc fns
	5 - interpolate fns
	6 - vecComplex, vecReal
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. integer lattice coord fns
		- c
		- intCoord
-------------------------------------------------------------------------------------------------------------------------*/

// c - inverse of intCoord
lint c(const uint& t, const uint& x, const Parameters& p);
lint c(const uint& t, const uint& x, const uint& Nt);

// intCoord
uint intCoord(const lint& loc, const uint& direc, const Parameters& p);
uint intCoord(const lint& loc, const uint& direc, const uint& Nt);

/*-------------------------------------------------------------------------------------------------------------------------
	2. real (or complex) lattice coord fns
		- simpleTime
		- simpleSpace
		- coord
		- coordB
-------------------------------------------------------------------------------------------------------------------------*/

// simpleTime
//static comp simpleTime(const uint& t, const Parameters& p);

// simpleSpace
//static double simpleSpace(const uint& c, const Parameters& p);

// coord
comp coord(const lint& loc, const int& direction, const Parameters& p);

// coordB
comp coordB(const lint& loc, const int& direction, const Parameters& p);

/*-------------------------------------------------------------------------------------------------------------------------
	3. neigh fns
		- periodic
		- spherical
		- neigh
-------------------------------------------------------------------------------------------------------------------------*/

//periodic
//static long periodic(const lint& loc, const uint& direction, const int& sign, const uint& xNt, const uint& xNx);

// spherical
//static long spherical(const lint& loc, const uint& direction, const int& sign, const uint& xNt, const uint& xNx);

// neigh - paramaters
long neigh (const lint& loc, const uint& direction, const int& sign, const Parameters p);

// neigh - Nt, Nx
long neigh (const lint& loc, const uint& direction, const int& sign, const uint& xNt, const uint& xNx);

/*-------------------------------------------------------------------------------------------------------------------------
	4. dx, dt etc
		- dtFn
		- DtFn
		- dxFn
		- DxFn
-------------------------------------------------------------------------------------------------------------------------*/

// dtFn
comp dtFn (const unsigned int& time, const Parameters& p);

// DtFn	
comp DtFn (const unsigned int& time, const Parameters& p);

// dxFn	
double dxFn (const unsigned int& space, const Parameters& p);

// DxFn	
double DxFn (const unsigned int& space, const Parameters& p);

/*-------------------------------------------------------------------------------------------------------------------------
	5. interpolate fns
		- interpolate (for real representation of complex 2d vector)
		- interpolateReal (for fundamentally real 2d vector)
		- interpolate1d (for 1d real vector)
-------------------------------------------------------------------------------------------------------------------------*/

//interpolate 2d (real representation of ) complex vector
vec interpolate(vec vec_old, const Parameters& p_old, const Parameters& p_new);
cVec interpolate(cVec vec_old, const Parameters& p_old, const Parameters& p_new);

//interpolate 2d real vector function
vec interpolateReal(vec vec_old, const Parameters& p_old, const Parameters& p_new);
cVec interpolateReal(cVec vec_old, const Parameters& p_old, const Parameters& p_new);

// interpolate1d	
vec interpolate1d(vec vec_old, const unsigned int & N_old, const unsigned int & N_new);
cVec interpolate1d(cVec vec_old, const unsigned int & N_old, const unsigned int & N_new);

/*-------------------------------------------------------------------------------------------------------------------------
	6. vecComplex etc
		- vecComplex
		- vecReal
-------------------------------------------------------------------------------------------------------------------------*/

//complexify a real vector - tDim
cVec vecComplex(vec realVec, const uint& tDim);

//complexify a real vector - parameters
cVec vecComplex(vec realVec, const Parameters& p);
	
//make a complex vector real - tDim
vec vecReal(cVec complexVec, const uint& tDim);

//make a complex vector real - parameters
vec vecReal(cVec complexVec, const Parameters& p);

#endif // __LATTICE_H_INCLUDED__

/*-------------------------------------------------------------------------------------------------------------------------
 	declarations for the functions to calculate omega
 -------------------------------------------------------------------------------------------------------------------------*/
 
#ifndef __OMEGA_H_INCLUDED__
#define __OMEGA_H_INCLUDED__

#include <Eigen/Dense>
#include "parameters.h"

using namespace std;

typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. h matrix
	2. analytic modes
	3. numerical modes
	4. constructing omegas
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. h matrix
		- hFn
-------------------------------------------------------------------------------------------------------------------------*/

// hFn
//static mat	hFn(const Parameters&);

/*-------------------------------------------------------------------------------------------------------------------------
	2. analytic modes
-------------------------------------------------------------------------------------------------------------------------*/

// analytic modes
void analyticModes(mat& modes, vec& freqs, vec& freqs_exp, const Parameters& p);

/*-------------------------------------------------------------------------------------------------------------------------
	3. numerical modes
-------------------------------------------------------------------------------------------------------------------------*/

// numerical modes
void numericalModes(mat& modes, vec& freqs, vec& freqs_exp, const Parameters& p);

/*-------------------------------------------------------------------------------------------------------------------------
	4. constructing omegas
-------------------------------------------------------------------------------------------------------------------------*/

// omegasFn
void omegasFn(const mat& modes, const mat& freqs, mat& omega_m1, mat& omega_0, mat& omega_1, mat& omega_2, const Parameters& p);

#endif // __OMEGA_H_INCLUDED__

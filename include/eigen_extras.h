/*-------------------------------------------------------------------------------------------------------------------------
	declarations for functions using eigen library functions
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __EIGEN_EXTRAS_H_INCLUDED__
#define __EIGEN_EXTRAS_H_INCLUDED__

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include "simple.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - printErrorInformation
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
 1. printErrorInformation
-------------------------------------------------------------------------------------------------------------------------*/

// printErrorInformation
void printErrorInformation(const vec& v, const string& name);

// printErrorInformation
void printErrorInformation(const vec& v, const string& name, const uint zm);

// printErrorInformation
void printErrorInformation(const mat& m, const string& name);

// printErrorInformation
void printErrorInformation(const spMat& m, const string& name);
	
#endif // __EIGEN_EXTRAS_H_INCLUDED__

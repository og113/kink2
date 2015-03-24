/*-------------------------------------------------------------------------------------------------------------------------
	declarations for global variables and such for main.cc
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __MAIN_H_INCLUDED__
#define __MAIN_H_INCLUDED__

#include <complex>
#include "check.h"
#include "error.h"
#include "fnptrs.h"
#include "folder.h"
#include "lattice.h"
#include "omega.h"
#include "parameters.h"
#include "potentials.h"
#include "print.h"
#include "simple.h"
#include "stepper.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - typedefs
	2 - definitions
	3 - global variables
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. typedefs
		- uint, lint, comp
		- vec, cVec, mat, cMat, spMat
-------------------------------------------------------------------------------------------------------------------------*/

typedef unsigned int uint;
typedef unsigned long lint;
typedef complex<double> comp;

typedef pair<string,string> StringPair;

typedef comp(*PotentialType)(const comp&, const struct params_for_V&);

typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;
typedef Eigen::SparseMatrix<double> spMat;

/*-------------------------------------------------------------------------------------------------------------------------
	2. definitions
		- pi
		- MIN_NUMBER
-------------------------------------------------------------------------------------------------------------------------*/

// pi
#ifndef pi
#define pi 3.14159265359
#endif

// MIN_NUMBER
#ifndef MIN_NUMBER
#define MIN_NUMBER 1.0e-16
#endif

/*-------------------------------------------------------------------------------------------------------------------------
	3. global variables
		- ii
-------------------------------------------------------------------------------------------------------------------------*/

comp ii(0.0,1.0);

#endif // __MAIN_H_INCLUDED__

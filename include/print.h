/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for functions to save, load and plot
 -------------------------------------------------------------------------------------------------------------------------*/
 
#ifndef __PRINT_H_INCLUDED__
#define __PRINT_H_INCLUDED__

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "parameters.h"
#include "lattice.h"

using namespace std;

typedef unsigned int uint;
typedef unsigned long lint;

typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. saveOptions
	2. save
	3. load
	4. plot
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. saveOptions
		- saveOptions
	
	n.b. saveOptions will also serve as loadOptions
-------------------------------------------------------------------------------------------------------------------------*/

// saveOptions
struct saveOptions {
	enum vectorType { simple=0, real=1, complex=2, realB=3, complexB=4 };
	enum extras { none=0, loc=1, coords=2};
	uint column;
	Parameters paramsIn;
	Parameters paramsOut;
};

/*-------------------------------------------------------------------------------------------------------------------------
	2. save
		- vec
		- cVec
		- mat
		- cMat
		- spMat
-------------------------------------------------------------------------------------------------------------------------*/

// save vec
void save(const string&, const saveOptions&, const vec&);

// save cVec
void save(const string&, const saveOptions&, const cVec&);

// save mat
void save(const string&, const saveOptions&, const mat&);

// save cMat
void save(const string&, const saveOptions&, const cMat&);

// save spMat
void save(const string&, const saveOptions&, const spMat&);

/*-------------------------------------------------------------------------------------------------------------------------
	3. load
		- vec
		- cVec
		- mat
		- cMat
		- spMat
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	4. plot
-------------------------------------------------------------------------------------------------------------------------*/

#endif // __PRINT_H_INCLUDED__

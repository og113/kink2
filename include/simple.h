/*-------------------------------------------------------------------------------------------------------------------------
declarations of some very simple functions and classes
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __SIMPLE_H_INCLUDED__
#define __SIMPLE_H_INCLUDED__

#include <string>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "simple.h"

using namespace std;

#define MIN_NUMBER 1.0e-16

typedef complex<double> comp;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - number to string and string to number
	2 - absDiff
	3 - factorial
	4 - currentDateTime
	5 - copyFile
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. number to string and string to number
		- number to string
		- string to number
-------------------------------------------------------------------------------------------------------------------------*/

template <class T>
string numberToString ( const T& Number );

template <class T>
T stringToNumber ( const string& Text );

/*-------------------------------------------------------------------------------------------------------------------------
	2. absDiff
		- double
		- comp
		- vec
		- cVec
		- mat
		- cMat
-------------------------------------------------------------------------------------------------------------------------*/

double absDiff (const double& A, const double& B);
double absDiff (const comp& A, const comp& B);
double absDiff (const vec& A, const vec& B);
double absDiff (const cVec& A, const cVec& B);
double absDiff (const mat& A, const mat& B);
double absDiff (const cMat& A, const cMat& B);

/*-------------------------------------------------------------------------------------------------------------------------
	3. factorial
-------------------------------------------------------------------------------------------------------------------------*/

//factorial
int factorial(const int& f_input);

/*-------------------------------------------------------------------------------------------------------------------------
	4. currentDateTime
-------------------------------------------------------------------------------------------------------------------------*/

//getting the date and time
string currentDateTime();

/*-------------------------------------------------------------------------------------------------------------------------
	5. copyFile
-------------------------------------------------------------------------------------------------------------------------*/

//copy a file
void copyFile(const string & inputFile, const string & outputFile);

#endif // __SIMPLE_H_INCLUDED__

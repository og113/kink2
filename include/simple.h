/*-------------------------------------------------------------------------------------------------------------------------
declarations of some very simple functions and classes
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __SIMPLE_H_INCLUDED__
#define __SIMPLE_H_INCLUDED__

#include <string>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "simple.h"

using namespace std;

#ifndef pi
#define pi 3.14159265359
#endif
#ifndef MIN_NUMBER
#define MIN_NUMBER 1.0e-16
#endif

typedef double number;
typedef unsigned int uint;
typedef unsigned long int lint;
typedef complex<double> comp;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;
typedef Eigen::SparseMatrix<double> spMat;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - number to string and string to number
	2 - absDiff
	3 - factorial
	4 - currentDateTime, currentPartSec
	5 - copyFile
	6 - count in files
	7 - smallestLoc
	8 - randDouble
	9 - mod
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
string nts ( const T& Number );

template <class T>
string nts ( const T& Number, const uint& prec);

template <class T>
T stringToNumber ( const string& Text );

template <class T>
T stn ( const string& Text );

bool isNumber( const string& Text );

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
	4. currentDateTime, currentPartSec
-------------------------------------------------------------------------------------------------------------------------*/

//getting the date and time
string currentDateTime();

// getting the part of the second
string currentPartSec();

/*-------------------------------------------------------------------------------------------------------------------------
	5. copyFile
-------------------------------------------------------------------------------------------------------------------------*/

//copy a file
void copyFile(const string & inputFile, const string & outputFile);

/*-------------------------------------------------------------------------------------------------------------------------
	6. count in file
		- countLines
		- countColumns
		- countDoubles
		- countComp
		- countType
-------------------------------------------------------------------------------------------------------------------------*/

// count lines
uint countLines(const string & f);

// count columns
uint countColumns(const string& f);

// count doubles (binary)
uint countDoubles(const string& f);

// count comp (binary)
uint countComp(const string& f);

// count type (binary)
template <class T>
uint countType(const string& f, const T& t);

/*-------------------------------------------------------------------------------------------------------------------------
	7. smallestLoc
-------------------------------------------------------------------------------------------------------------------------*/

// smallestLoc
template <class T>
uint smallestLoc(const vector<T>&);

/*-------------------------------------------------------------------------------------------------------------------------
	8. randDouble
	
	n.b. must first seed rand with srand(time(NULL)) or something similar
-------------------------------------------------------------------------------------------------------------------------*/

// randDouble
double randDouble(const double& min, const double& max);

/*-------------------------------------------------------------------------------------------------------------------------
	9. mod
-------------------------------------------------------------------------------------------------------------------------*/

// mod
double mod(const double& x, const double& min, const double& max);

#endif // __SIMPLE_H_INCLUDED__

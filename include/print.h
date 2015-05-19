/*-------------------------------------------------------------------------------------------------------------------------
 	declarations for functions to save, load and plot
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
typedef complex<double> comp;

typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;
typedef Eigen::SparseMatrix<double> spMat;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. SaveOptions
	2. save
	3. load
	4. plot
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. SaveOptions
		- SaveOptions
		- <<
	
	n.b. SaveOptions will also serve as loadOptions
-------------------------------------------------------------------------------------------------------------------------*/

// SaveOptions
struct SaveOptions {
	enum printTypeList { ascii=0, binary=1};
	enum vectorTypeList { simple=0, real=1, complex=2, realB=3, complexB=4, append=5 };
	enum extrasList { none=0, loc=1, coords=2, coordT=3, coordX=4};
	printTypeList printType;
	vectorTypeList vectorType;
	extrasList extras;
	uint column;					// for load this overrides the other options
	uint zeroModes;
	Parameters paramsIn;
	Parameters paramsOut;
	bool printMessage;
	bool good() const;					// checking if all required fields present
};

/*os.write(reinterpret_cast<const char*>(&opts.printType),sizeof(SaveOptions::printTypeList));
os.write(reinterpret_cast<const char*>(&opts.vectorType),sizeof(SaveOptions::vectorTypeList));
os.write(reinterpret_cast<const char*>(&opts.extras),sizeof(SaveOptions::extrasList));
os.write(reinterpret_cast<const char*>(&opts.column),sizeof(uint));
os.write(reinterpret_cast<const char*>(&opts.zeroModes),sizeof(uint));
os.write(reinterpret_cast<const char*>(&opts.paramsIn),sizeof(Parameters));
os.write(reinterpret_cast<const char*>(&opts.paramsOut),sizeof(Parameters));
os.write(reinterpret_cast<const char*>(&opts.printMessage),sizeof(bool));*/


// operator <<
ostream& operator<<(ostream& os, const SaveOptions& opts); 

/*-------------------------------------------------------------------------------------------------------------------------
	2. save
		- vec
		- cVec
		- mat
		- cMat
		- spMat
-------------------------------------------------------------------------------------------------------------------------*/

// save vec
void save(const string&, const SaveOptions&, const vec&);

// save cVec
void save(const string&, const SaveOptions&, const cVec&);

// save mat
void save(const string&, const SaveOptions&, const mat&);

// save cMat
void save(const string&, const SaveOptions&, const cMat&);

// save spMat
void save(const string&, const SaveOptions&, const spMat&);

/*-------------------------------------------------------------------------------------------------------------------------
	3. load
		- vec
		- cVec
		- mat
		- cMat
		- spMat
-------------------------------------------------------------------------------------------------------------------------*/

// load vec
void load(const string&, SaveOptions&, vec&);

// load cVec
void load(const string&, SaveOptions&, cVec&);

// load mat
void load(const string&, SaveOptions&, mat&);

// load cMat
void load(const string&, SaveOptions&, cMat&);

// load spMat
void load(const string&, SaveOptions&, spMat&);

/*-------------------------------------------------------------------------------------------------------------------------
	4. plot
		- PlotOptions
		- plot
-------------------------------------------------------------------------------------------------------------------------*/

// PlotOptions
struct PlotOptions {
	SaveOptions::printTypeList printType;
	string gp;
	uint column;
	uint column2;
	uint column3;
	uint column4;
	string style;
	string output;
	bool printMessage;
	vector<string> commands;
};

// plot
void plot (const string& filename, const PlotOptions& opts);

#endif // __PRINT_H_INCLUDED__

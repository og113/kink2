/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for functions to save, load and plot
 -------------------------------------------------------------------------------------------------------------------------*/
 
 #include <fstream>
 #include <iostream>
 #include <sstream>
 #include <iomanip>
 #include <cmath>
 #include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "folder.h"
#include "simple.h" // for countLines
#include "parameters.h"
#include "lattice.h"
#include "print.h"

using namespace std;

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
		- good
		- <<
-------------------------------------------------------------------------------------------------------------------------*/

// good
bool SaveOptions::good() const {
	//if (!paramsIn.empty() && paramsOut.empty())
	//	paramsOut = paramsIn;
	switch(vectorType) {
		case SaveOptions::simple:	if (extras==SaveOptions::none || extras==SaveOptions::loc)
										return true;
									else 
										return (!paramsIn.empty() && !paramsOut.empty());
									break;
		case SaveOptions::real:		return (!paramsIn.empty() && !paramsOut.empty());
									break;
		case SaveOptions::realB:	return (!paramsIn.empty() && !paramsOut.empty());
									break;
		case SaveOptions::complex:	return (!paramsIn.empty() && !paramsOut.empty());
									break;
		case SaveOptions::complexB:	return (!paramsIn.empty() && !paramsOut.empty());
									break;
		case SaveOptions::append:	return (printType!=SaveOptions::binary);
									break;
		default:					cerr << "SaveOptions error: vectorType " << vectorType << "not recognized" << endl;
									return false;
									break;
	}
}

// operator <<
ostream& operator<<(ostream& os, const SaveOptions& opts){
	os << opts.printType << endl;
	os << opts.vectorType << endl;
	os << opts.extras << endl;
	os << opts.column << endl;
	os << opts.zeroModes << endl;
	os << opts.paramsIn << endl;
	os << opts.paramsOut << endl;
	os << opts.printMessage << endl;
	return os;
}

// writeBinary
ostream& SaveOptions::writeBinary(ostream& os) const {
	os.write(reinterpret_cast<const char*>(&printType),sizeof(SaveOptions::printTypeList));
	os.write(reinterpret_cast<const char*>(&vectorType),sizeof(SaveOptions::vectorTypeList));
	os.write(reinterpret_cast<const char*>(&extras),sizeof(SaveOptions::extrasList));
	os.write(reinterpret_cast<const char*>(&column),sizeof(uint));
	os.write(reinterpret_cast<const char*>(&zeroModes),sizeof(uint));
	paramsIn.writeBinary(os);
	paramsOut.writeBinary(os);
	os.write(reinterpret_cast<const char*>(&printMessage),sizeof(bool));
	return os;
}

// readBinary
istream& SaveOptions::readBinary(istream& is) {
	is.read(reinterpret_cast<char*>(&printType),sizeof(SaveOptions::printTypeList));
	is.read(reinterpret_cast<char*>(&vectorType),sizeof(SaveOptions::vectorTypeList));
	is.read(reinterpret_cast<char*>(&extras),sizeof(SaveOptions::extrasList));
	is.read(reinterpret_cast<char*>(&column),sizeof(uint));
	is.read(reinterpret_cast<char*>(&zeroModes),sizeof(uint));
	paramsIn.readBinary(is);
	paramsOut.readBinary(is);
	is.read(reinterpret_cast<char*>(&printMessage),sizeof(bool));
	return is;
}	

/*-------------------------------------------------------------------------------------------------------------------------
	2. save
		- vec
			- static saveVec
			- static saveVecBinary
			- static saveVecSimpleAppend
			- static saveVecB
			- static saveVec
		- cVec
			- static savecVecSimple
			- static savecVecBinary
			- static savecVecSimpleAppend
			- static savecVecB
			- static savecVec
		- mat
		- cMat
		- spMat
-------------------------------------------------------------------------------------------------------------------------*/

// save vec - saveVecSimple
static void saveVecSimple(const string& f, const SaveOptions& opts, const vec& v) {
	fstream os;
	os.open(f.c_str(), ios::out);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	os.precision(16);
	os << left;
	uint length = v.size();
	if (opts.vectorType==SaveOptions::complex && !length%2) length = (uint)(length/2);
	for (uint j=0; j<length; j++) {
		if (opts.extras==SaveOptions::loc)			os << setw(25) << j;
		if (opts.vectorType==SaveOptions::complex)	os << setw(25) << v(2*j) << setw(25) << v(2*j+1);
		else										os << setw(25) << v(j);
													os << endl;
	}
	os.close();
}

// save vec - saveVecBinary
static void saveVecBinary(const string& f, const SaveOptions& opts,  const vec& v) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	const double* r;
	if (os.good()) {
		opts.writeBinary(os);
		for (uint j=0; j<v.size(); j++) {
			r = &v(j);
			os.write(reinterpret_cast<const char*>(r),sizeof(double));
		}
		os.close();
	}
	else {
		cerr << "save error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// save vec - simpleVecAppend
static void saveVecSimpleAppend(const string& f, const SaveOptions& opts, const vec& v) {
	unsigned int lengthOs = v.size();
	unsigned int lengthIs = countLines(f);
	ifstream is;
	is.open(f.c_str(),ios::in);
	if (!is.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	ofstream os;
	os.open("data/tempAppend",ios::out);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	os.precision(16);
	os << left;
	if (lengthOs!=lengthIs) {
		cerr << "save error: length of vector("<< lengthOs << ") to append not equal to file length("<< lengthIs << ")" << endl;
	}
	else {
		string lineIn;
		for (unsigned int j=0; j<lengthOs; j++){
		getline(is,lineIn);
		os << lineIn << setw(25) << v(j) << endl;
		}
		is.close();
		os.close();
		copyFile("data/tempAppend",f);
	}
}	

// save vec - saveVecB
static void saveVecB (const string& f, const SaveOptions& opts, const vec& v) {
	Parameters pin = opts.paramsIn, pout = opts.paramsOut;
	vec vo;
	if ((pin.N!=pout.N && pout.N!=0) || (pin.Nb!=pout.Nb && pout.Nb!=0)) {
		if (opts.vectorType==SaveOptions::realB)
			vo = interpolateReal(v,pin,pout);
		else if (opts.vectorType==SaveOptions::complexB)
			vo = interpolate(v,pin,pout);
		else {
			cerr << "save error: print vectorType option(" << opts.vectorType << ") not possible" << endl;
			return;
		}
	}
	else {
		vo = v;
	}
	
	fstream F;
	F.open(f.c_str(), ios::out);
	if (!F.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	uint x0 = intCoord(0,1,pout.Nb);
	F.precision(16);
	F << left;
	for (lint j=0; j<pout.N*pout.Nb; j++) {
		uint x = intCoord(j,1,pout.Nb);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		switch(opts.extras) {
			case SaveOptions::none:		break;
			case SaveOptions::loc: 		F << setw(25) << j;
										break;
			case SaveOptions::coords:	F << setw(25) << real(coordB(j,0,pout)) << setw(25) << imag(coordB(j,0,pout));
										F << setw(25) << real(coordB(j,1,pout)); 
										break;
			default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
										break;
		}
		switch(opts.vectorType) {
			case SaveOptions::realB:	F << setw(25) << vo(j) << endl;
										break;
			case SaveOptions::complexB:	F << setw(25) << vo(2*j) << setw(25) << vo(2*j+1)  << endl;
										break;
			default:					cerr << "save error: print vectorType option(" << opts.vectorType << ") not possible" << endl;
										return;
										break;
		}
	}
	if (vo.size()>2*pout.N*pout.Nb) {
		F << endl;
		for (unsigned int k=0; k<(vo.size()-2*pout.N*pout.Nb);k++) {
			F << setw(25) << vo(2*pout.N*pout.Nb+k) << endl;
		}
	}
	F.close();
}

// save vec - saveVec
static void saveVec(const string& f, const SaveOptions& opts, const vec& v) {
	Parameters pin = opts.paramsIn, pout = opts.paramsOut;
	vec vo;
	if ((pin.N!=pout.N && pout.N!=0) || (pin.NT!=pout.NT && pout.NT!=0)) {
		if (SaveOptions::real)
			vo = interpolateReal(v,pin,pout);
		else if (SaveOptions::complex)
			vo = interpolate(v,pin,pout);
	}
	else {
		vo = v;
	}
	fstream F;
	F.open(f.c_str(), ios::out);
	if (!F.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	unsigned int x0 = intCoord(0,1,pout.NT);
	F.precision(16);
	F << left;
	for (unsigned long int j=0; j<pout.N*pout.NT; j++) {
		unsigned int x = intCoord(j,1,pout.NT);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		switch(opts.extras) {
			case SaveOptions::none:		break;
			case SaveOptions::loc: 		F << setw(25) << j;
										break;
			case SaveOptions::coords:	F << setw(25) << real(coord(j,0,pout)) << setw(25) << imag(coord(j,0,pout));
										F << setw(25) << real(coord(j,1,pout)); 
										break;
			default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
										break;
		}
		switch(opts.vectorType) {
			case SaveOptions::real:		F << setw(25) << vo(j) << endl;
										break;
			case SaveOptions::complex:	F << setw(25) << vo(2*j) << setw(25) << vo(2*j+1)  << endl;
										break;
			default:					cerr << "save error: print vectorType option(" << opts.vectorType << ") not possible" << endl;
										break;
		}
	}
	if (vo.size()>2*pout.N*pout.NT) {
		F << endl;
		for (unsigned int j=0; j<(vo.size()-2*pout.N*pout.NT);j++) {
			F << setw(25) << vo(2*pout.N*pout.NT+j) << endl;
		}
	}
	F.close();
}

// save vec
void save(const string& f, const SaveOptions& opts, const vec& v) {
	if (!opts.good()) {
		cerr << "save error: SaveOptions not good, could not save " << f << endl;
		cerr << "SaveOptions: " << endl << opts << endl;
		return;
	}
	Filename F(f);
	if ((F.Suffix).compare(".dat")==0 && opts.printType==SaveOptions::ascii) {
		switch(opts.vectorType) {
				case SaveOptions::simple:	saveVecSimple(f,opts,v);
											break;
				case SaveOptions::real:		saveVec(f,opts,v);
											break;
				case SaveOptions::complex:	saveVec(f,opts,v);
											break;
				case SaveOptions::realB:	saveVecB(f,opts,v);
											break;
				case SaveOptions::complexB:	saveVecB(f,opts,v);
											break;
				case SaveOptions::append:	saveVecSimpleAppend(f,opts,v);
											break;
				default:					cerr << "save error: print vectorType option(" << opts.vectorType;
											cerr << ") not possible" << endl;
											break;
			}
	}
	else if ((F.Suffix).compare(".data")==0 && opts.printType==SaveOptions::binary) {
		saveVecBinary(f,opts,v);
	}
	else {
		cerr << "save error: printType option(" << opts.printType << ") not possible" << endl;
		cerr << "for file " << F << endl;
		return;
	}
	if (opts.printMessage) {
		printf("%12s%30s\n","saved: ",f.c_str());
	}
}

// save cVec - saveVecSimple
static void savecVecSimple(const string& f, const SaveOptions& opts, const cVec& v) {
	fstream F;
	F.open((f).c_str(), ios::out);
	if (!F.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	F.precision(16);
	F << left;
	uint length = v.size();
	for (uint j=0; j<length; j++) {
		if (opts.extras==SaveOptions::loc) 	F << setw(25) << j;
											F << setw(25) << real(v(j)) << setw(25) << imag(v(j)) << endl;
	}
	F.close();
}

// save cVec - savecVecBinary
static void savecVecBinary(const string& f, const SaveOptions& opts, const cVec& v) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	const comp* c;
	if (os.good()) {
		opts.writeBinary(os);
		for (uint j=0; j<v.size(); j++) {
			c = &v(j);
			os.write(reinterpret_cast<const char*>(c),sizeof(comp));
		}
		os.close();
	}
	else {
		cerr << "save error: cannot write to " << f << endl;
	}
}

// save cVec - simplecVecAppend
static void savecVecSimpleAppend(const string& f, const SaveOptions& opts, const cVec& v) {
	ifstream is;
	is.open(f.c_str(),ios::in);
	if (!is.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	ofstream os;
	os.open("data/tempAppend",ios::out);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	os.precision(16);
	os << left;
	unsigned int lengthOs = v.size();
	unsigned int lengthIs = countLines(f);
	if (lengthOs!=lengthIs) {
		cerr << "save error: length of vector("<< lengthOs << ") to append not equal to file length("<< lengthIs << ")" << endl;
		cerr << "in file " << f << endl;
		return;
	}
	else {
		string lineIn;
		for (unsigned int j=0; j<lengthOs; j++){
		getline(is,lineIn);
		os << lineIn << setw(25) << real(v(j)) << setw(25) << imag(v(j)) << endl;
		}
		is.close();
		os.close();
		copyFile("data/tempAppend",f);
	}
}	

// save cVec - saveVecB
static void savecVecB (const string& f, const SaveOptions& opts, const cVec& v) {
	Parameters pin = opts.paramsIn, pout = opts.paramsOut;
	cVec vo;
	if ((pin.N!=pout.N && pout.N!=0) || (pin.Nb!=pout.Nb && pout.Nb!=0)) {
		vo = interpolate(v,pin,pout);
	}
	else {
		vo = v;
	}
	fstream os;
	os.open(f.c_str(), ios::out);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	uint x0 = intCoord(0,1,pout.Nb);
	os.precision(16);
	os << left;
	for (lint j=0; j<pout.N*pout.Nb; j++) {
		uint x = intCoord(j,1,pout.Nb);
		if (x!=x0) { //this is put in for gnuplot
			os << endl;
			x0 = x;
		}
		switch(opts.extras) {
			case SaveOptions::none:		break;
			case SaveOptions::loc: 		os << setw(25) << j;
										break;
			case SaveOptions::coords:	os << setw(25) << real(coordB(j,0,pout)) << setw(25) << imag(coordB(j,0,pout));
										os << setw(25) << real(coordB(j,1,pout)); 
										break;
			default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
										break;
		}
		os << setw(25) << real(vo(j)) << setw(25) << imag(vo(j))  << endl;
	}
	if (vo.size()>pout.N*pout.Nb) {
		os << endl;
		for (unsigned int k=0; k<(vo.size()-pout.N*pout.Nb);k++) {
			os << setw(25) << vo(pout.N*pout.Nb+k) << endl;
		}
	}
	os.close();
}

// save cVec - saveVec
static void savecVec(const string& f, const SaveOptions& opts, const cVec& v) {
	Parameters pin = opts.paramsIn, pout = opts.paramsOut;
	cVec vo;
	if ((pin.N!=pout.N && pout.N!=0) || (pin.NT!=pout.NT && pout.NT!=0)) {
		vo = interpolate(v,pin,pout);
	}
	else {
		vo = v;
	}
	fstream os;
	os.open(f.c_str(), ios::out);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	uint x0 = intCoord(0,1,pout.NT);
	os.precision(16);
	os << left;
	for (lint j=0; j<pout.N*pout.NT; j++) {
		uint x = intCoord(j,1,pout.NT);
		if (x!=x0) { //this is put in for gnuplot
			os << endl;
			x0 = x;
		}
		switch(opts.extras) {
			case SaveOptions::none:		break;
			case SaveOptions::loc: 		os << setw(25) << j;
										break;
			case SaveOptions::coords:	os << setw(25) << real(coord(j,0,pout)) << setw(25) << imag(coord(j,0,pout));
										os << setw(25) << real(coord(j,1,pout)); 
										break;
			default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
										break;
		}
		os << setw(25) << real(vo(j)) << setw(25) << imag(vo(j))  << endl;
	}
	if (vo.size()>pout.N*pout.NT) {
		os << endl;
		for (uint j=0; j<(vo.size()-pout.N*pout.NT);j++) {
			os << setw(25) << vo(pout.N*pout.NT+j) << endl;
		}
	}
	os.close();
}

// save cVec
void save(const string& f, const SaveOptions& opts, const cVec& v) {
	if (!opts.good()) {
		cerr << "save error: SaveOptions not good, could not save " << f << endl;
		cerr << "SaveOptions: " << endl << opts << endl;
		return;
	}
	Filename F(f);
	if ((F.Suffix).compare(".dat")==0 && opts.printType==SaveOptions::ascii) {
		switch(opts.vectorType) {
				case SaveOptions::simple:	savecVecSimple(f,opts,v);
											break;
				case SaveOptions::real:		savecVec(f,opts,v);
											break;
				case SaveOptions::complex:	savecVec(f,opts,v);
											break;
				case SaveOptions::realB:	savecVecB(f,opts,v);
											break;
				case SaveOptions::complexB:	savecVecB(f,opts,v);
											break;
				case SaveOptions::append:	savecVecSimpleAppend(f,opts,v);
											break;
				default:					cerr << "save error: print vectorType option(" << opts.vectorType;
											cerr << ") not possible" << endl;
											break;
			}
	}
	else if ((F.Suffix).compare(".data")==0 && opts.printType==SaveOptions::binary) {
		savecVecBinary(f,opts,v);
	}
	else {
		cerr << "save error: printType option(" << opts.printType << ") not possible" << endl;
		cerr << "for file " << f << endl;
		return;
	}
	if (opts.printMessage) {
		printf("%12s%30s\n","saved: ",f.c_str());
	}
}

// save mat - binary
static void saveMatBinary(const string& f, const SaveOptions& opts, const mat& m) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	const double* d;
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			d = &m(j,k);
			os.write(reinterpret_cast<const char*>(d),sizeof(double));
		}
	}
	os.close();
}

// save mat - ascii
static void saveMatAscii(const string& f, const SaveOptions& opts, const mat& m) {
	fstream F;
	F.open(f.c_str(), ios::out);
	if (!F.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	F << left;
	F.precision(16);
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			if (opts.extras==SaveOptions::loc)
				F << setw(25) << j << setw(25) << k;
			F << setw(25) << m(j,k) << endl;
		}
	}
	F.close();
}

// save mat
void save(const string& f, const SaveOptions& opts, const mat& m) {
	if (!opts.good()) {
		cerr << "save error: SaveOptions not good, could not save " << f << endl;
		cerr << "SaveOptions: " << endl << opts << endl;
		return;
	}
	Filename F(f);
	if ((F.Suffix).compare(".dat")==0 && opts.printType==SaveOptions::ascii) {
		saveMatAscii(f,opts,m);
	}
	else if ((F.Suffix).compare(".data")==0 && opts.printType==SaveOptions::binary) {
		saveMatBinary(f,opts,m);
	}
	else {
		cerr << "save mat error: printType option(" << opts.printType << ") not possible" << endl;
		cerr << "for file " << f << endl;
		return;
	}
	if (opts.printMessage) {
		printf("%12s%30s\n","saved: ",f.c_str());
	}
}

// save cMat
static void savecMatBinary(const string& f, const SaveOptions& opts, const cMat& m) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	const comp* c;
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			c = &m(j,k);
			os.write(reinterpret_cast<const char*>(c),sizeof(comp));
		}
	}
	os.close();
}

// save cMat
void savecMatAscii(const string& f, const SaveOptions& opts, const cMat& m) {
	fstream F;
	F.open(f.c_str(), ios::out);
	if (!F.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	F << left;
	F.precision(16);
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			F << setw(25) << j << setw(25) << k;
			F << setw(25) << real(m(j,k)) << setw(25) << imag(m(j,k)) << endl;
		}
	}
	F.close();
}

// save cMat
void save(const string& f, const SaveOptions& opts, const cMat& m) {
	if (!opts.good()) {
		cerr << "save error: SaveOptions not good, could not save " << f << endl;
		cerr << "SaveOptions: " << endl << opts << endl;
		return;
	}
	Filename F(f);
	if ((F.Suffix).compare(".dat")==0 && opts.printType==SaveOptions::ascii) {
		savecMatAscii(f,opts,m);
	}
	else if ((F.Suffix).compare(".data")==0 && opts.printType==SaveOptions::binary) {
		savecMatBinary(f,opts,m);
	}
	else {
		cerr << "save cMat error: printType " << opts.printType << " not recognised" << endl;
		cerr << "for file " << f << endl;
		return;
	}
	if (opts.printMessage) {
		printf("%12s%30s\n","saved: ",f.c_str());
	}
}

// save spMat
void save(const string& f, const SaveOptions& opts, const spMat& m) {
	if (opts.printType!=SaveOptions::ascii) {
		cerr << "save spMat error: printType " << opts.printType << " not available" << endl;
		return;
	}
	if (!opts.good()) {
		cerr << "save error: SaveOptions not good, could not save " << f << endl;
		cerr << "SaveOptions: " << endl << opts << endl;
		return;
	}
	fstream F;
	F.open(f.c_str(), ios::out);
	if (!F.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	F << left;
	F.precision(16);
	for (int l=0; l<m.outerSize(); ++l) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(m,l); it; ++it) {
			F << setw(25) << it.row()+1 << setw(25) << it.col()+1;
			F << setw(25) << it.value() << endl;
		}
	}
	F.close();
	if (opts.printMessage) {
		printf("%12s%30s\n","saved: ",f.c_str());
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. load
		- vec
		- cVec
		- mat
		- cMat
		- spMat
-------------------------------------------------------------------------------------------------------------------------*/

// loadVecBinary
static void loadVecBinary(const string& f, SaveOptions& opts, vec& v) {
	ifstream is;
	is.open(f.c_str(),ios::binary);
	double d;
	if (is.good()) {
		opts.readBinary(is);
	}
	else {
		cerr << "load error: cannot read from " << f << endl;
		return;
	}
	uint pos = is.tellg();
	int lines = -1; // for some reason we should start on -1 not 0, see testBinaryPrint for verification
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&d),sizeof(double));
		lines++;
	}
	if (lines==-1) {
		cerr << "load error: no lines in file " << f << endl;
		return;
	}
	v = Eigen::VectorXd::Zero(lines);
	is.clear();
	is.seekg(pos);
	for (int j=0; j<lines; j++) {
		is.read(reinterpret_cast<char*>(&d),sizeof(double));
		v(j) = d;
	}
	is.close();
}


// load vec
void load(const string& f, SaveOptions& opts, vec& v) {
	Filename F(f);
	if (opts.printType==SaveOptions::binary && (F.Suffix).compare(".data")==0) {
		loadVecBinary(f,opts,v);
	}
	else if (opts.printType==SaveOptions::ascii && (F.Suffix).compare(".dat")==0) {
		uint col = opts.column;
		if (col==0) {
			switch(opts.extras) {
				case SaveOptions::none:		col = 1;
											break;
				case SaveOptions::loc:		col = 2;
											break;
				case SaveOptions::coords:	col = 4;
											break;
				default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
											return;
											break;
			}
		}
		uint fileLength = countLines(f);
		uint vLength;
		switch(opts.vectorType) {
			case SaveOptions::simple:	vLength = fileLength;
										break;
			case SaveOptions::real:		vLength = fileLength;
										break;
			case SaveOptions::complex:	vLength = 2*(fileLength-opts.zeroModes)+opts.zeroModes;
										break;
			case SaveOptions::realB:	vLength = fileLength;
										break;
			case SaveOptions::complexB:	vLength = 2*(fileLength-opts.zeroModes)+opts.zeroModes;
										break;
			default:					cerr << "save error: print vectorType option(" << opts.vectorType << ") not possible" << endl;
										return;
										break;
		}	
		vec vf = Eigen::VectorXd::Zero(vLength);
		fstream F;
		F.open(f.c_str(), ios::in);
		if (!F.good()) {
			cerr << "save error: stream not good for " << f << endl;
			return;
		}
		string line, temp;
		unsigned int j=0;
		while (getline(F, line)) {
			if (!line.empty()) {
				istringstream ss(line);
				if (col>1) for (unsigned int l=0; l<(col-1); l++) ss >> temp;
				if (opts.vectorType==SaveOptions::complex || opts.vectorType==SaveOptions::complexB) {
					if (j>=(fileLength-opts.zeroModes)) 	ss >> vf(fileLength-opts.zeroModes+j);
					else 									ss >> vf(2*j) >> vf(2*j+1);
				
				}
				else ss >> vf(j);
				j++;
			}
		}
		F.close();
	
		uint N1 = (opts.paramsIn).N, N2 = (opts.paramsOut).N;
		if (opts.vectorType!=SaveOptions::simple) {
			uint NT1 = (opts.paramsIn).NT, NT2 = (opts.paramsOut).NT;
			uint Nb1 = (opts.paramsIn).Nb, Nb2  = (opts.paramsOut).Nb;
			if (opts.vectorType==SaveOptions::complex && N1>0 && N2>0 && NT1>0 && NT2>0 && (N1!=N2 || NT1!=NT2)) {
				v = interpolate(vf,opts.paramsIn,opts.paramsOut);
				if (v.size()<(2*NT2*N2+opts.zeroModes)) {
					v.conservativeResize(2*NT2*N2+opts.zeroModes);
					for (uint zModes=0; zModes<opts.zeroModes; zModes++) {
						if (abs(v(2*NT2*N2+zModes))<MIN_NUMBER)  v(2*NT2*N2+zModes) = 0.5;
					}
				}
			}
			else if (opts.vectorType==SaveOptions::real && N1>0 && N2>0 && NT1>0 && NT2>0 && (N1!=N2 || NT1!=NT2)) {
				v = interpolateReal(vf,opts.paramsIn,opts.paramsOut);
				if (v.size()<(NT2*N2+opts.zeroModes)) {
					v.conservativeResize(NT2*N2+opts.zeroModes);
					for (uint zModes=0; zModes<opts.zeroModes; zModes++) {
						if (abs(v(NT2*N2+zModes))<MIN_NUMBER)  v(NT2*N2+zModes) = 0.5;
					}
				}
			}
			else if (opts.vectorType==SaveOptions::complexB && N1>0 && N2>0 && Nb1>0 && Nb2>0&& (N1!=N2 || Nb1!=Nb2)) {
				v = interpolate(vf,opts.paramsIn,opts.paramsOut);
				if (v.size()<(2*Nb2*N2+opts.zeroModes)) {
					v.conservativeResize(2*Nb2*N2+opts.zeroModes);
					for (uint zModes=0; zModes<opts.zeroModes; zModes++) {
						if (abs(v(2*Nb2*N2+zModes))<MIN_NUMBER)  v(2*Nb2*N2+zModes) = 0.5;
					}
				}
			}
			else if (opts.vectorType==SaveOptions::realB && N1>0 && N2>0 && Nb1>0 && Nb2>0 && (N1!=N2 || Nb1!=Nb2)) {
				v = interpolateReal(vf,opts.paramsIn,opts.paramsOut);
				if (v.size()<(Nb2*N2+opts.zeroModes)) {
					v.conservativeResize(Nb2*N2+opts.zeroModes);
					for (uint zModes=0; zModes<opts.zeroModes; zModes++) {
						if (abs(v(Nb2*N2+zModes))<MIN_NUMBER)  v(Nb2*N2+zModes) = 0.5;
					}
				}
			}
			else {
				v = vf;
			}
		}
		else {
			if (N1!=N2 && N2>0) {
				v = interpolate1d(vf,vf.size(),N2);
			}
			else {
				v = vf;
			}
		}
		if (opts.printMessage) {
			printf("%12s%30s\n","loaded: ",f.c_str());
		}
	}
	else {
		cerr << "load error: printType " << opts.printType << " not recognised" << endl;
		cerr << "for file " << f << endl;
		return;
	}
}

// loadcVecBinary
static void loadcVecBinary(const string& f, SaveOptions& opts, cVec& v) {
	ifstream is;
	is.open(f.c_str(),ios::binary);
	comp c;
	if (is.good()) {
		opts.readBinary(is);
	}
	else {
		cerr << "cannot read from " << f << endl;
		return;
	}
	uint pos = is.tellg();
	int lines = -1; // for some reason we should start on -1 not 0, see testBinaryPrint for verification
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&c),sizeof(comp));
		lines++;
	}
	if (lines==-1) {
		cerr << "load error: no lines in file " << f << endl;
		return;
	}
	is.clear();
	is.seekg(pos);
	v = Eigen::VectorXcd::Zero(lines);
	for (int j=0; j<lines; j++) {
		is.read(reinterpret_cast<char*>(&c),sizeof(comp));
		v(j) = c;
	}
	is.close();
}


// load cVec
void load(const string& f, SaveOptions& opts, cVec& v) {
	Filename F(f);
	if (opts.printType==SaveOptions::binary && (F.Suffix).compare(".data")==0) {
		loadcVecBinary(f,opts,v);
	}
	else if (opts.printType==SaveOptions::ascii && (F.Suffix).compare(".dat")==0) {
		uint col = opts.column;
		if (col==0) {
			switch(opts.extras) {
				case SaveOptions::none:		col = 1;
											break;
				case SaveOptions::loc:		col = 2;
											break;
				case SaveOptions::coords:	col = 4;
											break;
				default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
											break;
			}
		}
		uint fileLength = countLines(f);
		uint vLength = fileLength;
		cVec vf = Eigen::VectorXcd::Zero(vLength);
		fstream F;
		F.open(f.c_str(), ios::in);
		if (!F.good()) {
			cerr << "save error: stream not good for " << f << endl;
			return;
		}
		string line, temp;
		uint j=0;
		double realPart, imagPart;
		while (getline(F, line)) {
			if (!line.empty()) {
				istringstream ss(line);
				if (col>1) for (unsigned int l=0; l<(col-1); l++) ss >> temp;
				ss >> realPart >> imagPart;
				vf(j) = realPart + comp(0.0,1.0)*imagPart;
				j++;
			}
		}
		F.close();
	
		uint N1 = (opts.paramsIn).N, N2 = (opts.paramsOut).N;
		if (opts.vectorType!=SaveOptions::simple) {
			uint NT1 = (opts.paramsIn).NT, NT2 = (opts.paramsOut).NT;
			uint Nb1 = (opts.paramsIn).Nb, Nb2  = (opts.paramsOut).Nb;
			if (opts.vectorType==SaveOptions::complex || opts.vectorType==SaveOptions::real) {
				if ((N1!=N2 || NT1!=NT2) && N1>0 && N2>0 && NT1>0 && NT2>0) {
					v = interpolate(vf,opts.paramsIn,opts.paramsOut);
				}
				else {
					v = vf;
				}
			}
			else if (opts.vectorType==SaveOptions::complexB || opts.vectorType==SaveOptions::realB) {
				if ((N1!=N2 || Nb1!=Nb2) && N1>0 && N2>0 && Nb1>0 && Nb2>0) {
					v = interpolate(vf,opts.paramsIn,opts.paramsOut);
				}
				else {
					v = vf;
				}
			}
			else {
				v = vf;
			}
		}
		else {
			if (N1!=N2 && N1>0 && N2>0) {
				v = interpolate1d(vf,N1,N2);
			}
			else {
				v = vf;
			}
		}
		if (opts.printMessage) {
			printf("%12s%30s\n","loaded: ",f.c_str());
		}
	}
	else {
		cerr << "load error: printType " << opts.printType << " not recognised" << endl;
		cerr << "for file " << f << endl;
		return;
	}
}

// load mat  - binary - only works for square matrices
void loadMatBinary(const string& f, SaveOptions& opts, mat& m) {
	Filename F(f);
	if ((F.Suffix).compare(".data")!=0) {
		cerr << "load error: printType " << opts.printType << " not recognised" << endl;
		cerr << "for file " << f << endl;
		return;
	}
	ifstream is;
	is.open(f.c_str(),ios::binary);
	if (!is.good()) {
		cerr << "cannot read from " << f << endl;
	}
	uint pos = is.tellg();
	int lines = -1;
	double d;
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&d),sizeof(double));
		lines++;
	}
	if (lines==-1) {
		cerr << "load error: no lines in file " << f << endl;
		return;
	}
	is.clear();
	is.seekg(pos);
	uint rows = (uint)sqrt(lines);
	if (abs((double)rows-sqrt(lines))>MIN_NUMBER*1.0e2) {
		cerr << "load mat error: matrix in " << f << " not square" << endl; 
	}
	m = Eigen::MatrixXd::Zero(rows,rows);
	for (uint j=0; j<rows; j++) {
		for (uint k=0; k<rows; k++) {
			is.read(reinterpret_cast<char*>(&d),sizeof(double));
			m(j,k) = d;
		}
	}
	is.close();
}

// load mat - ascii- assumes square matrix
void loadMatAscii(const string& f, SaveOptions& opts, mat& m) {
	Filename F(f);
	if ((F.Suffix).compare(".dat")!=0) {
		cerr << "load error: printType " << opts.printType << " not recognised" << endl;
		cerr << "for file " << f << endl;
		return;
	}
	uint rowsF = countLines(f), rows;
	uint columnsF = countColumns(f), columns;
	if (opts.extras==SaveOptions::loc && columnsF==3) {
		rows = (uint)sqrt(rowsF); // this option could certainly be improved to allow for non-square matrices if required
		columns = (uint)sqrt(rowsF);
	}
	else if (opts.extras==SaveOptions::none && columnsF==1) {
		rows = (uint)sqrt(rowsF);
		columns = (uint)sqrt(rowsF);
	}
	else if (opts.extras==SaveOptions::none) {
		rows = rowsF;
		columns = columnsF;
	}
	else {
		cerr << "SaveOptions::extras choice not allowed: " << opts.extras << endl;
		return;
	}
	m = Eigen::MatrixXd::Zero(rows,columns);
	fstream is;
	is.open(f.c_str(), ios::in);
	string line;
	uint j = 0, k = 0;
	double v;
	while (getline(is, line)) {
		if (!line.empty()) {
			istringstream ss(line);
			if (opts.extras==SaveOptions::loc && columnsF==3) {
				ss >> j >> k >> v;
				if (j>rows || k>columns)
					cerr << "load error: matrix index(" << j << "," << k << ") > " << columns << endl;
				m(j,k) = v;
			}
			else if (opts.extras==SaveOptions::none && columnsF==1) {
				ss >> v;
				m(j,k) = v;
				k++;
				if (k==columns) {
					k = 0;
					j++;
				}
			}
			else if (opts.extras==SaveOptions::none) {
				while (ss >> v) {
					m(j,k) = v;
					k++;
					if (k==columns) {
						k = 0;
						j++;
					}
				}
			}
		}
	}
	is.close();
	if (opts.printMessage) {
		printf("%12s%30s\n","loaded: ",f.c_str());
	}
}

// load mat
void load(const string& f, SaveOptions& opts, mat& m) {
	switch(opts.printType) {
		case SaveOptions::binary:
									loadMatBinary(f,opts,m);
									break;
		case SaveOptions::ascii:	
									loadMatAscii(f,opts,m);
									break;
		default:
									cerr << "load mat error: printType " << opts.printType << " not recognised" << endl;
									break;
		}
}

// load cMat - binary - only works for square matrices
void loadcMatBinary(const string& f, SaveOptions& opts, cMat& m) {
	Filename F(f);
	if ((F.Suffix).compare(".data")!=0) {
		cerr << "load error: printType " << opts.printType << " not recognised" << endl;
		cerr << "for file " << f << endl;
		return;
	}
	ifstream is;
	is.open(f.c_str(),ios::binary);
	if (!is.good()) {
		cerr << "cannot read from " << f << endl;
	}
	uint pos = is.tellg();
	int lines = -1;
	double d;
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&d),sizeof(double));
		lines++;
	}
	if (lines==-1) {
		cerr << "load error: no lines in file " << f << endl;
		return;
	}
	is.clear();
	is.seekg(pos);
	uint rows = (uint)sqrt(lines);
	if (abs((double)rows-sqrt(lines))>MIN_NUMBER*1.0e2) {
		cerr << "load mat error: matrix in " << f << " not square" << endl; 
	}
	m = Eigen::MatrixXcd::Zero(rows,rows);
	comp c;
	for (uint j=0; j<rows; j++) {
		for (uint k=0; k<rows; k++) {
			is.read(reinterpret_cast<char*>(&c),sizeof(comp));
			m(j,k) = c;
		}
	}
	is.close();
}

// load cMat - ascii
void loadcMatAscii(const string& f, SaveOptions& opts, cMat& m) {
	Filename F(f);
	if ((F.Suffix).compare(".dat")!=0) {
		cerr << "load error: printType " << opts.printType << " not recognised" << endl;
		cerr << "for file " << f << endl;
		return;
	}
	uint fileLength = countLines(f);
	uint matLength = (uint)sqrt(fileLength);
	fstream is;
	is.open(f.c_str(), ios::in);
	string line;
	uint j, k;
	double realPart, imagPart;
	m = Eigen::MatrixXcd::Zero(matLength,matLength);
	while (getline(is, line)) {
		if (!line.empty()) {
			istringstream ss(line);
			ss >> j >> k >> realPart >> imagPart;
			if (j>matLength || k>matLength) cerr << "load error: matrix index > sqrt(fileLength)" << endl;
			m(j,k) = realPart + comp(0.0,1.0)*imagPart;;
		}
	}
	is.close();
	if (opts.printMessage) {
		printf("%12s%30s\n","loaded: ",f.c_str());
	}
}

// load cMat
void load(const string& f, SaveOptions& opts, cMat& m) {
	switch(opts.printType) {
		case SaveOptions::binary:
									loadcMatBinary(f,opts,m);
									break;
		case SaveOptions::ascii:	
									loadcMatAscii(f,opts,m);
									break;
		default:
									cerr << "load cMat error: printType " << opts.printType << " not recognised" << endl;
									break;
		}
}

// load spMat
void load(const string& f , SaveOptions& opts, spMat& m) {
	uint fileLength = countLines(f);
	Eigen::VectorXi to_reserve(fileLength); //an overestimate
	to_reserve.setZero(fileLength);
	fstream F;
	F.open(f.c_str(), ios::in);
	string line;
	unsigned int nnz = 0, length = 0, count = 0, row, column;
	double value;	
	Eigen::VectorXi rowVec(fileLength), columnVec(fileLength);
	vec valueVec(fileLength);
	while (getline(F, line)) {
		if (!line.empty()) {
			istringstream ss(line);
			ss >> row >> column >> value;
			if (abs(value)>MIN_NUMBER) {
				rowVec(count) = row-1;
				columnVec(count) = column-1;
				valueVec(count) = value;
				to_reserve(row-1) = to_reserve(row-1) + 1;
				nnz++;
				count++;
				if (row>length) {
					length = row;
				}
			}
		}
	}
	if (nnz==0) {
		cerr << "loadSpMat failed, no data in file: " << f << endl;
		}
	to_reserve.conservativeResize(length);
	rowVec.conservativeResize(count);
	columnVec.conservativeResize(count);
	valueVec.conservativeResize(count);
	spMat M(length,length);
	M.setZero(); //just making sure
	M.reserve(to_reserve);
	for (unsigned int l=0;l<count;l++) {
		M.insert(rowVec(l),columnVec(l)) = valueVec(l);
	}
	M.makeCompressed();
	m = M;
	if (opts.printMessage) {
		printf("%12s%30s\n","loaded: ",f.c_str());
	}
}
/*-------------------------------------------------------------------------------------------------------------------------
	4. plot
-------------------------------------------------------------------------------------------------------------------------*/

// plot
void plot(const string& f, const PlotOptions& opts) {
	string style = ((opts.style).empty()? "linespoints": opts.style);
	string output = ((opts.output).empty()? "pics/pic.png": opts.output);
	uint col = (opts.column==0? 1: opts.column);
	string columns1 = (opts.column2==0? numberToString<uint>(col): (numberToString<uint>(col)+":"+numberToString<uint>(opts.column2)));
	string columns2 = (opts.column3==0? "": (numberToString<uint>(col)+":"+numberToString<uint>(opts.column3)));
	string columns3 = (opts.column4==0? "": (numberToString<uint>(col)+":"+numberToString<uint>(opts.column4)));
	if ((opts.gp).empty()) {		
		string commandOpenStr = "gnuplot -persistent";
		const char * commandOpen = commandOpenStr.c_str();
		FILE * gnuplotPipe = popen (commandOpen,"w");
		string command1Str = "set term png size 1600,800";
		string command2Str = "set output \""+output+"\"";
		string command3Str = "plot \"" + f + "\" using " + columns1 + " with " + style;
		if (!columns2.empty()) command3Str += ", \"" + f + "\" using " + columns2 + " with " + style;
		if (!columns3.empty()) command3Str += ", \"" + f + "\" using " + columns3 + " with " + style;
		string command4Str = "pause -1";
		if (output.compare("gui")!=0) {
			fprintf(gnuplotPipe, "%s \n",command1Str.c_str());
			fprintf(gnuplotPipe, "%s \n",command2Str.c_str());
		}
		fprintf(gnuplotPipe, "%s \n",command3Str.c_str());
		if (output.compare("gui")==0)
			fprintf(gnuplotPipe, "%s \n",command4Str.c_str());
		pclose(gnuplotPipe);
	}
	else {
		string commandStr = "gnuplot -e \"inFile='"+f+"'; outFile='"+output+"'; "+" style='"+style+"'\" "+"'"+opts.gp+"'";
		if (output.compare("gui")==0) commandStr += " -persistent";
		FILE * gnuplotPipe = popen(commandStr.c_str(),"w");
		fprintf(gnuplotPipe,"%s\n"," ");
		pclose(gnuplotPipe);
	}
	if (opts.printMessage) {
		printf("%12s%30s\n","plotted: ",output.c_str());
	}
}

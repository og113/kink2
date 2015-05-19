/*
	an attempt at printing binary rather than ascii files
*/

#include<iostream>
#include<fstream>
#include<string>
#include<complex>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include "parameters.h"
#include "print.h"
#include "simple.h"

using namespace std;

typedef complex<double> comp;
typedef unsigned int uint;
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXcd cVec;
typedef Eigen::MatrixXd mat;
typedef Eigen::MatrixXcd cMat;

#ifndef pi
#define pi 3.14159265359
#endif

// save vec
void saveBinary(const string& f, const vec& v) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	const double* d;
	for (uint j=0; j<v.size(); j++) {
		d = &v(j);
		os.write(reinterpret_cast<const char*>(d),sizeof(double));
	}
	os.close();
}

// save cVec
void saveBinary(const string& f, const cVec& v) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	const comp* c;
	for (uint j=0; j<v.size(); j++) {
		c = &v(j);
		os.write(reinterpret_cast<const char*>(c),sizeof(comp));
	}
	os.close();
}

// save mat
void saveBinary(const string& f, const mat& m) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	const double* d;
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			d = &m(j,k);
			os.write(reinterpret_cast<const char*>(d),sizeof(double));
		}
	}
	os.close();
}

// save cMat
void saveBinary(const string& f, const cMat& m) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	const comp* c;
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			c = &m(j,k);
			os.write(reinterpret_cast<const char*>(c),sizeof(comp));
		}
	}
	os.close();
}

// load vec
void loadBinary(const string& f, vec& v) {
	ifstream is;
	is.open(f.c_str(),ios::binary);
	uint lines = countDoubles(f);
	v = Eigen::VectorXd::Zero(lines);
	double d;
	for (uint j=0; j<lines; j++) {
		is.read(reinterpret_cast<char*>(&d),sizeof(double));
		v(j) = d;
	}
	is.close();
}

// load cVec
void loadBinary(const string& f, cVec& v) {
	ifstream is;
	is.open(f.c_str(),ios::binary);
	comp c;
	uint lines = countType(f,c);
	v = Eigen::VectorXcd::Zero(lines);
	for (uint j=0; j<lines; j++) {
		is.read(reinterpret_cast<char*>(&c),sizeof(comp));
		v(j) = c;
	}
	is.close();
}

// load mat, only works for square matrices
void loadBinary(const string& f, mat& m) {
	ifstream is;
	is.open(f.c_str(),ios::binary);
	uint lines = countDoubles(f);
	uint rows = (uint)sqrt(lines);
	if (abs((double)rows-sqrt(lines))>MIN_NUMBER*1.0e2) {
		cerr << "load mat error: matrix in " << f << " not square" << endl; 
	}
	m = Eigen::MatrixXd::Zero(rows,rows);
	double d;
	for (uint j=0; j<rows; j++) {
		for (uint k=0; k<rows; k++) {
			is.read(reinterpret_cast<char*>(&d),sizeof(double));
			m(j,k) = d;
		}
	}
	is.close();
}

// load cMat - binary, only works for square matrices
void loadBinary(const string& f, cMat& m) {
	ifstream is;
	is.open(f.c_str(),ios::binary);
	uint lines = countDoubles(f);
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

int main() {

cout << "test print binary: " << endl;

uint length = 4;

// vec
vec v, w;
v = Eigen::VectorXd::Random(length);

string f = "tests/data/v.dat";

saveBinary(f,v);

cout << "v saved to " << f << endl;

uint lines = countDoubles(f);

cout << lines << " lines in " << f << endl;

loadBinary(f,w);

cout << f << " loaded to w" << endl;

cout << "w = " << endl;
for (uint j=0; j<w.size(); j++) {
	cout << w[j] << endl;
}
cout << endl;


// cVec
cVec cv, cw;
cv = Eigen::VectorXcd::Random(length);

string cf = "tests/data/cv.dat";

saveBinary(cf,cv);

cout << "cv saved to " << cf << endl;

comp c;
lines = countType(cf,c);

cout << lines << " lines in " << cf << endl;

loadBinary(cf,cw);

cout << cf << " loaded to cw" << endl;

cout << "cw = " << endl;
for (uint j=0; j<cw.size(); j++) {
	cout << cw[j] << endl;
}
cout << endl;

// mat
mat m, n;
m = Eigen::MatrixXd::Random(length,length);

string mf = "tests/data/m.dat";

saveBinary(mf,m);

cout << "m saved to " << mf << endl;

double d;
lines = countType(mf,d);

cout << lines << " lines in " << mf << endl;

loadBinary(mf,n);

cout << mf << " loaded to n" << endl;

cout << "n = " << endl;
for (uint j=0; j<n.rows(); j++) {
	for (uint k=0; k<n.cols(); k++) {
		cout << n(j,k) << " ";
	}
	cout << endl;
}
cout << endl;

// cMat
cMat mc, nc;
mc = Eigen::MatrixXcd::Random(length,length);

string mcf = "tests/data/mc.dat";

saveBinary(mcf,mc);

cout << "mc saved to " << mcf << endl;

lines = countType(mcf,c);

cout << lines << " lines in " << mcf << endl;

loadBinary(mcf,nc);

cout << mcf << " loaded to nc" << endl;

cout << "nc = " << endl;
for (uint j=0; j<nc.rows(); j++) {
	for (uint k=0; k<nc.cols(); k++) {
		cout << nc(j,k) << " ";
	}
	cout << endl;
}
cout << endl;

// spMat - not yet done
cout << "n.b. no binary method for printing spMat yet written" << endl;

// parameters
cout << "saving parameters in binary ";;

PrimaryParameters p, q;
p.load("inputsP");

string pf = "tests/data/p.dat";
cout << "to " << pf << endl;
ofstream os;
os.open(pf.c_str(), ios::binary);
if (os.good()) {
	os.write(reinterpret_cast<char*>(&p),sizeof(p));
}
else {
	cerr << "cannot write to " << pf << endl;
}
os.close();

ifstream is;
is.open(pf.c_str(), ios::binary);
string dross;
if (is.good()) {
	is.read(reinterpret_cast<char*>(&q),sizeof(q));
}
else {
	cerr << "cannot read from " << pf << endl;
}
cout << q << endl;

cout << "parameters print test = " << 1-(p==q) << endl;

return 0;
}

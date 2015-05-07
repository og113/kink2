/*
	an attempt at printing binary rather than ascii files
*/

#include<iostream>
#include<fstream>
#include<string>
#include<Eigen/Dense>

using namespace std;

typedef Eigen::VectorXd vec;
typedef unsigned int uint;

#ifndef pi
#define pi 3.14159265359
#endif

uint countDoubles(const string& f) {
	uint lines = -1;
	ifstream is;
	is.open(f.c_str(),ios::binary);
	double dross;
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&dross),sizeof(double));
		lines++;
	}
	is.close();
	return lines;
}

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

int main() {

cout << "test print binary: " << endl;

uint length = 10;
vec v, w;
v = Eigen::VectorXd::Random(length);

string f = "tests/data/v.dat";

saveBinary(f,v);

cout << "v saved to " << f << endl;

uint lines = countDoubles(f);

cout << lines << " lines in f" << endl;

loadBinary(f,w);

cout << f << " loaded to w" << endl;

cout << "w = " << endl;
for (uint j=0; j<w.size(); j++) {
	cout << w[j] << endl;
}
cout << endl;

return 0;
}

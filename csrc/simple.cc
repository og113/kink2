/*-------------------------------------------------------------------------------------------------------------------------
definitions of some very simple functions and classes
-------------------------------------------------------------------------------------------------------------------------*/

#include <cstring>
#include <cstdlib> //for rand, srand
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <sys/time.h>
#include <Eigen/Dense>
#include <glob.h>
#include "simple.h"

using namespace std;

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
	10 - explicit instantiation
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. number to string and string to number
		- number to string
		- string to number
-------------------------------------------------------------------------------------------------------------------------*/

//to convert number to string, usage is string str = NumberToString<number type>(x);
template <class T>
string numberToString ( const T& Number ) {
	stringstream ss;
	ss << Number;
	return ss.str();
}
	
//shorthand version;
template <class T>
string nts ( const T& Number ) {
	stringstream ss;
	ss << Number;
	return ss.str();
}

//shorthand version;
template <class T>
string nts ( const T& Number, const uint& prec) {
	stringstream ss;
	ss << fixed << setprecision(prec) << Number;
	return ss.str();
}

//to convert string to number, usage is (number type) x = StringToNumber<number type>(str);
template <class T>
T stringToNumber ( const string& Text ) {
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}
	
//shorthand version;
template <class T>
T stn ( const string& Text ) {
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

// is number?
bool isNumber( const string& Text ) {
	return( strspn( Text.c_str(), "0123456789" ) == Text.size() );
}
	
/*-------------------------------------------------------------------------------------------------------------------------
	2. absDiff
		- double
		- comp
		- vec
		- cVec
		- mat
		- cMat
-------------------------------------------------------------------------------------------------------------------------*/
	
// absDiff double
double absDiff(const double& numA, const double& numB) {
	if (abs(numA)>MIN_NUMBER && abs(numB)>MIN_NUMBER) return 2.0*abs(numA-numB)/sqrt(numA*numA+numB*numB);
	else if (abs(numA)>MIN_NUMBER || abs(numB)>MIN_NUMBER) return abs(numA-numB);
	else return 0.0;
}

// absDiff comp
double absDiff (const comp& numA, const comp& numB) {
	if (abs(numA)>MIN_NUMBER && abs(numB)>MIN_NUMBER) return 2.0*abs(numA-numB)/sqrt(norm(numA)+norm(numB));
	else if (abs(numA)>MIN_NUMBER || abs(numB)>MIN_NUMBER) return abs(numA-numB);
	else return 0.0;
}

// absDiff Eigen::VectorXd
double absDiff(const vec& A, const vec& B) {
	double normSqrdA = A.squaredNorm();
	double normSqrdB = B.squaredNorm();
	vec diff = A-B;
	double normDiff = diff.norm();
	if (normSqrdA>MIN_NUMBER && normSqrdB>MIN_NUMBER ) return 2.0*normDiff/sqrt(normSqrdA+normSqrdB);
	else if (normSqrdA>MIN_NUMBER || normSqrdB>MIN_NUMBER ) return normDiff;
	else return 0.0;
}

// absDiff Eigen::VectorXcd
double absDiff(const cVec& A, const cVec& B) {
	double normSqrdA = A.squaredNorm();
	double normSqrdB = B.squaredNorm();
	cVec diff = A-B;
	double normDiff = diff.norm();
	if (normSqrdA>MIN_NUMBER && normSqrdB>MIN_NUMBER ) return 2.0*normDiff/sqrt(normSqrdA+normSqrdB);
	else if (normSqrdA>MIN_NUMBER || normSqrdB>MIN_NUMBER ) return normDiff;
	else return 0.0;
}

// absDiff Eigen::MatrixXd
double absDiff(const mat& A, const mat& B) {
	double normSqrdA = A.squaredNorm();
	double normSqrdB = B.squaredNorm();
	mat diff = A-B;
	double normDiff = diff.norm();
	if (normSqrdA>MIN_NUMBER && normSqrdB>MIN_NUMBER ) return 2.0*normDiff/sqrt(normSqrdA+normSqrdB);
	else if (normSqrdA>MIN_NUMBER || normSqrdB>MIN_NUMBER ) return normDiff;
	else return 0.0;
}

// absDiff Eigen::MatrixXcd
double absDiff(const cMat& A, const cMat& B) {
	double normSqrdA = A.squaredNorm();
	double normSqrdB = B.squaredNorm();
	cMat diff = A-B;
	double normDiff = diff.norm();
	if (normSqrdA>MIN_NUMBER && normSqrdB>MIN_NUMBER ) return 2.0*normDiff/sqrt(normSqrdA+normSqrdB);
	else if (normSqrdA>MIN_NUMBER || normSqrdB>MIN_NUMBER ) return normDiff;
	else return 0.0;
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. factorial
-------------------------------------------------------------------------------------------------------------------------*/

//factorial
int factorial(const int& f_input){
	int f_result = 1;
	for (int l=0; l<f_input; l++) f_result *= (l+1);
	return f_result;
}

/*-------------------------------------------------------------------------------------------------------------------------
	4. currentDateTime, currentPartSec
	
N.B. Visit http://en.cppreference.com/w/cpp/chrono/c/strftime for more information about date/time format
-------------------------------------------------------------------------------------------------------------------------*/

//getting the date and time
string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%y%m%d%H%M%S", &tstruct);
    return buf;
}

//getting the part of the second, a number between 0 and 1e6, divide by 1e6 to get the fraction
string currentPartSec() {
    timeval tim;
    gettimeofday(&tim, NULL);
    return numberToString<double>(tim.tv_usec);
}

/*-------------------------------------------------------------------------------------------------------------------------
	5. copyFile
-------------------------------------------------------------------------------------------------------------------------*/

//copy a file
void copyFile(const string & inputFile, const string & outputFile) {
	ifstream  is(inputFile.c_str(), ios::binary);
	ofstream  os(outputFile.c_str(), ios::binary);
	os << is.rdbuf();
}

/*-------------------------------------------------------------------------------------------------------------------------
	6. countLines, countColumns, countDoubles
-------------------------------------------------------------------------------------------------------------------------*/

// count lines
//count non-empty lines of a file
uint countLines(const string & file_to_count) {
	ifstream fin;
	fin.open(file_to_count.c_str());
	if (!fin.good()) cerr << "countLines error: " << file_to_count << " not opened properly." << endl;
	string line;
	unsigned int counter = 0;
	while(!fin.eof()) {
		getline(fin,line);
		if(line.empty()) continue;
		counter++;
	}		
	fin.close();
    return counter;
}

// countColumns
uint countColumns(const string & file_to_count) {
	ifstream fin;
	fin.open(file_to_count.c_str());
	if (!fin.good()) cerr << "countRows error: " << file_to_count << " not opened properly." << endl;
	string line;
	unsigned int counter = 0;
	while(!fin.eof()) {
		getline(fin,line);
		if(line.empty()) continue;
		else {
			stringstream ss(line);
			string temp;
			while (ss >> temp) {
				counter++;
			}
			break;
		}
	}		
	fin.close();
    return counter;
}

// countColumns
uint countColumns(const string & file_to_count, const char& sep) {
	ifstream fin;
	fin.open(file_to_count.c_str(), ios::in);
	if (!fin.good()) cerr << "countColumns error: " << file_to_count << " not opened properly." << endl;
	string line;
	unsigned int counter = 0;
	while(!fin.eof()) {
		getline(fin,line);
		if(line.empty()) continue;
		else {
			stringstream ss(line);
			string temp;
			while (getline(ss,temp,sep)) {
				counter++;
			}
			break;
		}
	}		
	fin.close();
    return counter;
}

// countDoubles
uint countDoubles(const string& f) {
	uint lines = -1; // for some reason we should start on -1 not 0, see testBinaryPrint for verification
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

// count comp (binary)
uint countComp(const string& f) {
	uint lines = -1; 
	ifstream is;
	is.open(f.c_str(),ios::binary);
	comp dross;
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&dross),sizeof(comp));
		lines++;
	}
	is.close();
	return lines;
}

// count type (binary)
template <class T>
uint countType(const string& f, const T& t) {
	uint lines = -1; 
	ifstream is;
	is.open(f.c_str(),ios::binary);
	T dross;
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&dross),sizeof(T));
		lines++;
	}
	is.close();
	return lines;
}

/*-------------------------------------------------------------------------------------------------------------------------
	7. smallestLoc
-------------------------------------------------------------------------------------------------------------------------*/

// smallestLoc
// function giving location of smallest element of a vector of type T
template <typename T>
uint smallestLoc(const vector<T>& inVector) {
	uint loc = 0;
	for(uint l=1;l<inVector.size();l++) {
		if (abs(inVector[l])<abs(inVector[loc])) {
			loc = l;
		}
	}
	return loc;
}

/*-------------------------------------------------------------------------------------------------------------------------
	8. randDouble

	n.b. must first seed rand with srand(time(NULL)) or something similar
-------------------------------------------------------------------------------------------------------------------------*/

// randDouble
double randDouble(const double& min, const double& max) {
	double f = (double)rand() / RAND_MAX;
    return min + f*(max - min);
}

/*-------------------------------------------------------------------------------------------------------------------------
	9. mod
-------------------------------------------------------------------------------------------------------------------------*/

// mod
double mod(const double& x, const double& min, const double& max) {
	double Min, Max;
	if (min<max) {
		Min = min;
		Max = max;
	}
	else if (min>max) {
		Min = max;
		Max = min;
	}
	else {
		cerr << "mod error, range of size zero" << endl;
		return 1.0;
	}
		
	if (x>=Min && x<=Max)
		return x;
	else if (x>Max) {
		int ranges = (int)((x-Min)/(Max-Min));
		return x-(double)ranges*(Max-Min);
	}
	else {
		int ranges = (int)((Max-x)/(Max-Min));
		return x+(double)ranges*(Max-Min);
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	10. delta
-------------------------------------------------------------------------------------------------------------------------*/

double delta(const uint& i, const uint& j) {
	return (i==j? 1.0: 0.0);
}

/*-------------------------------------------------------------------------------------------------------------------------
	11. sigFig
-------------------------------------------------------------------------------------------------------------------------*/

double sigFig(double num, double N) {
    double d = log10(num);
    double power;
    if (num > 0) {
        d = ceil(d);
        power = -(d-N);
    }
    else {
        d = floor(d); 
        power = -(d-N);
    }

    return (int)(num * pow(10.0, power) + 0.5) * pow(10.0, -power);
}

/*-------------------------------------------------------------------------------------------------------------------------
	12. sign
-------------------------------------------------------------------------------------------------------------------------*/

double sign(double num) {
	return (num>0.0? 1.0: -1.0);
}

/*-------------------------------------------------------------------------------------------------------------------------
	13. splitString
-------------------------------------------------------------------------------------------------------------------------*/

vector<string> splitString(const string & str, const string & delimiters) {
    vector<string> v;
    size_t start = 0;
    size_t pos = str.find_first_of(delimiters, start);
    while(pos != string::npos) {
        if(pos != start) // ignore empty tokens
            v.push_back(str.substr(start, pos - start));
        start = pos + 1;
        pos = str.find_first_of(delimiters, start);
    }
    if(start < str.length()) // ignore trailing delimiter
        v.push_back(str.substr(start,str.length() - start)); // add what's left of the string
    return v;
}


/*-------------------------------------------------------------------------------------------------------------------------
	14. explicit instantiation
		- numberToString, stringToNumber
		- smallestLoc
		- countType
-------------------------------------------------------------------------------------------------------------------------*/


template string numberToString<int>(const int&);
template string numberToString<uint>(const uint&);
template string numberToString<lint>(const lint&);
template string numberToString<long long unsigned>(const long long unsigned&);
template string numberToString<double>(const double&);
template string numberToString<comp>(const comp&);

template string nts<int>(const int&);
template string nts<uint>(const uint&);
template string nts<lint>(const lint&);
template string nts<long long unsigned>(const long long unsigned&);
template string nts<double>(const double&);
template string nts<double>(const double&, const uint&);
template string nts<comp>(const comp&);

template int stringToNumber<int>(const string&);
template uint stringToNumber<uint>(const string&);
template lint stringToNumber<lint>(const string&);
template long long unsigned stringToNumber<long long unsigned>(const string&);
template double stringToNumber<double>(const string&);
template comp stringToNumber<comp>(const string&);

template int stn<int>(const string&);
template uint stn<uint>(const string&);
template lint stn<lint>(const string&);
template long long unsigned stn<long long unsigned>(const string&);
template double stn<double>(const string&);
template comp stn<comp>(const string&);

template uint smallestLoc<int>(const vector<int>&);
template uint smallestLoc<double>(const vector<double>&);
template uint smallestLoc<comp>(const vector<comp>&);

template uint countType<int>(const string&,const int&);
template uint countType<uint>(const string&,const uint&);
template uint countType<double>(const string&,const double&);
template uint countType<comp>(const string&,const comp&);


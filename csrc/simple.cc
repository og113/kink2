/*-------------------------------------------------------------------------------------------------------------------------
definitions of some very simple functions and classes
-------------------------------------------------------------------------------------------------------------------------*/
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "simple.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - number to string and string to number
	2 - absDiff
	3 - factorial
	4 - currentDateTime
	5 - copyFile
	6 - countLines
	7 - smallestLoc
	8 - explicit instantiation
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. number to string and string to number
		- number to string
		- string to number
-------------------------------------------------------------------------------------------------------------------------*/

//to convert number to string, usage is string str = NumberToString<number type>(x);
template <class T>
string numberToString ( const T& Number )
	{
	stringstream ss;
	ss << Number;
	return ss.str();
	}

//to convert string to number, usage is (number type) x = StringToNumber<number type>(str);
template <class T>
T stringToNumber ( const string& Text )
	{
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
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
	if (abs(numA)>MIN_NUMBER || abs(numB)>MIN_NUMBER) return 2.0*abs(numA-numB)/sqrt(numA*numA+numB*numB);
	else return 0.0;
}

// absDiff comp
double absDiff (const comp& numA, const comp& numB) {
	if (abs(numA)>MIN_NUMBER || abs(numB)>MIN_NUMBER) return 2.0*abs(numA-numB)/sqrt(norm(numA)+norm(numB));
	else return 0.0;
}

// absDiff Eigen::VectorXd
double absDiff(const vec& A, const vec& B) {
	double normSqrdA = A.squaredNorm();
	double normSqrdB = B.squaredNorm();
	vec diff = A-B;
	double normDiff = diff.norm();
	if (normSqrdA>MIN_NUMBER || normSqrdB>MIN_NUMBER ) return 2.0*normDiff/pow(normSqrdA+normSqrdB,0.5);
	else return 0.0;
}

// absDiff Eigen::VectorXcd
double absDiff(const cVec& A, const cVec& B) {
	double normSqrdA = A.squaredNorm();
	double normSqrdB = B.squaredNorm();
	cVec diff = A-B;
	double normDiff = diff.norm();
	if (normSqrdA>MIN_NUMBER || normSqrdB>MIN_NUMBER ) return 2.0*normDiff/pow(normSqrdA+normSqrdB,0.5);
	else return 0.0;
}

// absDiff Eigen::MatrixXd
double absDiff(const mat& A, const mat& B) {
	double normSqrdA = A.squaredNorm();
	double normSqrdB = B.squaredNorm();
	mat diff = A-B;
	double normDiff = diff.norm();
	if (normSqrdA>MIN_NUMBER || normSqrdB>MIN_NUMBER ) return 2.0*normDiff/pow(normSqrdA+normSqrdB,0.5);
	else return 0.0;
}

// absDiff Eigen::MatrixXcd
double absDiff(const cMat& A, const cMat& B) {
	double normSqrdA = A.squaredNorm();
	double normSqrdB = B.squaredNorm();
	cMat diff = A-B;
	double normDiff = diff.norm();
	if (normSqrdA>MIN_NUMBER || normSqrdB>MIN_NUMBER ) return 2.0*normDiff/pow(normSqrdA+normSqrdB,0.5);
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
	4. currentDateTime
	
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
	6. countLines
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
	8. explicit instantiation
		- numberToString, stringToNumber
		- smallestLoc
-------------------------------------------------------------------------------------------------------------------------*/

template string numberToString<int>(const int&);
template string numberToString<uint>(const uint&);
template string numberToString<double>(const double&);
template string numberToString<comp>(const comp&);

template int stringToNumber<int>(const string&);
template uint stringToNumber<uint>(const string&);
template double stringToNumber<double>(const string&);
template comp stringToNumber<comp>(const string&);

template uint smallestLoc<int>(const vector<int>&);
template uint smallestLoc<uint>(const vector<uint>&);
template uint smallestLoc<double>(const vector<double>&);
template uint smallestLoc<comp>(const vector<comp>&);


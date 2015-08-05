/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for Check class
 -------------------------------------------------------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include "check.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. Check class
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. Check class
		- constructor
		- add
		- good
		- checkMessage
		- back
		- tests
		- closeness
-------------------------------------------------------------------------------------------------------------------------*/

// constructor
Check::Check(const string& m, const double& c): Message(m), Closeness(c), Tests() {}

// add
void Check::add(const double& test) {
	double t = test;
	if (t<0) t = -t;
	Tests.push_back(t);
}

// good
bool Check::good() const {
	if (Tests.empty()) return false;
	return Tests.back()<Closeness;
}

// checkMessage
void Check::checkMessage() const {
	checkMessage(cout);
}

// checkMessage
void Check::checkMessage(ostream& os) const {
	if (!good()) {
		if (Tests.empty())
			os << "Check " << Message << ": no tests done yet" << endl;
		else
			os << "Check " << Message << ": test(" << Tests.back() << ") > closeness(" << Closeness << ")" << endl;
	}
}
	
// back
double Check::back() const {
	if (Tests.empty()) return 1.0;
	else return Tests.back();
}

// tests
vector<double> Check::tests() const {
	return Tests;
}

// closeness
double Check::closeness() const {
	return Closeness;
}

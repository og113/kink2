/*-------------------------------------------------------------------------------------------------------------------------
 	declarations for Check class
 -------------------------------------------------------------------------------------------------------------------------*/
 
#ifndef __CHECK_H_INCLUDED__
#define __CHECK_H_INCLUDED__

#include <iostream>
#include <vector>

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. Check class
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. Check class
-------------------------------------------------------------------------------------------------------------------------*/

class Check {
public:
	Check(const string& m, const double& c);		// constructor
	~Check() {}										// destructor
	void add(const double& t);						// add test to Tests
	bool good() const;								// carry out check
	void checkMessage() const;						// carry out check and print message if !good()
	double back() const;							// gives result of last test
	vector<double> tests() const;					// gives Tests
	double closeness() const;						// gives closeness
private:
	string 			Message;
	double 			Closeness;
	vector<double>	Tests;
};

#endif // __CHECK_H_INCLUDED__

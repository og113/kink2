/*
	analysis
		for analysis of monte carlo data, and newton-raphson data
*/

#ifndef __ANALYSIS_H_INCLUDED__
#define __ANALYSIS_H_INCLUDED__

#include <iostream>
#include <string>
#include <vector>
#include "simple.h"
#include "parameters.h"

/*-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. NewtonRaphsonDatum class
	2. NewtonRaphsonData class
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. NewtonRaphsonDatum class
-------------------------------------------------------------------------------------------------------------------------*/

class NewtonRaphsonDatum {
public:
	// constructors, destructor
	NewtonRaphsonDatum(const uint& idsize, const uint& datumsize);
	NewtonRaphsonDatum(const vector<string>&, const Parameters&, const vector<number>&);
	~NewtonRaphsonDatum();
	
	// save and load ascii
	void save(const string&) const;
	void load(const string&);
	void load(const vector<string>&);
	void load(const vector<string>&, const Parameters&, const vector<number>&);
	
	// get
	uint size() const;
	vector<string> id() const;
	Parameters parameters() const;
	vector<number> datum() const;
	vector<string> strings() const;
	
	// checks
	bool checkID(const vector<string>&) const;
	bool checkParameters(const Parameters&) const;
	
private:
	uint IDSize;
	vector<string> ID;
	Parameters P;
	uint DatumSize;
	vector<number> Datum;
	
	vector<string> stringVector() const;
};

// operator==
bool operator==(const NewtonRaphsonDatum& lhs, const NewtonRaphsonDatum& rhs);

// operator<<
ostream& operator<<(ostream&, const NewtonRaphsonDatum&);

// operator>>
istream& operator>>(istream&, NewtonRaphsonDatum&);

/*-------------------------------------------------------------------------------------------------------------------------
	2. NewtonRaphsonData class
		
-------------------------------------------------------------------------------------------------------------------------*/

class NewtonRaphsonData {
public:
	// constructors, destructor
	NewtonRaphsonData(const uint& idsize, const uint& datumsize);
	NewtonRaphsonData(const string&, const uint& idsize, const uint& datumsize);
	NewtonRaphsonData(const vector<NewtonRaphsonDatum>&);
	~NewtonRaphsonData();
	
	// save and load ascii
	void save(const string&) const;
	void saveAppend(const string&) const;
	void load(const string&);
	
	// get
	uint size() const;
	
	// add
	void add(const NewtonRaphsonDatum&);
	
	// search
	bool find(const vector<string>& id) const;
	bool find(const Parameters&) const;
	bool find(const vector<string>& id, const Parameters&) const;
	bool find(NewtonRaphsonDatum&) const;
	
private:
	uint IDSize;
	uint DatumSize;
	uint DataSize;
	vector<NewtonRaphsonDatum> DataArray;
};

#endif // __ANALYSIS_H_INCLUDED__

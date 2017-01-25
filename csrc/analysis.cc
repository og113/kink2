/*
	analysis
		for analysis of monte carlo data
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <gsl/gsl_sf_log.h>
#include "analysis.h"
#include "folder.h"
#include "print.h"
#include "simple.h"

/*-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. NewtonRaphsonDatum class
	2. NewtonRaphsonDatum classes
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. NewtonRaphsonDatum class
		
-------------------------------------------------------------------------------------------------------------------------*/

// constructor
NewtonRaphsonDatum::NewtonRaphsonDatum(const uint& idsize, const uint& datumsize): 
	IDSize(idsize), ID(), P(), DatumSize(datumsize), Datum() {
	ID.resize(IDSize);
	Datum.resize(DatumSize);
	}

// constructor
NewtonRaphsonDatum::NewtonRaphsonDatum(const vector<string>& id, const Parameters& p, const vector<number>& d): 
	ID(id), P(p), Datum(d) {
	IDSize = ID.size();
	DatumSize = Datum.size();
}

// destructor
NewtonRaphsonDatum::~NewtonRaphsonDatum() {}

// stringVector
vector<string> NewtonRaphsonDatum::stringVector() const {
	uint totalSize = IDSize + Parameters::Size + DatumSize;
	vector<string> PVector = P.valueVector();
	vector<string> v(totalSize);
	
	for (uint j=0; j<IDSize; j++)
		v[j] = ID[j]; 
	for (uint j=0; j<Parameters::Size; j++)
		v[IDSize+j] = PVector[j]; 
	for (uint j=0; j<DatumSize; j++)
		v[IDSize+Parameters::Size+j] = nts(Datum[j],16);
		
	return v;
}

// save
void NewtonRaphsonDatum::save(const string& f) const {	
	saveVectorCsvAppend(f,stringVector());
}

// load
void NewtonRaphsonDatum::load(const vector<string>& v) {
	uint totalSize = IDSize + Parameters::Size + DatumSize;
	if (v.size()!=totalSize) {
		cerr << "NewtonRaphsonDatum::load error: vector of wrong length, " << v.size() << "!=" << totalSize << endl;
		return;
	}
	for (uint j=0; j<IDSize; j++)
		ID[j] = v[j];
	vector<string> PVector(Parameters::Size);
	for (uint j=0; j<Parameters::Size; j++)
		PVector[j] = v[j+IDSize];
	P.load(PVector);
	for (uint j=0; j<DatumSize; j++)
		Datum[j] = stn<number>(v[IDSize+Parameters::Size+j]);
}

// load
void NewtonRaphsonDatum::load(const string& f) {
	vector<string> v;
	loadVectorCsvAppend(f,v);
	load(v);
}

// load
void NewtonRaphsonDatum::load(const vector<string>& id, const Parameters& p, const vector<number>& d) {
	ID = id;
	P = p;
	Datum = d;
	IDSize = ID.size();
	DatumSize = Datum.size();
}

// size
uint NewtonRaphsonDatum::size() const {
	return IDSize+Parameters::Size+DatumSize;
}

// id
vector<string> NewtonRaphsonDatum::id() const {
	return ID;
}

// parameters
Parameters NewtonRaphsonDatum::parameters() const {
	return P;
}

// datum
vector<number> NewtonRaphsonDatum::datum() const {
	return Datum;
}

// strings
vector<string> NewtonRaphsonDatum::strings() const {
	return stringVector();
}

// checkID
bool NewtonRaphsonDatum::checkID(const vector<string>& id) const {
	if (id.size()==IDSize) {
		for (uint j=0; j<IDSize; j++) {
			if (!(id[j].empty())) {
				if ((id[j]).compare(ID[j])!=0)
					return false;
			}
		}
		return true;
	}
	else
		return false;
}

// checkParameters
bool NewtonRaphsonDatum::checkParameters(const Parameters& p) const {
	return P==p;
}

// ==
bool operator==(const NewtonRaphsonDatum& lhs, const NewtonRaphsonDatum& rhs) {
	return (lhs.checkID(rhs.id()) && lhs.checkParameters(rhs.parameters()));
}

// operator<<
ostream& operator<<(ostream& os, const NewtonRaphsonDatum& d) {
	vector<string> toPrint = d.strings();
	for (uint j=0; j<toPrint.size(); j++) {
		os << toPrint[j];
		if (j<(toPrint.size()-1))
			os << ",";
	}
	return os;
}

// operator>>
istream& operator>>(istream& is, NewtonRaphsonDatum& d) {
	string line;
	getline(is,line);
	stringstream ss(line);
	vector<string> v(d.size());
	for (uint k=0; k<d.size(); k++)
		getline(ss,v[k],',');
	d.load(v);
	return is;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. NewtonRaphsonData class
		
-------------------------------------------------------------------------------------------------------------------------*/

// constructor
NewtonRaphsonData::NewtonRaphsonData(const string& f, const uint& idsize, const uint& datumsize): \
		IDSize(idsize), DatumSize(datumsize), DataSize(), DataArray() {
	load(f);
}

// constructor
NewtonRaphsonData::NewtonRaphsonData(const vector<NewtonRaphsonDatum>& v) {
	DataArray = v;
	DataSize = v.size();
}

// destructor
NewtonRaphsonData::~NewtonRaphsonData() {}

// save
void NewtonRaphsonData::save(const string& f) const {
	ofstream os;
	os.open(f.c_str());
	if (os.good()) {
		for (uint j=0; j<DataSize; j++) {
			os << DataArray[j] << endl;
		}
		os.close();
	}
	else {
		cerr << "NewtonRaphsonData::save error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// saveAppend
void NewtonRaphsonData::saveAppend(const string& f) const {
	ofstream os;
	os.open(f.c_str());
	if (os.good()) {
		for (uint j=0; j<DataSize; j++) {
			os << DataArray[j] << endl;
		}
		os.close();
	}
	else {
		cerr << "NewtonRaphsonData::save error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// load
void NewtonRaphsonData::load(const string& f) {
	Filename filename = f;
	bool fexists = filename.exists();
	if (fexists) {
		uint rows = countColumns(f,',');
		if (rows!=(IDSize+Parameters::Size+DatumSize)) {
			cerr << "NewtonRaphsonData::load error: wrong number of rows in " << f << ", "\
			 << rows << "!=" << (IDSize+Parameters::Size+DatumSize) << endl;
			return;
		}
		uint dataSize = countLines(f);
		ifstream is;
		is.open(f.c_str());
		if (is.good()) {
			NewtonRaphsonDatum temp(IDSize,DatumSize);
			for (uint j=0; j<dataSize; j++) {
				is >> temp;
				DataArray.push_back(temp);
			}
			is.close();
			DataSize = dataSize;
		}
		else {
			cerr << "NewtonRaphsonData::load error: cannot write to " << f << endl;
			is.close();
			return;
		}
	}
	else {
		return; // file doesn't exist
	}
}

// add
void NewtonRaphsonData::add(const NewtonRaphsonDatum& d) {
	if (IDSize == (d.id()).size() && DatumSize == (d.datum()).size()) {
		DataArray.push_back(d);
		DataSize++;
	}
	else {
		cerr << "NewtonRaphsonData::add error: cannot add datum as sizes wrong" << endl;
		cerr << IDSize << "," << (d.id()).size() << endl;
		cerr << DatumSize << "," << (d.datum()).size() << endl;
		return;
	}
}

// size
uint NewtonRaphsonData::size() const {
	return DataSize;
}

// find
bool NewtonRaphsonData::find(const vector<string>& id) const {
	for (uint j=0; j<DataSize; j++) {
		if ((DataArray[j]).checkID(id))
			return true;
	}
	return false;
}

// find
bool NewtonRaphsonData::find(const Parameters& p) const {
	for (uint j=0; j<DataSize; j++) {
		if ((DataArray[j]).checkParameters(p))
			return true;
	}
	return false;
}

// find
bool NewtonRaphsonData::find(const vector<string>& id, const Parameters& p) const {
	for (uint j=0; j<DataSize; j++) {
		if ((DataArray[j]).checkID(id) && (DataArray[j]).checkParameters(p))
			return true;
	}
	return false;
}

// find
bool NewtonRaphsonData::find(NewtonRaphsonDatum& d) const {
	for (uint j=0; j<DataSize; j++) {
		if (DataArray[j]==d)
			return true;
	}
	return false;
}

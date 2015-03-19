/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for the Folder class and dependencies
 -------------------------------------------------------------------------------------------------------------------------*/

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <utility> // for pair
#include <cstdlib> // for system
#include <algorithm> // for sort
#include <iterator>
#include "simple.h"
#include "folder.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. FilenameAttributes
	2. Errors
	3. Filename
	4. FilenameComparator
	5. Folder
	6. functions (reduceTo, getLastInt)
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. defintions for the FilenameAttributes class, a base class to be inherited from.
		- copy
		- copy constructor
		- operator=
		- operator<<
		- operator==
-------------------------------------------------------------------------------------------------------------------------*/

// pair
template <class T, class U>
ostream& operator<<(ostream& os, const pair<T,U>& p) {
	os << "(" << p.first << "," << p.second << ")";
	return os;
}

// copy
void FilenameAttributes::copy(const FilenameAttributes& fa) {
	Directory 	= fa.Directory;
	Timenumber 	= fa.Timenumber;
	ID			= fa.ID;
	Suffix		= fa.Suffix;
	Extras		= fa.Extras;
}

// copy constructor
FilenameAttributes::FilenameAttributes(const FilenameAttributes& fa) {
	copy(fa);
}

//operator=
FilenameAttributes& FilenameAttributes::operator=(const FilenameAttributes& rhs) {
	copy(rhs);
	return *this;
}

// operator<<
ostream& operator<<(ostream& os, const FilenameAttributes& fa) {
	os << "Directory:  " << fa.Directory << endl;
	os << "Timenumber: " << fa.Timenumber << endl;
	os << "ID:         " << fa.ID << endl;
	os << "Extras:     ";
	if ((fa.Extras).size()>0) {
		for (unsigned int l=0; l<(fa.Extras).size(); l++) {
			os << (fa.Extras[l]).first << ", " << (fa.Extras[l]).second << endl;
		}
	}
	os << endl;
	os << "Suffix:     " << fa.Suffix << endl;
	return os;
}

// operator==
bool operator==(const FilenameAttributes& lhs, const FilenameAttributes& rhs) {
	if ((lhs.Directory).compare(rhs.Directory)!=0) return false;
	if ((lhs.Timenumber).compare(rhs.Timenumber)!=0) return false;
	if ((lhs.ID).compare(rhs.ID)!=0) return false;	
	if ((lhs.Suffix).compare(rhs.Suffix)!=0) return false;
	if ((lhs.Extras).size()!=(rhs.Extras).size()) return false;
	bool ExtraOK;
	for (unsigned int n=0; n<(lhs.Extras).size(); n++) {
		ExtraOK = false;
		for (unsigned int m=0; m<(lhs.Extras).size(); m++) {
			if (((lhs.Extras[m]).first).compare(((rhs.Extras[n]).first))==0) ExtraOK = true;
		}
		if (!ExtraOK) return false;
	}
	return true;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. definitions for Filename etc errors Errors
		- FilenameError::Extras
		- FilenameComparatorError::LU
		- FolderError::System
-------------------------------------------------------------------------------------------------------------------------*/

string FilenameError::Extras::message() const {
	return "Filename error: Extras not in pairs in " + Filename;
}

string FilenameComparatorError::LU::message() const {
	return "FilenameComparator error: Lower." + Property + " = " + Lower + ", Upper." + Property + " = " + Upper;
}

string FolderError::System::message() const {
	return "Folder error: system call failure, finding files in data";
}

string FolderError::Add::message() const {
	return "Folder error: cannot add " + Filename + " as not consistent with Comparator";
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. declarations for the Filename class, publicly inherited from FilenameAttributes.
		- set
		- operator=
		- constructor(const string& filename)
		- operator string() - conversion
		- operator()
		- operator<<
		- operator<
		- operator>
-------------------------------------------------------------------------------------------------------------------------*/

// set
void Filename::set(const string& f) {
	Extras.clear();
	string temp = f;
	size_t stop;
	stop = temp.find_last_of("/");
	if (stop!=string::npos) {
		Directory = temp.substr(0,stop);
		temp = temp.substr(stop+1);
	}
	if (temp.find_first_of("0123456789")==0) {
		stop = temp.find_first_not_of("0123456789");
		Timenumber = temp.substr(0,stop);
		temp = temp.substr(stop);
	}
	if (temp.find_first_not_of("_")==0) {
	 stop = temp.find_first_of("_.");
	 ID = temp.substr(0,stop);
	 temp = temp.substr(stop);
	}
	if (temp[0]=='_') {
		temp = temp.substr(1);
		while (stop!=string::npos && temp[stop]!='.') {
			stop = temp.find("_");
			if (stop==string::npos) {
				FilenameError::Extras e(f);
				cerr << e;
			}
			StringPair sp;
			sp.first = temp.substr(0,stop);
			temp = temp.substr(stop+1);
			stop = temp.find_first_of("_.");
			sp.second = temp.substr(0,stop);
			Extras.push_back(sp);
		}
		temp = temp.substr(stop);
	}
	if (stop!=string::npos && temp[0]=='.') {
		Suffix = temp;
	}
}

// operator=
Filename& Filename::operator=(const Filename& rhs) {
	FilenameAttributes::operator=(rhs);
	return *this;
}

// operator=
Filename& Filename::operator=(const string& rhs) {
	set(rhs);
	return *this;
}


// constructor(const string& filename)
Filename::Filename(const string& f): FilenameAttributes() {
	set(f);
}

// operator string() - conversion
Filename::operator string() const {
	string filename = Directory + "/" + Timenumber + ID;
	for (unsigned int l=0; l<Extras.size(); l++) {
		filename += "_" + Extras[l].first + "_" + Extras[l].second;
	}
	filename += Suffix;
	return filename;
}

// operator()
string Filename::operator()() const {
	string filename = Directory + "/" + Timenumber + ID;
	for (unsigned int l=0; l<Extras.size(); l++) {
		filename += "_" + Extras[l].first + "_" + Extras[l].second;
	}
	filename += Suffix;
	return filename;
}

// operator<<
ostream& operator<<(ostream& os, const Filename& f) {
	os << f();
	return os;
}

// operator >>
istream& operator>>(istream& is, Filename& f) {
	string filename;
	is >> filename;
	f = filename;
	return is;
}

// operator<
bool operator<(const Filename& lhs, const Filename& rhs) {
	return lhs()<rhs();
}

// operator>
bool operator>(const Filename& lhs, const Filename& rhs) {
	return lhs()>rhs();
}

/*-------------------------------------------------------------------------------------------------------------------------
	4. defintions for the FilenameComparator class, which is used by the Folder class. FilenameComparator sees the ugly details of Filename.
		- copy
		- copy constructor
		- check
		- constructor(lower,upper)
		- operator=
		- set
		- setLower
		- setUpper
		- operator(Filename)
		- <<
-------------------------------------------------------------------------------------------------------------------------*/

// copy
void FilenameComparator::copy(const FilenameComparator& fc) {
	Lower = fc.Lower;
	Upper = fc.Upper;
}

// copy constructor
FilenameComparator::FilenameComparator(const FilenameComparator& fc) {
	copy(fc);
}

// check
bool FilenameComparator::check(const FilenameAttributes& low, const FilenameAttributes& u) const {
	if ((low.Directory).compare(u.Directory)!=0) {
		FilenameComparatorError::LU e("Directory",low.Directory,u.Directory);
		cerr << e;
		return false;
		}
	if ((low.ID).compare(u.ID)!=0) {
		FilenameComparatorError::LU e("ID",low.ID,u.ID);
		cerr << e;
		return false;
	}
	if ((low.Suffix).compare(u.Suffix)!=0) {
		FilenameComparatorError::LU e("Suffix",low.Suffix,u.Suffix);
		cerr << e;
		return false;
	}
	if ((low.Extras).size()!=(u.Extras).size()) {
		FilenameComparatorError::LU e("Extras.size()",numberToString<unsigned int>((low.Extras).size())\
								,numberToString<unsigned int>((u.Extras).size()));
		cerr << e;
		return false;
	}
	bool ExtraOK;
	for (unsigned int n=0; n<(low.Extras).size(); n++) {
		ExtraOK = false;
		for (unsigned int m=0; m<(low.Extras).size(); m++) {
			if (((low.Extras[m]).first).compare(((u.Extras[n]).first))==0) ExtraOK = true;
		}
		if (!ExtraOK) {
			FilenameComparatorError::LU e("Extras",(low.Extras[n]).first,u.Extras[n].first);
			cerr << e;
			return false;
		}
	}
	return true;
}

// constructor(lower,upper)
FilenameComparator::FilenameComparator(const FilenameAttributes& l, const FilenameAttributes& u): Lower(l), Upper(u) {
	check(l,u);
}

// operator=
FilenameComparator& FilenameComparator::operator=(const FilenameComparator& rhs) {
	copy(rhs);
	return *this;
}

// set
void FilenameComparator::set(const FilenameAttributes& l, const FilenameAttributes& u) {
	if (check(l,u)) {
		Lower = l;
		Upper = u;
	}
}

// setLower
void FilenameComparator::setLower(const FilenameAttributes& l) {
	if (check(l,Upper)) {
		Lower = l;
	}
}

// setUpper
void FilenameComparator::setUpper(const FilenameAttributes& u) {
	if (check(Lower,u)) {
		Upper = u;
	}
}

// operator(Filename)
bool FilenameComparator::operator()(const Filename& f) const{
	if (!(Lower.Directory).empty()) {
		if ((f.Directory).compare(Lower.Directory)!=0) return false;
	}
	if (!(Lower.Timenumber).empty()) {
		if (f.Timenumber<Lower.Timenumber || f.Timenumber > Upper.Timenumber) return false;
	}
	if (!(Lower.ID).empty()) {
		if ((f.ID).compare(Lower.ID)!=0) return false;
	}
	if (!(Lower.Suffix).empty()) {
		if ((f.Suffix).compare(Lower.Suffix)!=0) return false;
	}
	size_t NumExtras = (Lower.Extras).size();
	if (NumExtras>0) {
		if ((f.Extras).size()!=NumExtras) return false;
		bool ExtraOK;
		for (unsigned int n=0; n<NumExtras; n++) {
			ExtraOK = false;
			for (unsigned int m=0; m<NumExtras; m++) {
				if (((Lower.Extras[m]).first).compare(((f.Extras[n]).first))==0) {
					if (((f.Extras[n]).second)>=((Lower.Extras[m]).second) && ((f.Extras[n]).second)<=((Upper.Extras[m]).second))
						ExtraOK = true;
				}
			}
			if (!ExtraOK) return false;
		}
	}
	return true;
}

// operator<<
ostream& operator<<(ostream& os, const FilenameComparator& fc){
	os << "Lower: " << endl << fc.Lower << endl << "Upper: " << fc.Upper << endl;
	return os;
}

/*-------------------------------------------------------------------------------------------------------------------------
	5. definitions for the Folder class.
		- isPresent(Filename)
		- refresh
		- sort
		- order
		- update
		- begin
		- end
		- add
		- erase
		- copy
		- copy constructor
		- constructor(FilenameComparator)
		- operator=
		- set
		- size
		- operator[]
		- <<
-------------------------------------------------------------------------------------------------------------------------*/

// isPresent(Filename)
bool Folder::isPresent(const Filename& f) {
	for (unsigned int l=0; l<Filenames.size(); l++) {
		if (f==Filenames[l]) return true;
	}
	return false;
}

// begin
FolderIterator Folder::begin() {
	return Filenames.begin();
}

// begin
ConstFolderIterator Folder::begin() const{
	return Filenames.begin();
}

// end
FolderIterator Folder::end() {
	return Filenames.end();
}

// end
ConstFolderIterator Folder::end() const{
	return Filenames.end();
}

// sort
void Folder::sort() {
	std::sort(Filenames.begin(),Filenames.end());
}

// order
void Folder::order() {
	Folder::sort();
}

// refresh
void Folder::refresh() {
	int systemCall = system("find data/* -type f > data/dataFiles");
	if (systemCall==-1) {
		FolderError::System e;
		cerr << e;
	}
	ifstream is;
	Filename f;
    is.open ("data/dataFiles");
	while ( !is.eof() ){
		is >> f;
		if (!isPresent(f) && Comparator(f)) Filenames.push_back(f);
	}
    is.close();
    sort();
}

// update
void Folder::update() {
	refresh();
}

// add
void Folder::add(const Filename& f) {
	if(Comparator(f)) {
		Filenames.push_back(f);
		sort();
	}
	else {
		FolderError::Add e(f());
		cerr << e;
	}
}

// erase
void Folder::erase(FolderIterator it) {
	Filenames.erase(it);
}

// copy
void Folder::copy(const Folder& f) {
	Comparator = f.Comparator;
	Filenames = f.Filenames;
}

// copy constructor
Folder::Folder(const Folder& f): Comparator(), Filenames() {
	copy(f);
}

// constructor(FilenameComparator)
Folder::Folder(const FilenameComparator& fc): Comparator(fc), Filenames() {
	refresh();
}

// constructor(FilenameAttributes)
Folder::Folder(const FilenameAttributes& l): Comparator(l), Filenames() {
	refresh();
}

// constructor(FilenameAttributes, FilenameAttributes)
Folder::Folder(const FilenameAttributes& l, const FilenameAttributes& u): Comparator(l,u), Filenames() {
	refresh();
}

// operator=
Folder& Folder::operator=(const Folder& f) {
	copy(f);
	return *this;
}

// set
void Folder::set(const FilenameComparator& fc) {
	Comparator = fc;
	refresh();
}

// size
unsigned int Folder::size() const{
	return Filenames.size();
}

// operator[]
Filename Folder::operator[](const int& index) const {
	if (index<0 || index>(size()-1)) {
		IndexError::OutOfBounds e(index,size(),0);
		cerr << e;
	}
	return Filenames[index];
}

// operator<<
ostream& operator<<(ostream& os, const Folder& f) {
	for (unsigned int l=0; l<f.size(); l++)
		os << f[l] << endl;
	return os;
}

/*-------------------------------------------------------------------------------------------------------------------------
	6. functions acting on Filenames and Folders
		- removeUnshared
		- getLastInt (used in getting ints from zmx and zmt)
-------------------------------------------------------------------------------------------------------------------------*/

// removeUnshared
void removeUnshared(Folder& f1,Folder& f2) {
	if (f1.size()==0) 		f2 = f1;
	else if (f2.size()==0) 	f1 = f2;
	else{
		for (unsigned int j=0;j<f1.size();j++) {
			if(find(f2.begin(), f2.end(), f1[j]) == f2.end()) {
				f1.erase(f1.begin()+j);
			}
		}
		for (unsigned int j=0;j<f2.size();j++) {
			if(find(f1.begin(), f1.end(), f2[j]) == f1.end()) {
				f2.erase(f2.begin()+j);
			}
		}
	}
	f1.order();
	f2.order();
}

// getLastInt - returns last integer in string
uint getLastInt(const string& str) {
	size_t first_index;
	size_t last_index;
	if (str.find_last_not_of("0123456789")!=string::npos) first_index = str.find_last_not_of("0123456789");
	else {
		cerr << "getLastInt error, no non-numeric characters found in " << str << endl;
		return -1;
	}
	string temp = str.substr(first_index+1);
	if (isdigit(temp[0])==0) { // zero for false
		cerr << "getLastInt error, first character of " << temp << " in " << str << " not numeric" << endl;
		return -1;
	}
	if (temp.find_last_of("0123456789")!=string::npos) last_index = temp.find_last_of("0123456789");
	else {
		cerr << "getLastInt error, last numeric character not found in " << temp << " in " << str << endl;
		return -1;
	}
	temp = temp.substr(0,last_index+1);
	return stringToNumber<uint>(temp);
}

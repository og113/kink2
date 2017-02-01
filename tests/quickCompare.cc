/*
	quick comparison of some output numbers from main3 and main3old
*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "simple.h"

using namespace std;

int main() {
cout << "quick comparison of main3 and main3old outputs: " << endl;

string file = "temp/main3_out.txt";
string fileOld = "temp/main3old_out.txt";
string line, lineOld, object, dross;

double d, dOld;

ifstream is, isOld;
is.open(file.c_str());
isOld.open(fileOld.c_str());

cout << left << setprecision(20) << endl;
while (!is.eof() && !isOld.eof()) {
	getline(is,line);
	getline(isOld,lineOld);
	istringstream iss(line);
	istringstream issOld(lineOld);
	
	iss >> object >> dross >> d;
	issOld >> object >> dross >> dOld;
	
	if (dross.compare("=")==0) {
		cout << setw(30) << object << setw(30) << d << setw(30) << dOld;
		cout << setw(30) << d-dOld;
		cout << setw(30) << (d-dOld)/(0.5*(abs(d)+abs(dOld)))/MIN_NUMBER << endl;
	}
}

is.close();
isOld.close();

return 0;
}

/*
	main for program to test simple.cc
*/

#include <iostream>
#include "simple.h"

using namespace std;

int main() {
cout << "test simple: " << endl;

cout << "test of mod: " << endl;
double a=2.5*pi, b=2.0*pi, c=1.5*pi, d=pi, e=0.5*pi, f=0.0, g=-0.5*pi, h=-pi, i=-1.5*pi, j=-2.0*pi, k=-2.5*pi;
cout << mod(a,-pi,pi) << endl;
cout << mod(b,-pi,pi) << endl;
cout << mod(c,-pi,pi) << endl;
cout << mod(d,-pi,pi) << endl;
cout << mod(e,-pi,pi) << endl;
cout << mod(f,-pi,pi) << endl;
cout << mod(g,-pi,pi) << endl;
cout << mod(h,-pi,pi) << endl;
cout << mod(i,-pi,pi) << endl;
cout << mod(j,-pi,pi) << endl;
cout << mod(k,-pi,pi) << endl;

cout << "test of numberToString and stringToNumber:" << endl;
string strl = "12345678";
lint numl = stringToNumber<lint>(strl);
cout << "numl = " << numl << endl;
string strll = "123456789101";
unsigned long long numll = stringToNumber<unsigned long long>(strll);
cout << "numll = " << numll << endl;


return 0;
}

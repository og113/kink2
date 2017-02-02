/*
	test complex
*/

#include <cmath>
#include <iostream>
#include <complex>

using namespace std;

int main() {

complex<double> z(2.75859e-313,2.75859e-313);

cout << z << endl;
cout << z*z << endl;
cout << pow(z,2) << endl;
cout << pow(z,2.0) << endl;

return 0;
}

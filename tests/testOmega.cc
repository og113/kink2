/*
	main for program to test omega.cc
*/

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "simple.h"
#include "omega.h"

using namespace std;

int main() {
cout << "test omega: " << endl;

// loading parameters
Parameters ps;
ps.load("inputsM");
ps.print();

// analytic
vec freqsA(ps.N), freqs_expA(ps.N);
mat modesA(ps.N,ps.N);
mat omega_m1A(ps.N,ps.N), omega_0A(ps.N,ps.N), omega_1A(ps.N,ps.N), omega_2A(ps.N,ps.N);
analyticModes(modesA,freqsA,freqs_expA,ps);
omegasFn(true,modesA,freqsA,omega_m1A,omega_0A,omega_1A,omega_2A,ps);

// numerical
vec freqsN(ps.N), freqs_expN(ps.N);
mat modesN(ps.N,ps.N);
mat omega_m1N(ps.N,ps.N), omega_0N(ps.N,ps.N), omega_1N(ps.N,ps.N), omega_2N(ps.N,ps.N);
numericalModes(modesN,freqsN,freqs_expN,ps);
omegasFn(false,modesN,freqsN,omega_m1N,omega_0N,omega_1N,omega_2N,ps);

// difference
double freqsD = absDiff(freqsA,freqsN);
double freqs_expD = absDiff(freqs_expA,freqs_expN);
double modesD = absDiff(modesA,modesN);
double omega_m1D = absDiff(omega_m1A,omega_m1N);
double omega_0D = absDiff(omega_0A,omega_0N);
double omega_1D = absDiff(omega_1A,omega_1N);
double omega_2D = absDiff(omega_2A,omega_2N);

cout << "freqsA = " << endl << freqsA << endl;
cout << "freqsN = " << endl << freqsN << endl;

// printing resutls
cout << "difference in freqs = " << freqsD << endl;
cout << "difference in freqs_exp = " << freqs_expD << endl;
cout << "difference in modes = " << modesD << endl;
cout << "difference in omega_m1 = " << omega_m1D << endl;
cout << "difference in omega_0 = " << omega_0D << endl;
cout << "difference in omega_1 = " << omega_1D << endl;
cout << "difference in omega2D = " << omega_2D << endl;
cout << endl;

return 0;
}

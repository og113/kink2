/*
	main for program to test omega.cc
*/

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "simple.h"
#include "omega.h"

using namespace std;

// for sorting vectors
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

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

// numerical
vec freqsN(ps.N), freqs_expN(ps.N);
mat modesN(ps.N,ps.N);
mat omega_m1N(ps.N,ps.N), omega_0N(ps.N,ps.N), omega_1N(ps.N,ps.N), omega_2N(ps.N,ps.N);
numericalModes(modesN,freqsN,freqs_expN,ps);

// sorting vectors
vector<double> freqsNv(ps.N), freqsAv(ps.N);
for (size_t j=0; j<ps.N; j++) {
	freqsNv[j] = freqsN(j);
	freqsAv[j] = freqsA(j);
}
vector<size_t> orderN = sort_indexes(freqsNv);
vector<size_t> orderA = sort_indexes(freqsAv);
vec freqs_expNc = freqs_expN, freqs_expAc = freqs_expA;
mat modesNc = modesN, modesAc = modesA;
for (size_t j=0; j<ps.N; j++) {
	freqsN(j) = freqsNv[orderN[j]];
	freqsA(j) = freqsAv[orderA[j]];
	freqs_expN(j) = freqs_expNc(orderN[j]);
	freqs_expA(j) = freqs_expAc(orderA[j]);
	modesN.col(j) = modesNc.col(orderN[j]);
		if ((modesN.col(j))(1)<0) modesN.col(j) *= -1.0;
	modesA.col(j) = modesAc.col(orderA[j]);
		if ((modesA.col(j))(1)<0) modesA.col(j) *= -1.0;
}

// omegas
omegasFn(false,modesN,freqsN,omega_m1N,omega_0N,omega_1N,omega_2N,ps);
omegasFn(true,modesA,freqsA,omega_m1A,omega_0A,omega_1A,omega_2A,ps);

// difference
double freqsD = absDiff(freqsA,freqsN);
double freqs_expD = absDiff(freqs_expA,freqs_expN);
double modesD = absDiff((Eigen::MatrixXd)modesA.block(2,2,ps.N-2,ps.N-2),(Eigen::MatrixXd)modesN.block(2,2,ps.N-2,ps.N-2));
double omega_m1D = absDiff((Eigen::MatrixXd)omega_m1A.block(2,2,ps.N-2,ps.N-2),(Eigen::MatrixXd)omega_m1N.block(2,2,ps.N-2,ps.N-2));
double omega_0D = absDiff((Eigen::MatrixXd)omega_0A.block(2,2,ps.N-2,ps.N-2),(Eigen::MatrixXd)omega_0N.block(2,2,ps.N-2,ps.N-2));
double omega_1D = absDiff((Eigen::MatrixXd)omega_1A.block(2,2,ps.N-2,ps.N-2),(Eigen::MatrixXd)omega_1N.block(2,2,ps.N-2,ps.N-2));
double omega_2D = absDiff((Eigen::MatrixXd)omega_2A.block(2,2,ps.N-2,ps.N-2),(Eigen::MatrixXd)omega_2N.block(2,2,ps.N-2,ps.N-2));

// printing resutls
cout << "difference in freqs = " << freqsD << endl;
cout << "difference in freqs_exp = " << freqs_expD << endl;
cout << "difference in modes = " << modesD << endl;
cout << "difference in omega_m1 = " << omega_m1D << endl;
cout << "difference in omega_0 = " << omega_0D << endl;
cout << "difference in omega_1 = " << omega_1D << endl;
cout << "difference in omega_2 = " << omega_2D << endl;
cout << "n.b. differences in matrices exclude edges" << endl;
cout << endl;

// some quick extra checks
mat modesAred = modesA.block(2,2,ps.N-2,ps.N-2);
mat modesNred = modesN.block(2,2,ps.N-2,ps.N-2);
mat modesDmat = modesAred-modesNred;
cout << "some extra quick checks: " << endl;
cout << (modesDmat.col(0)).norm() << endl;
cout << (modesDmat.col(1)).norm() << endl;
cout << (modesDmat.col(2)).norm() << endl;
cout << (modesDmat.col(3)).norm() << endl;
cout << (modesDmat.col(4)).norm() << endl;

return 0;
}

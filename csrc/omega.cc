/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for the functions to calculate omega
 -------------------------------------------------------------------------------------------------------------------------*/
 
 #include <iostream>
 #include <cmath>
 #include <complex>
#include <Eigen/Dense>
#include "parameters.h"
#include "lattice.h" // for DxFn

using namespace std;

typedef unsigned int uint;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. h matrix
	2. analytic modes
	3. numerical modes
	4. constructing omegas
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. h matrix
		- hFn
-------------------------------------------------------------------------------------------------------------------------*/

// hFn
static mat	hFn(const Parameters& p) {
	mat xh(p.N,p.N);	xh = Eigen::MatrixXd::Zero(p.N,p.N);
	double diag = p.mass2 + 2.0/pow(p.a,2.0);
	double offDiag1 = -1.0/pow(p.a,2.0);
	double offDiag2 = (p.pot==3? -pow(2.0,0.5)/pow(p.a,2.0): -1.0/pow(p.a,2.0) );
	for (unsigned int l=0; l<p.N; l++) {
		if (l==0) {
			if (p.pot==3) 	xh(l,l) = 1.0; 					// taking into account boundary conditions
			else {
							xh(l,l) = diag;
							xh(l,l+1) = offDiag2;	
			}		
		}
		else if (l==(p.N-1)) {
			if (p.pot==3) 	xh(l,l) = 1.0; 					// taking into account boundary conditions
			else {
							xh(l,l) = diag;
							xh(l,l-1) = offDiag2;
			}
		}
		else {
			xh(l,l) = diag;
			xh(l,l+1) = ((l+1)==(p.N-1)? 	offDiag2: offDiag1);
			xh(l,l-1) = ((l-1)==0?			offDiag2: offDiag1);
		}
	}
	return xh;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. analytic modes
-------------------------------------------------------------------------------------------------------------------------*/

// analytic modes
void analyticModes(mat& modes, vec& freqs, vec& freqs_exp, const Parameters& p) {
	double normalisation = sqrt(2.0/(p.N-1.0));
	for (unsigned int l=0; l<p.N; l++) {
		freqs(l) = sqrt(1.0+pow(2.0*sin(PI*l/(p.N-1.0)/2.0)/p.a,2.0));
		freqs_exp(l) = (2.0/p.b)*asin(p.b*freqs(l)/2.0);
		for (unsigned int m=0; m<p.N; m++) {
			if (p.pot==3) 	modes(l,m) = normalisation*sin(PI*l*m/(p.N-1.0));
			else			modes(l,m) = normalisation*cos(PI*l*m/(p.N-1.0));
		}
	}
	freqs(p.N-1) = 1.0;
	freqs_exp(p.N-1) = (2.0/p.b)*asin(p.b*freqs(p.N-1)/2.0);
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. numerical modes
-------------------------------------------------------------------------------------------------------------------------*/

// numerical modes
void numericalModes(mat& modes, vec& freqs, vec& freqs_exp, const Parameters& p) {
	mat h = hFn(p);
	Eigen::EigenSolver<mat> eigensolver(h);
	cVec cFreqs2(p.N);
	cMat cModes(p.N,p.N);
	if (eigensolver.info() != Eigen::Success) {
		cerr << "h eigensolver failed" << endl;
		cerr << "N = " << p.N << ", a = " << p.a << ", mass2 = " << p.mass2 << endl;
	}
	else {
		cFreqs2 = eigensolver.eigenvalues();
		cModes = eigensolver.eigenvectors(); //automatically normalised to have unit norm
	}
	for (unsigned int j=0; j<p.N; j++) {
		freqs(j) = sqrt(real(cFreqs2(j)));
		freqs_exp(j) = (2.0/p.b)*asin(p.b*freqs(j)/2.0);
		for (unsigned int k=0; k<p.N; k++) {
			modes(j,k) = real(cModes(j,k));
			if (abs(imag(cModes(j,k)))>MIN_NUMBER)
				cerr << "omega error: imaginary part of modes(" << j << "," << k << ") = " << imag(cModes(j,k)) << endl;
		}
	}
}
/*-------------------------------------------------------------------------------------------------------------------------
	4. constructing omegas
-------------------------------------------------------------------------------------------------------------------------*/

// omegasFn
void omegasFn(const bool& analytic, const mat& modes, const vec& freqs, mat& omega_m1, mat& omega_0, mat& omega_1, mat& omega_2, const Parameters& p) {
	omega_m1 = Eigen::MatrixXd::Zero(p.N,p.N);
	omega_0 = Eigen::MatrixXd::Zero(p.N,p.N);
	omega_1 = Eigen::MatrixXd::Zero(p.N,p.N);
	omega_2 = Eigen::MatrixXd::Zero(p.N,p.N);
	double djdk;
	for (unsigned int j=0; j<p.N; j++) {
		for (unsigned int k=0; k<p.N; k++) {
			for (unsigned int l=0; l<p.N; l++) {
				if (analytic) 	 djdk = p.a;
				else 			 djdk = sqrt(DxFn(j,p)*DxFn(k,p));
				omega_m1(j,k) += djdk*pow(freqs(l),-1.0)*modes(j,l)*modes(k,l);
				omega_0(j,k)  += djdk*modes(j,l)*modes(k,l);
				omega_1(j,k)  += djdk*freqs(l)*modes(j,l)*modes(k,l);
				omega_2(j,k)  += djdk*pow(freqs(l),2.0)*modes(j,l)*modes(k,l);
			}
		}
	}
}

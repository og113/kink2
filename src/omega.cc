/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for the functions to calculate omega
 -------------------------------------------------------------------------------------------------------------------------*/
 
#include <Eigen/Dense>
#include "parameters.h"
#include "lattice.h" // for DxFn

using namespace std;

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
mat	hFn(const Parameters& p) {
	mat xh(p.N,p.N);	xh = Eigen::MatrixXd::Zero(p.N,p.N);
	double diag = p.mass2 + 2.0/pow(p.a,2.0);
	double offDiag1 = -1.0/pow(p.a,2.0);
	double offDiag2 = (p.pot[0]=='3'? -pow(2.0,0.5)/pow(p.a,2.0): -1.0/pow(p.a,2.0) );
	for (unsigned int l=0; l<p.N; l++)
		{
		if (l==0)
			{
			xh(l,l) = 1.0; 					// taking into account boundary conditions
			//xh(l,l) = diag;
			//xh(l,l+1) = offDiag2;			
			}
		else if (l==(p.N-1))
			{
			xh(l,l) = 1.0; 					// taking into account boundary conditions
			//xh(l,l) = diag;
			//xh(l,l-1) = offDiag2;
			}
		else
			{
			xh(l,l) = diag;
			xh(l,l+1) = ((l+1)==(p.N-1)? 	offDiag2: offDiag1);
			xh(l,l-1) = ((l-1)==0?			offDiag2: offDiag1);
			}
	return xh;
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. analytic modes
-------------------------------------------------------------------------------------------------------------------------*/

// analytic modes
void analyticModes(const mat& h, mat& modes, vec& freqs, const Parameters& p) {
	double normalisation = sqrt(2.0/(p.N-1.0));
	for (unsigned int l=0; l<p.N; l++) {
		freqs(l) = 1.0+pow(2.0*sin(pi*l/(p.N-1.0)/2.0)/p.a,2.0);
		for (unsigned int m=0; m<p.N; m++) {
			if (pot[0]=='3') modes(l,m) = normalisation*sin(pi*l*m/(p.N-1.0));
			else			 modes(l,m) = normalisation*cos(pi*l*m/(p.N-1.0));
		}
	}
	freqs(p.N-1) = 1.0;		
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. numerical modes
-------------------------------------------------------------------------------------------------------------------------*/

// numerical modes
void numericalModes(const mat& h, mat& modes, vec& freqs, const Parameters& p) {
	Eigen::EigenSolver<mat> eigensolver(h);
	cVec cFreqs(p.N);
	cMat cModes(p.N,p.N);
	if (eigensolver.info() != Eigen::Success) {
		cerr << "h eigensolver failed" << endl;
		cerr << "N = " << p.N << ", a = " << p.a << ", mass2 = " << p.mass2 << endl;
	}
	else {
		cFreqs = eigensolver.eigenvalues();
		cModes = eigensolver.eigenvectors(); //automatically normalised to have unit norm
	}
	for (unsigned int j=0; j<p.N; j++) {
		freqs(j) = real(cFreqs(j));
		for (unsigned int k=0; k<p.N; k++) {
			modes(j,k) = real(cModes(j,k));
		}
	}
}
/*-------------------------------------------------------------------------------------------------------------------------
	4. constructing omegas
-------------------------------------------------------------------------------------------------------------------------*/

// omegasFn
void omegasFn(const mat& modes, const mat& freqs, mat& omega_m1, mat& omega_0, mat& omega_1, mat& omega_2, const Parameters& p) {
	double djdk;
	for (unsigned int j=0; j<p.N; j++) {
		w_n_exp(j) = (2.0/p.b)*asin(p.b*sqrt(real(freqs(j)))/2.0);
		for (unsigned int k=0; k<p.N; k++) {
			for (unsigned int l=0; l<p.N; l++) {
				if (approxOmega) djdk = p.a;
				else 			 djdk = sqrt(DxFn(j,p)*DxFn(k,p));
				omega_m1(j,k) += djdk*pow(real(freqs(l)),-0.5)*real(modes(j,l))*real(modes(k,l));
				omega_0(j,k)  += djdk*real(modes(j,l))*real(modes(k,l));
				omega_1(j,k)  += djdk*pow(real(freqs(l)),0.5)*real(modes(j,l))*real(modes(k,l));
				omega_2(j,k)  += djdk*real(freqs(l))*real(modes(j,l))*real(modes(k,l));
			}
		}
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	definitions of functions for newton-raphson quantities
-------------------------------------------------------------------------------------------------------------------------*/

#include <cmath>
#include <complex>
#include "lattice.h"
#include "nr.h"
#include "simple.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - 
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1 - accessory functions
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	2 - nr scalar functions
-------------------------------------------------------------------------------------------------------------------------*/

// Kinetic_nr
void Kinetic_nr (const uint& j, const uint& dir, const vec& p, const Parameters& pr, const cVec& f, comp& result) {
	
	long pj = neigh(j,dir,1,pr); // assuming one spatial direction (the 1-direction)
	if (pj!=-1) {
		comp cp_j(p(2*j), p(2*j+1));
		comp cp_pj(p(2*pj), p(2*pj+1));
		result += f(j)*pow(cp_pj-cp_j,2.0)/2.0;
	}
	
}

// Potential_nr
void Potential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& v, const cVec& f, comp& result) {
	comp cp_j(p(2*j), p(2*j+1));
	result += f(j)*v(cp_j);
}

/*-------------------------------------------------------------------------------------------------------------------------
	3 - nr vector functions
-------------------------------------------------------------------------------------------------------------------------*/

// mdKinetic_nr
void mdKinetic_nr (const uint& j, const uint& dir, const vec& p, const Parameters& pr, const cVec& f, vec& mds) {
	long pj = neigh(j,dir,1,pr);
	long mj = neigh(j,dir,-1,pr);
	if (pj!=-1 && mj!=-1) {
		mds(2*j) -=  real(f(mj))*(-p(2*mj) + p(2*j)) + real(f(j))*(p(2*j) - p(2*pj)) \
					+ imag(f(mj))*(p(2*mj + 1) - p(2*j + 1)) + imag(f(j))*(-p(2*j + 1) + p(2*pj + 1));
		mds(2*j+1) -=  imag(f(mj))*(-p(2*mj) + p(2*j)) + imag(f(j))*(p(2*j) - p(2*pj)) \
					+ real(f(mj))*(-p(2*mj + 1) + p(2*j + 1)) + real(f(j))*(p(2*j + 1) - p(2*pj + 1));
	}
	else if (pj==-1) {
		mds(2*j) -=  real(f(mj))*(-p(2*mj) + p(2*j)) + imag(f(mj))*(p(2*mj + 1) - p(2*j + 1));
		mds(2*j+1) -=  imag(f(mj))*(-p(2*mj) + p(2*j)) + real(f(mj))*(-p(2*mj + 1) + p(2*j + 1));
	}
	else if (mj==-1) {
		mds(2*j) -= real(f(j))*(p(2*j) - p(2*pj)) + imag(f(j))*(-p(2*j + 1) + p(2*pj + 1));
		mds(2*j+1) -= imag(f(j))*(p(2*j) - p(2*pj)) + real(f(j))*(p(2*j + 1) - p(2*pj + 1));
	}
}

// mdPotential_nr
void mdPotential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& v, const cVec& f, vec& mds) {

}

/*-------------------------------------------------------------------------------------------------------------------------
	4 - nr matrix functions
-------------------------------------------------------------------------------------------------------------------------*/

// ddKinetic_nr
void ddKinetic_nr (const uint& j, const uint& dir, const uint& k, const vec& p, const Parameters& pr, const cVec& f, spMat& dds) {
	long mj = neigh(j,dir,-1,pr);
	long mk = neigh(k,dir,-1,pr);
	long pj = neigh(j,dir,1,pr);
	long pk = neigh(k,dir,1,pr);
	
	if (j==k) {
		if (mj!=-1) {
			dds.coeffRef(2*j,2*k) += + real(f(mj)) + real(f(j)); 
			dds.coeffRef(2*j,2*k+1) += - imag(f(mj)) - imag(f(j));
			dds.coeffRef(2*j+1,2*k) += + imag(f(mj)) + imag(f(j));
			dds.coeffRef(2*j+1,2*k+1) += + real(f(mj)) + real(f(j));
		}
		else {
			dds.coeffRef(2*j,2*k) +=  real(f(j)); 
			dds.coeffRef(2*j,2*k+1) += - imag(f(j));
			dds.coeffRef(2*j+1,2*k) +=  imag(f(j));
			dds.coeffRef(2*j+1,2*k+1) +=  real(f(j));
		}
		
	}

	if (j==pk && mj==k) { // not sure this double check is needed
		dds.coeffRef(2*j,2*pk) += - real(f(mj));
		dds.coeffRef(2*j,2*pk+1) += imag(f(mj));
		dds.coeffRef(2*j+1,2*pk) += - imag(f(mj));
		dds.coeffRef(2*j+1,2*pk+1) += - real(f(mj));
	}
	if (j==mk && pj==k) {
		dds.coeffRef(2*j,2*mk) += - real(f(j));
		dds.coeffRef(2*j,2*mk+1) += + imag(f(j));
		dds.coeffRef(2*j+1,2*mk) += - imag(f(j));
		dds.coeffRef(2*j+1,2*mk+1) += - real(f(j));
	}
}

// ddPotential_nr
void ddPotential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& v, const cVec& f, spMat& dds) {

}

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

// KineticS_nr
void KineticS_nr (const uint& j, const vec& p, const Parameters& pr, const cVec& f, comp& result) {
	
	long pxj = neigh(j,1,1,pr); // assuming one spatial direction (the 1-direction)
	if (pxj!=-1) {
		comp cp_j(p(2*j), p(2*j+1));
		comp cp_pxj(p(2*pxj), p(2*pxj+1));
		result += f(j)*pow(cp_pxj-cp_j,2.0)/2.0;
	}
	
}

// KineticT_nr
void KineticT_nr (const uint& j, const vec& p, const Parameters& pr, const cVec& f, comp& result) {
	
	long ptj = neigh(j,0,1,pr); // taking time to be the 0 direction
	if (ptj!=-1) {
		comp cp_j(p(2*j), p(2*j+1));
		comp cp_ptj(p(2*ptj), p(2*ptj+1));
		result += f(j)*pow(cp_ptj - cp_j,2.0)/2.0;
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

// mdKineticS_nr
void mdKineticS_nr (const uint& j, const vec& p, const Parameters& pr, const cVec& f, vec& mds) {
	long pxj = neigh(j,1,1,pr);
	long mxj = neigh(j,1,-1,pr);
	if (pxj!=-1 && mxj!=-1) {
		mds(2*j) -=  real(f(mxj))*(-p(2*mxj) + p(2*j)) + real(f(j))*(p(2*j) - p(2*pxj)) \
					+ imag(f(mxj))*(p(2*mxj + 1) - p(2*j + 1)) + imag(f(j))*(-p(2*j + 1) + p(2*pxj + 1));
		mds(2*j+1) -=  imag(f(mxj))*(-p(2*mxj) + p(2*j)) + imag(f(j))*(p(2*j) - p(2*pxj)) \
					+ real(f(mxj))*(-p(2*mxj + 1) + p(2*j + 1)) + real(f(j))*(p(2*j + 1) - p(2*pxj + 1));
	}
	else if (pxj==-1) {
		mds(2*j) -=  real(f(mxj))*(-p(2*mxj) + p(2*j)) + imag(f(mxj))*(p(2*mxj + 1) - p(2*j + 1));
		mds(2*j+1) -=  imag(f(mxj))*(-p(2*mxj) + p(2*j)) + real(f(mxj))*(-p(2*mxj + 1) + p(2*j + 1));
	}
	else if (mxj==-1) {
		mds(2*j) -= real(f(j))*(p(2*j) - p(2*pxj)) + imag(f(j))*(-p(2*j + 1) + p(2*pxj + 1));
		mds(2*j+1) -= imag(f(j))*(p(2*j) - p(2*pxj)) + real(f(j))*(p(2*j + 1) - p(2*pxj + 1));
	}
}


// mdKineticT_nr
void mdKineticT_nr (const uint& j, const vec& p, const Parameters& pr, const cVec& f, vec& mds) {
	long ptj = neigh(j,0,1,pr);
	long mtj = neigh(j,0,-1,pr);
	if (ptj!=-1 && mtj!=-1) {
		mds(2*j) -=  real(f(mtj))*(-p(2*mtj) + p(2*j)) + real(f(j))*(p(2*j) - p(2*ptj)) \
					+ imag(f(mtj))*(p(2*mtj + 1) - p(2*j + 1)) + imag(f(j))*(-p(2*j + 1) + p(2*ptj + 1));
		mds(2*j+1) -=  imag(f(mtj))*(-p(2*mtj) + p(2*j)) + imag(f(j))*(p(2*j) - p(2*ptj)) \
					+ real(f(mtj))*(-p(2*mtj + 1) + p(2*j + 1)) + real(f(j))*(p(2*j + 1) - p(2*ptj + 1));
	}
	else if (ptj==-1) {
		mds(2*j) -=  real(f(mtj))*(-p(2*mtj) + p(2*j)) + imag(f(mtj))*(p(2*mtj + 1) - p(2*j + 1));
		mds(2*j+1) -=  imag(f(mtj))*(-p(2*mtj) + p(2*j)) + real(f(mtj))*(-p(2*mtj + 1) + p(2*j + 1));
	}
	else if (mtj==-1) {
		mds(2*j) -= real(f(j))*(p(2*j) - p(2*ptj)) + imag(f(j))*(-p(2*j + 1) + p(2*ptj + 1));
		mds(2*j+1) -= imag(f(j))*(p(2*j) - p(2*ptj)) + real(f(j))*(p(2*j + 1) - p(2*ptj + 1));
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	4 - nr matrix functions
-------------------------------------------------------------------------------------------------------------------------*/

// ddKineticS_nr
void ddKineticS_nr (const uint& j, const uint& k, const vec& p, const Parameters& pr, const cVec& f, spMat& dds) {
	long mxj = neigh(j,1,-1,pr);
	long mxk = neigh(k,1,-1,pr);
	long pxj = neigh(j,1,1,pr);
	long pxk = neigh(k,1,1,pr);
	
	if (j==k) {
		if (mxj!=-1) {
			dds.coeffRef(2*j,2*k) += + real(f(mxj)) + real(f(j)); 
			dds.coeffRef(2*j,2*k+1) += - imag(f(mxj)) - imag(f(j));
			dds.coeffRef(2*j+1,2*k) += + imag(f(mxj)) + imag(f(j));
			dds.coeffRef(2*j+1,2*k+1) += + real(f(mxj)) + real(f(j));
		}
		else {
			dds.coeffRef(2*j,2*k) +=  real(f(j)); 
			dds.coeffRef(2*j,2*k+1) += - imag(f(j));
			dds.coeffRef(2*j+1,2*k) +=  imag(f(j));
			dds.coeffRef(2*j+1,2*k+1) +=  real(f(j));
		}
		
	}

	if (j==pxk && mxj==k) { // not sure this double check is needed
		dds.coeffRef(2*j,2*pxk) += - real(f(mxj));
		dds.coeffRef(2*j,2*pxk+1) += imag(f(mxj));
		dds.coeffRef(2*j+1,2*pxk) += - imag(f(mxj));
		dds.coeffRef(2*j+1,2*pxk+1) += - real(f(mxj));
	}
	if (j==mxk && pxj==k) {
		dds.coeffRef(2*j,2*mxk) += - real(f(j));
		dds.coeffRef(2*j,2*mxk+1) += + imag(f(j));
		dds.coeffRef(2*j+1,2*mxk) += - imag(f(j));
		dds.coeffRef(2*j+1,2*mxk+1) += - real(f(j));
	}
}

// ddKineticT_nr
void ddKineticT_nr (const uint& j, const uint& k, const vec& p, const Parameters& pr, const cVec& f, spMat& dds) {
long ptj = neigh(j,0,1,pr);
long ptk = neigh(k,0,1,pr);
long mtj = neigh(j,0,-1,pr);
long mtk = neigh(k,0,-1,pr);

	if (j==k) {
		if (mtj!=-1) {
			dds.coeffRef(2*j,2*k) += + real(f(mtj)) + real(f(j)); 
			dds.coeffRef(2*j,2*k+1) += - imag(f(mtj)) - imag(f(j));
			dds.coeffRef(2*j+1,2*k) += + imag(f(mtj)) + imag(f(j));
			dds.coeffRef(2*j+1,2*k+1) += + real(f(mtj)) + real(f(j));
		}
		else {
			dds.coeffRef(2*j,2*k) +=  real(f(j)); 
			dds.coeffRef(2*j,2*k+1) += - imag(f(j));
			dds.coeffRef(2*j+1,2*k) +=  imag(f(j));
			dds.coeffRef(2*j+1,2*k+1) +=  real(f(j));
		}
		
	}

	if (j==ptk && mtj==k) { // not sure this double check is needed
		dds.coeffRef(2*j,2*ptk) += - real(f(mtj));
		dds.coeffRef(2*j,2*ptk+1) += imag(f(mtj));
		dds.coeffRef(2*j+1,2*ptk) += - imag(f(mtj));
		dds.coeffRef(2*j+1,2*ptk+1) += - real(f(mtj));
	}
	if (j==mtk && ptj==k) {
		dds.coeffRef(2*j,2*mtk) += - real(f(j));
		dds.coeffRef(2*j,2*mtk+1) += + imag(f(j));
		dds.coeffRef(2*j+1,2*mtk) += - imag(f(j));
		dds.coeffRef(2*j+1,2*mtk+1) += - real(f(j));
	}
}

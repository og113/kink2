/*-------------------------------------------------------------------------------------------------------------------------
	definitions of functions for newton-raphson quantities
-------------------------------------------------------------------------------------------------------------------------*/

#include <cmath>
#include <complex>
#include "potentials.h"

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

// KineticS
void KineticS (const uint& j, const vec& p, const Parameters& p, const cVec& f, comp& result) {
	
	long pxj = neigh(j,1,1,p); // assuming one spatial direction (the 1-direction)
	if (pxj!=-1) {
		comp cp_j(p(2*j), p(2*j+1));
		comp cp_pxj(p(2*pxj), p(2*pxj+1));
		result += f(j)*pow(cp_pxj-cp_j,2.0)/2.0;
	}
	
}

// KineticT
void KineticT (const uint& j, const vec& p, const Parameters& p, const cVec& f, comp& result) {
	
	long ptj = neigh(j,o,1,p); // taking time to be the 0 direction
	if (ptj!=-1) 
		comp cp_j(p(2*j), p(2*j+1));
		comp cp_ptj(p(2*ptj), p(2*pxj+1));
		result += f(j)*pow(cp_ptj-cp_j,2.0)/2.0;
	}
	
}

/*-------------------------------------------------------------------------------------------------------------------------
	3 - nr vector functions
-------------------------------------------------------------------------------------------------------------------------*/

// mdKineticS
void mdKineticS (const uint& j, const vec& p, const Parameters& p, const cVec& f, vec& mds) {
	long pxj = neigh(j,1,1,p);
	long mxj = neigh(j,1,-1,p);
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


// mdKineticT
void mdKineticT (const uint& j, const vec& p, const Parameters& p, const cVec& f, vec& mds) {
	long ptj = neigh(j,0,1,p);
	long mtj = neigh(j,0,-1,p);
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

// ddKineticS
void ddKineticS (const uint& j, const uint& k, const vec& p, const Parameters& p, const cVec& f, spMat& dds) {
	long mxj = neigh(j,1,-1,p);
	long mxk = neigh(k,1,-1,p);
	long pxj = neigh(j,1,1,p);
	long pxk = neigh(k,1,1,p);
	
	DfkDul
	  real(f(mk))*(-delta(mk,j) + delta(k,j))\
	   + real(f(k))*(delta(k,j) - delta(pk,j))
	 
	DfkDvl 
	  imag(f(mj))*(delta(mj,k) - delta(j,k)) 
	  + imag(f(j))*(-delta(j,k) + delta(pj,k))
	
	dds.coeffRef(2*j,2*k) += real(f(j)) + real(f(mxj)); 
	dds.coeffRef(2*j,2*k+1) += - imag(f(mj)) - imag(f(j));
	dds.coeffRef(2*j+1,2*k) += ;
	dds.coeffRef(2*j+1,2*k+1) += ;
	if (pxk!=-1) {
		dds.coeffRef(2*j,2*pxk) += - real(f(mxj));
		dds.coeffRef(2*j,2*pxk+1) += ;
		dds.coeffRef(2*j+1,2*pxk) += ;
		dds.coeffRef(2*j+1,2*pxk+1) += ;
	}
	else {
	
	}
	if (mxk!=-1) {
		dds.coeffRef(2*j,2*mxk) += - real(f(j))
		dds.coeffRef(2*j,2*mxk+1) += ;
		dds.coeffRef(2*j+1,2*mxk) += ;
		dds.coeffRef(2*j+1,2*mxk+1) += ;
	}
	else {
	
	}
}

// ddKineticT
void ddKineticT (const uint& j, const uint& k, const vec& p, const Parameters& p, const cVec& f, spMat& dds) {
long ptj = neigh(j,0,1,p);
long ptk = neigh(k,0,1,p);
	if (ptj!=-1) {
		

	}
}

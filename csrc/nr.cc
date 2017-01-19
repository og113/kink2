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
	if (pxj!=-1) { // NOT AT ALL SURE THAT THIS IS THE RIGHT CONDITION, NEED TO THINK ABOUT BOUNDARY CONDITIONS
		comp cp_j(p(2*j), p(2*j+1));
		comp cp_pxj(p(2*pxj), p(2*pxj+1));
		
		mds(2*j) -= real(f(j)*(cp_j-cp_pxj));
		mds(2*j+1) -= imag(f(j)*(cp_j-cp_pxj));
	}
	if (mxj!=-1) { // NOT AT ALL SURE THAT THIS IS THE RIGHT CONDITION, NEED TO THINK ABOUT BOUNDARY CONDITIONS
		comp cp_j(p(2*j), p(2*j+1));
		comp cp_mxj(p(2*mxj), p(2*mxj+1));
		
		mds(2*j) -= real(f(mxj)*(cp_j-cp_mxj));
		mds(2*j+1) -= imag(f(mxj)*(cp_j-cp_mxj));
	}
}


// mdKineticT
void mdKineticT (const uint& j, const vec& p, const Parameters& p, const cVec& f, vec& mds) {
	long ptj = neigh(j,1,1,p);
	if (ptj!=-1) {
		comp cp_j(p(2*j), p(2*j+1));
		comp cp_ptj(p(2*ptj), p(2*ptj+1));
		
		mds(2*j) -= real(f(j)*(cp_j-cp_ptj));
		mds(2*j+1) -= imag(f(j)*(cp_j-cp_ptj));
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	4 - nr matrix functions
-------------------------------------------------------------------------------------------------------------------------*/

// ddKineticS
void ddKineticS (const uint& j, const vec& p, const Parameters& p, const cVec& f, spMat& dds) {
	long pxj = neigh(j,1,1,p);
	if (pxj!=-1) {
		comp cp_j(p(2*j), p(2*j+1));
		comp cp_pxj(p(2*pxj), p(2*pxj+1));
		
		dds.coeffRef(2*j,2*j) += ; // is this really what i want? all of these terms?
		dds.coeffRef(2*j,2*j+1) += ;
		dds.coeffRef(2*j+1,2*j) += ;
		dds.coeffRef(2*j+1,2*j) += ;
		dds.coeffRef(2*j,2*pxj) += ;
		dds.coeffRef(2*j,2*pxj+1) += ;
		dds.coeffRef(2*j+1,2*pxj) += ;
		dds.coeffRef(2*j+1,2*pxj) += ;
		dds.coeffRef(2*pxj,2*j) += ;
		dds.coeffRef(2*pxj,2*j+1) += ;
		dds.coeffRef(2*pxj+1,2*j) += ;
		dds.coeffRef(2*pxj+1,2*j) += ;
		dds.coeffRef(2*pxj,2*pxj) += ;
		dds.coeffRef(2*pxj,2*pxj+1) += ;
		dds.coeffRef(2*pxj+1,2*pxj) += ;
		dds.coeffRef(2*pxj+1,2*pxj) += ;
	}
}

// ddKineticT
void ddKineticT (const uint& j, const vec& p, const Parameters& p, const cVec& f, spMat& dds);

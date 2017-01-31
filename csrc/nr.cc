/*-------------------------------------------------------------------------------------------------------------------------
	definitions of functions for newton-raphson quantities
-------------------------------------------------------------------------------------------------------------------------*/

#include <cmath>
#include <complex>
#include "lattice.h"
#include "nr.h"
#include "simple.h"
#include "main.h"

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
void mdKinetic_nr (const uint& j, const uint& dir, const vec& p, const Parameters& pr, const cVec& f, vec& mds\
		, const Complex_nr::Option& opt) {
	long pj = neigh(j,dir,1,pr);
	long mj = neigh(j,dir,-1,pr);
	if (pj!=-1 && mj!=-1) {
		if (opt==Complex_nr::real || opt==Complex_nr::both)
			mds(2*j) -=  real(f(mj))*(-p(2*mj) + p(2*j)) + real(f(j))*(p(2*j) - p(2*pj)) \
						+ imag(f(mj))*(p(2*mj + 1) - p(2*j + 1)) + imag(f(j))*(-p(2*j + 1) + p(2*pj + 1));
		if (opt==Complex_nr::imaginary || opt==Complex_nr::both)
			mds(2*j+1) -=  imag(f(mj))*(-p(2*mj) + p(2*j)) + imag(f(j))*(p(2*j) - p(2*pj)) \
						+ real(f(mj))*(-p(2*mj + 1) + p(2*j + 1)) + real(f(j))*(p(2*j + 1) - p(2*pj + 1));
	}
	else if (pj==-1) {
		if (opt==Complex_nr::real || opt==Complex_nr::both)
			mds(2*j) -=  real(f(mj))*(-p(2*mj) + p(2*j)) + imag(f(mj))*(p(2*mj + 1) - p(2*j + 1));
		if (opt==Complex_nr::imaginary || opt==Complex_nr::both)
			mds(2*j+1) -=  imag(f(mj))*(-p(2*mj) + p(2*j)) + real(f(mj))*(-p(2*mj + 1) + p(2*j + 1));
	}
	else if (mj==-1) {
		if (opt==Complex_nr::real || opt==Complex_nr::both)
			mds(2*j) -= real(f(j))*(p(2*j) - p(2*pj)) + imag(f(j))*(-p(2*j + 1) + p(2*pj + 1));
		if (opt==Complex_nr::imaginary || opt==Complex_nr::both)
			mds(2*j+1) -= imag(f(j))*(p(2*j) - p(2*pj)) + real(f(j))*(p(2*j + 1) - p(2*pj + 1));
	}
}

// mdPotential_nr
void mdPotential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& dv, const cVec& f, vec& mds\
		, const Complex_nr::Option& opt) {
	  if (opt==Complex_nr::real || opt==Complex_nr::both)
	  	mds(2*j) -= real(f(j)*dv(p(2*j) + ii*p(2*j + 1)));
	  if (opt==Complex_nr::imaginary || opt==Complex_nr::both)
	  	mds(2*j+1) -= imag(f(j)*dv(p(2*j) + ii*p(2*j + 1)));
}

/*-------------------------------------------------------------------------------------------------------------------------
	4 - nr matrix functions
-------------------------------------------------------------------------------------------------------------------------*/

// ddKinetic_nr
void ddKinetic_nr (const uint& j, const uint& dir, const vec& p, const Parameters& pr, const cVec& f, spMat& dds\
		, const Complex_nr::Option& opt) {
	long mj = neigh(j,dir,-1,pr);
	long pj = neigh(j,dir,1,pr);
	
	// coincident points
	if (mj!=-1) {
		if (opt==Complex_nr::real || opt==Complex_nr::both) {
			dds.coeffRef(2*j,2*j) += real(f(mj)) + real(f(j)); 
			dds.coeffRef(2*j,2*j+1) += - imag(f(mj)) - imag(f(j));
		}
		if (opt==Complex_nr::imaginary || opt==Complex_nr::both) {
			dds.coeffRef(2*j+1,2*j) += imag(f(mj)) + imag(f(j));
			dds.coeffRef(2*j+1,2*j+1) += real(f(mj)) + real(f(j));
		}
	}
	else {
		if (opt==Complex_nr::real || opt==Complex_nr::both) {
			dds.coeffRef(2*j,2*j) += real(f(j)); 
			dds.coeffRef(2*j,2*j+1) += - imag(f(j));
		}
		if (opt==Complex_nr::imaginary || opt==Complex_nr::both) {
			dds.coeffRef(2*j+1,2*j) += imag(f(j));
			dds.coeffRef(2*j+1,2*j+1) += real(f(j));
		}
	}
	
	// forward neighbour
	if (pj!=-1) {
		if (opt==Complex_nr::real || opt==Complex_nr::both) {
			dds.coeffRef(2*j,2*pj) += - real(f(j));
			dds.coeffRef(2*j,2*pj+1) += imag(f(j));
		}
		if (opt==Complex_nr::imaginary || opt==Complex_nr::both) {
			dds.coeffRef(2*j+1,2*pj) += - imag(f(j));
			dds.coeffRef(2*j+1,2*pj+1) += - real(f(j));
		}
	}
	
	// backward neighbour
	if (mj!=-1) {
		if (opt==Complex_nr::real || opt==Complex_nr::both) {
			dds.coeffRef(2*j,2*mj) += - real(f(mj));
			dds.coeffRef(2*j,2*mj+1) += imag(f(mj));
		}
		if (opt==Complex_nr::imaginary || opt==Complex_nr::both) {
			dds.coeffRef(2*j+1,2*mj) += - imag(f(mj));
			dds.coeffRef(2*j+1,2*mj+1) += - real(f(mj));
		}
	}
}

// ddPotential_nr
void ddPotential_nr (const uint& j, const vec& p, const Parameters& pr, const Potential<comp>& ddv, const cVec& f, spMat& dds\
		, const Complex_nr::Option& opt) {
	
	if (opt==Complex_nr::real || opt==Complex_nr::both) {
		dds.coeffRef(2*j,2*j) += real(f(j)*ddv(p(2*j) + ii*p(2*j + 1))); 
		dds.coeffRef(2*j,2*j+1) += - imag(f(j)*ddv(p(2*j) + ii*p(2*j + 1)));
	}
	if (opt==Complex_nr::imaginary || opt==Complex_nr::both) {
		dds.coeffRef(2*j+1,2*j) += imag(f(j)*ddv(p(2*j) + ii*p(2*j + 1)));
		dds.coeffRef(2*j+1,2*j+1) += real(f(j)*ddv(p(2*j) + ii*p(2*j + 1)));
	}
	
}


/*-------------------------------------------------------------------------------------------------------------------------
	5 - filename functions
-------------------------------------------------------------------------------------------------------------------------*/

// filenameMain
Filename filenameMain(const Parameters& pr, const string& base, const string& subfolder, const string& id, const string& suffix) {
	Filename f = base+"data/pot_"+nts(pr.Pot)+"/"+subfolder+"/"+id\
					+"_N_"+nts(pr.N)\
					+"_Na_"+nts(pr.Na)\
					+"_Nb_"+nts(pr.Nb)\
					+"_Nc_"+nts(pr.Nc)\
					+"_LoR_"+nts(pr.LoR)\
					+"_DE_"+nts(pr.DE)\
					+"_Tb_"+nts(pr.Tb)\
					+"_Theta_"+nts(pr.Theta)\
					+"_Reg_"+nts(pr.Reg)\
					+ suffix;
	return f;
}

// filenameSpatial
Filename filenameSpatial(const Parameters& pr, const string& base, const string& subfolder, const string& id, const string& suffix) {
	Filename f = base+"data/pot_"+nts(pr.Pot)+"/"+subfolder+"/"+id\
					+"_N_"+nts(pr.N)\
					+"_LoR_"+nts(pr.LoR)\
					+"_DE_"+nts(pr.DE)\
					+"_Reg_"+nts(pr.Reg)\
					+ suffix;
	return f;
}

// filenamePeriodic
Filename filenamePeriodic(const Parameters& pr, const string& base, const string& subfolder, const string& id, const string& suffix) {
	Filename f = base+"data/pot_"+nts(pr.Pot)+"/"+subfolder+"/"+id\
					+"_N_"+nts(pr.N)\
					+"_Nb_"+nts(pr.Nb)\
					+"_LoR_"+nts(pr.LoR)\
					+"_DE_"+nts(pr.DE)\
					+"_Tb_"+nts(pr.Tb)\
					+"_Reg_"+nts(pr.Reg)\
					+ suffix;
	return f;
}

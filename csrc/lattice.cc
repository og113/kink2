/*-------------------------------------------------------------------------------------------------------------------------
	definitions for functions and classes for dealing with the representation of the lattice
-------------------------------------------------------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <complex>
#include "error.h"
#include "lattice.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - integer lattice coord fns
	2 - real (or complex) lattice coord fns
	3 - neigh fns
	4 - dx, dt etc fns
	5 - interpolate fns
	6 - vecComplex and vecReal
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. integer lattice coord fns
		- c
		- intCoord
-------------------------------------------------------------------------------------------------------------------------*/

// c - inverse of intCoord - (t,x,p)
lint c(const uint& t, const uint& x, const Parameters& p ) {
	return t+x*p.NT;
}

// c - inverse of intCoord - (t,x,Nt)
lint c(const uint& t, const uint& x, const uint& Nt) {
	return t+x*Nt;
}

// intCoord - (t,x,p)
uint intCoord(const lint& loc, const uint& direc, const Parameters& p) {
	uint x = floor(loc/p.NT);
	if (direc==1) 		return x;
	else 				return loc - x*p.NT;
}

// intCoord - (t,x,Nt)
uint intCoord(const lint& loc, const uint& direc, const uint& Nt) {
	uint x = floor(loc/Nt);
	if (direc==1) 		return x;
	else 				return loc - x*Nt;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. real (or complex) lattice coord fns
		- simpleTime
		- simpleSpace
		- coord
		- coordB
-------------------------------------------------------------------------------------------------------------------------*/

// simpleTime
static comp simpleTime(const uint& t, const Parameters& p) {
	if ( t < p.Na) {
		double temp = (double)t;
		temp -= (double)p.Na;
		return p.b*temp + comp(0.0,1.0)*p.Tb;
	}
	else if (t < (p.Na+p.Nb)) {
		double temp = (double)t;
		temp -= (double)p.Na; //as complex doesn't support (complex double)*integer (though it does support double*integer added to a complex double) - and as int to double seems to cock up here (perhaps because the integers are unsigned)
		return comp(0.0,1.0)*(p.Tb - p.b*temp);
	}
	else {
		double temp = (double)t;
		temp -= (double)p.Na;
		temp -= (double)p.Nb;
		return p.b*(temp+1.0); //the 1.0 is because the corner is part of the vertical contour
	}
}

// simpleSpace
static double simpleSpace (const uint& x, const Parameters& p) {
	if (p.pot==0) {
		cerr << "simpleSpace error: p.pot = 0" << endl;
		return 0;
	}
	if (p.pot!=3) {
		return -p.L/2.0 + (double)x*p.a;
	}
	else {
		return p.r0 + (double)x*p.a;
	}
}

// coord
comp coord(const lint& loc, const int& direction, const Parameters& p) {
	if (direction==0)		return simpleTime(intCoord(loc,0,p),p);
	else if (direction==1)	return simpleSpace(intCoord(loc,1,p),p);
	else return 0.0;
}

// coordB
//gives values of coordinates on section BC
comp coordB(const lint& loc, const int& direction, const Parameters& p)
	{
	comp XcoordB;
	if (direction==0) {
		uint t = intCoord(loc,0,p.Nb);
		double temp = (double)t;
		XcoordB = comp(0.0,1.0)*(p.Tb - p.b*temp);
	}
	if (direction==1) {
		uint x = intCoord(loc,1,p.Nb);
		XcoordB = simpleSpace(x,p);
	}
	return XcoordB;
	}

/*-------------------------------------------------------------------------------------------------------------------------
	3. neigh fns
		- periodic
		- spherical
		- neigh
-------------------------------------------------------------------------------------------------------------------------*/

static long periodic(const lint& locNum, const uint& direction, const int& sign, const uint& xNt, const uint& xNx) //periodic in space but not time, degree refers to the number of neighbours, 1 is for just positive neighbours, 2 is for both
	{
	unsigned int c = intCoord(locNum,direction,xNt);
	if (direction==0) {
		if (sign==1 and c!=(xNt-1))			return locNum+1;
		else if (sign==-1 and c!=0)			return locNum-1;
		else 								return -1;
	}
	else if (c==0 and sign==-1)				return locNum+(xNx-1)*xNt;
	else if (c==(xNx-1) and sign==1)		return locNum-(xNx-1)*(int)xNt;
	else									return locNum+sign*(int)xNt;
	}
	
//spherical system, reflective at r=0, nothing at r=R
static long spherical(const lint& locNum, const uint& direction, const int& sign, const uint& xNt, const uint& xNx) { //periodic in space but not time, degree refers to the number of neighbours, 1 is for just positive neighbours, 2 is for both
	
	unsigned int c = intCoord(locNum,direction,xNt);
	if (direction==0) {
		if (sign==1 and c!=(xNt-1))		return locNum+1;
		else if (sign==-1 and c!=0)		return locNum-1;
		else 							return -1;
	}
	else if (c==0 and sign==-1)			return -1;
	else if (c==(xNx-1) and sign==1)	return -1;
	else								return locNum+sign*(int)xNt;
}

// neigh	
long neigh (const lint& loc, const uint& direction, const int& sign, const Parameters p) {
	if (p.pot==0) {
		cerr << "neigh error: p.pot = 0" << endl;
		return 0;
	}
	else if (p.pot!=3) {
		return periodic(loc,direction,sign,p.NT,p.N);
	}
	else {
		return spherical(loc,direction,sign,p.NT,p.N);
	}
}

// neigh	
long neigh (const lint& loc, const uint& direction, const int& sign, const uint& xNt, const uint& xNx, const uint& pot) {
	if (pot!=3) {
		return periodic(loc,direction,sign,xNt,xNx);
	}
	else {
		return spherical(loc,direction,sign,xNt,xNx);
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	4. dx, dt etc
		- dtFn
		- DtFn
		- dxFn
		- DxFn
-------------------------------------------------------------------------------------------------------------------------*/

// dtFn
comp dtFn (const unsigned int& time, const Parameters& p) {
	return (time<(p.NT-1)? simpleTime(time+1,p)-simpleTime(time,p): 0.0);
}

// DtFn	
comp DtFn (const unsigned int& time, const Parameters& p) {
	if (time==(p.NT-1)) return (simpleTime(p.NT-1,p)-simpleTime(p.NT-2,p))/2.0;
	else if (time==0) return (simpleTime(1,p)-simpleTime(0,p))/2.0;
	else return (simpleTime(time+1,p)-simpleTime(time-1,p))/2.0;
}

// dxFn	
double dxFn (const unsigned int& space, const Parameters& p) {
	if (p.pot==0) {
		cerr << "dxFn error: p.pot = 0" << endl;
		return 0;
	}
	if (p.pot!=3) return p.a;
	return (space<(p.N-1)? p.a: 0.0);
}

// DxFn	
double DxFn (const unsigned int& space, const Parameters& p) {
	if (p.pot==0) {
		cerr << "DxFn error: p.pot = 0" << endl;
		return 0;
	}
	if (p.pot!=3) return p.a;
	return ((space==(p.N-1) || space==0)? p.a/2.0: p.a);
}

/*-------------------------------------------------------------------------------------------------------------------------
	5. interpolate fns
		- interpolate (for real representation of complex 2d vector)
		- interpolateReal (for fundamentally real 2d vector)
		- interpolate1d (for 1d real vector)
-------------------------------------------------------------------------------------------------------------------------*/

// interpolate -  2d (real representation of ) complex vector
vec interpolate(vec vec_old, const Parameters& p_old, const Parameters& p_new) {
	uint N_old = p_old.N, Nt_old = p_old.NT, N_new = p_new.N, Nt_new = p_new.NT;
	uint old_size = vec_old.size();
	if (N_old == N_new && Nt_old == Nt_new) return vec_old;
	if (Nt_new==0 || Nt_old==0) {
		Nt_new = p_new.Nb;
		Nt_old = p_old.Nb;
	}
	if (old_size<2*N_old*Nt_old) {
		if (old_size>=2*N_old*p_old.Nb) {
			Nt_old = p_old.Nb;
			Nt_new = p_new.Nb;
		}
		else {
			cerr << "interpolate error: vec_old.size() = " << old_size << " , 2*N_old*Nt_old = " << 2*N_old*Nt_old << endl;
		}
	}
	uint zero_modes = old_size - 2*N_old*Nt_old;
	vec vec_new (2*N_new*Nt_new+zero_modes);
	
	uint x_new, t_new, x_old, t_old;
	double exact_x_old, exact_t_old, rem_x_old, rem_t_old;
	lint pos;
	int neigh_t, neigh_x, neigh_tx;
	
	for (lint l=0;l<N_new*Nt_new;l++) {
		t_new = intCoord(l,0,Nt_new);
		x_new = intCoord(l,1,Nt_new);
		exact_t_old = t_new*(Nt_old-1.0)/(Nt_new-1.0);
		exact_x_old = x_new*(N_old-1.0)/(N_new-1.0);
		t_old = (unsigned int)exact_t_old;
		x_old = (unsigned int)exact_x_old;
		rem_t_old = exact_t_old;
		rem_t_old -= (double)(t_old);
		rem_x_old = exact_x_old;
		rem_x_old -= (double)(x_old);
		pos = t_old + Nt_old*x_old;
		neigh_t = neigh(pos,0,1,Nt_old,N_old,p_old.pot);
		neigh_x = neigh(pos,1,1,Nt_old,N_old,p_old.pot);
		neigh_tx = neigh(neigh_t,1,1,Nt_old,N_old,p_old.pot);
		if  (t_old<(Nt_old-1) ) {
			if (neigh_t!=-1 && neigh_x!=-1) {
				vec_new(2*l) = (1.0-rem_t_old)*(1.0-rem_x_old)*vec_old(2*pos) \
							+ (1.0-rem_t_old)*rem_x_old*vec_old(2*neigh_x) \
							+ rem_t_old*(1.0-rem_x_old)*vec_old(2*neigh_t) \
							+ rem_t_old*rem_x_old*vec_old(2*neigh_tx);
				vec_new(2*l+1) = (1.0-rem_t_old)*(1.0-rem_x_old)*vec_old(2*pos+1)\
			 				+ (1.0-rem_t_old)*rem_x_old*vec_old(2*neigh_x+1) \
							+ rem_t_old*(1.0-rem_x_old)*vec_old(2*neigh_t+1)\
						 	+ rem_t_old*rem_x_old*vec_old(2*neigh_tx+1);
			}
			else if (neigh_t!=-1) {
				vec_new(2*l) = (1.0-rem_t_old)*vec_old(2*pos) \
							+ rem_t_old*vec_old(2*neigh_t);
				vec_new(2*l+1) = (1.0-rem_t_old)*vec_old(2*pos+1)\
							+ rem_t_old*vec_old(2*neigh_t+1);
			}
			else if(neigh_x!=-1) {
				vec_new(2*l) = (1.0-rem_x_old)*vec_old(2*pos) \
							+ rem_x_old*vec_old(2*neigh_x);
				vec_new(2*l+1) = (1.0-rem_x_old)*vec_old(2*pos+1)\
							+ rem_x_old*vec_old(2*neigh_x+1);
			}
			else {
				vec_new(2*l) = vec_old(2*pos);
				vec_new(2*l+1) = vec_old(2*pos+1);				
			}
		}
	}
	for (uint l=0; l<zero_modes; l++) {
		vec_new(2*N_new*Nt_new+l) = vec_old(2*N_old*Nt_old+l);
	}
	return vec_new;
}

// interpolate -  2d (real representation of ) complex vector
cVec interpolate(cVec vec_old, const Parameters& p_old, const Parameters& p_new) {
	uint N_old = p_old.N, Nt_old = p_old.NT, N_new = p_new.N, Nt_new = p_new.NT;
	uint old_size = vec_old.size();
	if (N_old == N_new && Nt_old == Nt_new) return vec_old;
	if (Nt_new==0 || Nt_old==0) {
		Nt_new = p_new.Nb;
		Nt_old = p_old.Nb;
	}
	if (old_size<N_old*Nt_old) {
		if (old_size>=N_old*p_old.Nb) {
			Nt_old = p_old.Nb;
			Nt_new = p_new.Nb;
		}
		else {
			cerr << "interpolate error: vec_old.size() = " << old_size << " , N_old*Nt_old = " << N_old*Nt_old << endl;
		}
	}
	uint zero_modes = old_size - N_old*Nt_old;
	cVec vec_new (N_new*Nt_new+zero_modes);
	
	uint x_new, t_new, x_old, t_old;
	double exact_x_old, exact_t_old, rem_x_old, rem_t_old;
	lint pos;
	int neigh_t, neigh_x, neigh_tx;
	
	for (lint l=0;l<N_new*Nt_new;l++) {
		t_new = intCoord(l,0,Nt_new);
		x_new = intCoord(l,1,Nt_new);
		exact_t_old = t_new*(Nt_old-1.0)/(Nt_new-1.0);
		exact_x_old = x_new*(N_old-1.0)/(N_new-1.0);
		t_old = (unsigned int)exact_t_old;
		x_old = (unsigned int)exact_x_old;
		rem_t_old = exact_t_old;
		rem_t_old -= (double)(t_old);
		rem_x_old = exact_x_old;
		rem_x_old -= (double)(x_old);
		pos = t_old + Nt_old*x_old;
		neigh_t = neigh(pos,0,1,Nt_old,N_old,p_old.pot);
		neigh_x = neigh(pos,1,1,Nt_old,N_old,p_old.pot);
		neigh_tx = neigh(neigh_t,1,1,Nt_old,N_old,p_old.pot);
		if  (t_old<(Nt_old-1) ) {
			if (neigh_t!=-1 && neigh_x) {
				vec_new(l) = (1.0-rem_t_old)*(1.0-rem_x_old)*vec_old(pos) \
							+ (1.0-rem_t_old)*rem_x_old*vec_old(neigh_x) \
							+ rem_t_old*(1.0-rem_x_old)*vec_old(neigh_t) \
							+ rem_t_old*rem_x_old*vec_old(neigh_tx);
			}
			else if (neigh_t!=-1) {
				vec_new(l) = (1.0-rem_t_old)*vec_old(pos) \
							+ rem_t_old*vec_old(neigh_t);
			}
			else if(neigh_x!=-1) {
				vec_new(l) = (1.0-rem_x_old)*vec_old(pos) \
							+ rem_x_old*vec_old(neigh_x);
			}
			else {
				vec_new(l) = vec_old(pos);				
			}
		}
	}
	for (uint l=0; l<zero_modes; l++) {
		vec_new(N_new*Nt_new+l) = vec_old(N_old*Nt_old+l);
	}
	return vec_new;
}
	
// interpolateReal -  2d real vector function
vec interpolateReal(vec vec_old, const Parameters& p_old, const Parameters& p_new) {
	uint N_old = p_old.N, Nt_old = p_old.NT, N_new = p_new.N, Nt_new = p_new.NT;
	uint old_size = vec_old.size();
	if (N_old == N_new && Nt_old == Nt_new) return vec_old;
	if (Nt_new==0 || Nt_old==0) {
		Nt_new = p_new.Nb;
		Nt_old = p_old.Nb;
	}
	if (old_size<N_old*Nt_old) {
		if (old_size>=N_old*p_old.Nb) {
			Nt_old = p_old.Nb;
			Nt_new = p_new.Nb;
		}
		else {
			cerr << "interpolateReal error: vec_old.size() = " << old_size << " , N_old*Nt_old = " << N_old*Nt_old << endl;
		}
	}
	vec vec_new(N_new*Nt_new);
	uint x_new, t_new, x_old, t_old;
	double exact_x_old, exact_t_old, rem_x_old, rem_t_old;
	lint pos;
	int neigh_t, neigh_x, neigh_tx;	
	
	for (lint l=0;l<N_new*Nt_new;l++) {
		t_new = intCoord(l,0,Nt_new);
		x_new = intCoord(l,1,Nt_new);
		exact_t_old = t_new*(Nt_old-1.0)/(Nt_new-1.0);
		exact_x_old = x_new*(N_old-1.0)/(N_new-1.0);
		t_old = (unsigned int)exact_t_old;
		x_old = (unsigned int)exact_x_old;
		rem_t_old = exact_t_old;
		rem_t_old -= (double)(t_old);
		rem_x_old = exact_x_old;
		rem_x_old -= (double)(x_old);
		pos = t_old + Nt_old*x_old;
		neigh_t = neigh(pos,0,1,Nt_old,N_old,p_old.pot);
		neigh_x = neigh(pos,1,1,Nt_old,N_old,p_old.pot);
		neigh_tx = neigh(neigh_t,1,1,Nt_old,N_old,p_old.pot);
		if  (t_old<(Nt_old-1) ) {
			if (neigh_x!=-1 && neigh_t!=-1) {
				vec_new(l) = (1.0-rem_t_old)*(1.0-rem_x_old)*vec_old(pos) \
							+ (1.0-rem_t_old)*rem_x_old*vec_old(neigh_x) \
							+ rem_t_old*(1.0-rem_x_old)*vec_old(neigh_t) \
							+ rem_t_old*rem_x_old*vec_old(neigh_tx);
			}
			else if (neigh_t!=-1) {
				vec_new(l) = (1.0-rem_t_old)*vec_old(pos) \
							+ rem_t_old*vec_old(neigh_t);
			}
			else if (neigh_x!=-1) {
				vec_new(l) = (1.0-rem_x_old)*vec_old(pos) \
							+ rem_x_old*vec_old(neigh_x);
			}
			else {
				vec_new(l) = vec_old(pos);
			}
		}
	}
	return vec_new;
}
	
// interpolate1d	
vec interpolate1d(vec vec_old, const unsigned int & N_old, const unsigned int & N_new) {
	if (N_old == N_new) return vec_old;
	uint old_size = vec_old.size();
	if (old_size<N_old) cerr << "interpolate error, vec_old.size() = " << old_size << " , N_old = " << N_old << endl;
	vec vec_new (N_new);
	
	uint x_old;
	double exact_x_old, rem_x_old;
	
	for (lint l=0;l<N_new;l++) {
		exact_x_old = l*(N_old-1.0)/(N_new-1.0);
		x_old = (unsigned int)exact_x_old;
		rem_x_old = exact_x_old - (double)(x_old);
		if  (x_old<(N_old-1)) vec_new[l] = (1.0-rem_x_old)*vec_old[x_old] + rem_x_old*vec_old[x_old+1];
		else vec_new[l] = vec_old[x_old];
	}
	return vec_new;
}

// interpolate1d	
cVec interpolate1d(cVec vec_old, const unsigned int & N_old, const unsigned int & N_new) {
	if (N_old == N_new) return vec_old;
	uint old_size = vec_old.size();
	if (old_size<N_old) cerr << "interpolate error, vec_old.size() = " << old_size << " , N_old = " << N_old << endl;
	cVec vec_new (N_new);
	
	uint x_old;
	double exact_x_old, rem_x_old;
	
	for (lint l=0;l<N_new;l++) {
		exact_x_old = l*(N_old-1.0)/(N_new-1.0);
		x_old = (unsigned int)exact_x_old;
		rem_x_old = exact_x_old - (double)(x_old);
		if  (x_old<(N_old-1)) vec_new[l] = (1.0-rem_x_old)*vec_old[x_old] + rem_x_old*vec_old[x_old+1];
		else vec_new[l] = vec_old[x_old];
	}
	return vec_new;
}

/*-------------------------------------------------------------------------------------------------------------------------
	6. vecComplex etc
		- vecComplex
		- vecReal
-------------------------------------------------------------------------------------------------------------------------*/

//complexify a real vector - tDim
cVec vecComplex(vec realVec, const uint & tDim) {
	cVec complexVec(tDim);
	if (realVec.size() >= (2*tDim)) {
		for (uint l=0; l<tDim; l++) {
			complexVec(l) = realVec(2*l) + comp(0.0,1.0)*realVec(2*l+1);
		}
	}
	else cerr << "vecComplex error";
	return complexVec;
}

//complexify a real vector - parameters
cVec vecComplex(vec realVec, const Parameters& ps) {
	cVec complexVec(ps.N*ps.NT);
	if (realVec.size() >= (2*ps.N*ps.NT)) {
		for (uint l=0; l<ps.N*ps.NT; l++) {
			complexVec(l) = realVec(2*l) + comp(0.0,1.0)*realVec(2*l+1);
		}
	}
	else if (realVec.size() >= (2*ps.N*ps.Nb)) {
		complexVec.resize(ps.N*ps.Nb);
		for (uint l=0; l<ps.N*ps.Nb; l++) {
			complexVec(l) = realVec(2*l) + comp(0.0,1.0)*realVec(2*l+1);
		}
	}
	else cerr << "vecComplex error";
	return complexVec;
}
	
//make a complex vector real - tDim
vec vecReal(cVec complexVec, const uint &  tDim) {
	vec realVec(2*tDim);
	if (complexVec.size() == tDim) {
		for (uint l=0; l<tDim; l++) {
			realVec(2*l) = real(complexVec(l));
			realVec(2*l+1) = imag(complexVec(l));
		}
	}
	else cerr << "vecReal error";
	return realVec;
}

//make a complex vector real - parameters
vec vecReal(cVec complexVec, const Parameters&  p) {
	vec realVec(2*p.N*p.NT);
	if (complexVec.size() == p.N*p.NT) {
		for (uint l=0; l<p.N*p.NT; l++) {
			realVec(2*l) = real(complexVec(l));
			realVec(2*l+1) = imag(complexVec(l));
		}
	}
	else if (complexVec.size() == p.N*p.Nb) {
		realVec.resize(2*p.N*p.Nb);
		for (uint l=0; l<p.N*p.Nb; l++) {
			realVec(2*l) = real(complexVec(l));
			realVec(2*l+1) = imag(complexVec(l));
		}
	}
	else cerr << "vecReal error";
	return realVec;
}

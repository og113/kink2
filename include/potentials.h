/*-------------------------------------------------------------------------------------------------------------------------
	declarations for functions giving potentials and their derivatives
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __POTENTIALS_H_INCLUDED__
#define __POTENTIALS_H_INCLUDED__

#include <complex>

using namespace std;

typedef complex<double> comp;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - small parameters structs
	2 - Potential class
	3 - Vs
	4 - dVs
	5 - ddVs
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. small structs containing required parameters
		- params_for_V
	
	also defined some extern instances
		- params_for_V - paramsV, paramsV0
	N.B. the extern word is to prevent problems with multiple definitions.
-------------------------------------------------------------------------------------------------------------------------*/

struct params_for_V {
	double epsi;
	double aa;
};

/*-------------------------------------------------------------------------------------------------------------------------
	 2. potential class
		- Potential
-------------------------------------------------------------------------------------------------------------------------*/

typedef comp(*PotentialType)(const comp&, const struct params_for_V&);

class Potential {
public:
	Potential(const params_for_V&, PotentialType);
	~Potential() {}
	comp operator()(const comp&) const;
	comp operator()(const comp&, const params_for_V&) const;
	void setParams(const params_for_V&);
private:
	params_for_V parameters;
	PotentialType potential;
};

/*-------------------------------------------------------------------------------------------------------------------------
	3. potential functions
		- V1
		- Z
		- V2
		- V3
		- VrFn
	
	NB - gsl doesn't like const T &; default parameters set
	template <class T>
-------------------------------------------------------------------------------------------------------------------------*/

//V1
template <class T> T V1 (const T&, const struct params_for_V&);

//Z, for V2
template <class T> T Z (const T&);
	
//V2
template <class T> T V2 (const T&, const struct params_for_V&);

//V3
template <class T> T V3 (const T&, const struct params_for_V&);

//Vr	
comp VrFn (const comp & phi, const double & minimaL, const double & minimaR);

/*-------------------------------------------------------------------------------------------------------------------------
	4. first derivatives of potential functions
		- dV1
		- dZ
		- dV2
		- dV3
		- dVrFn

-------------------------------------------------------------------------------------------------------------------------*/
//dV1
template <class T> T dV1 (const T&, const struct params_for_V&);

//dZ for dV2
template <class T> T dZ (const T&);
	
//dV2
template <class T> T dV2 (const T&, const struct params_for_V&);

//dV3
template <class T> T dV3 (const T&, const struct params_for_V&);
	
//dVr
comp dVrFn (const comp & phi, const double & minimaL, const double & minimaR);
	
/*-------------------------------------------------------------------------------------------------------------------------
	5. second derivatives of potential functions
		- ddV1
		- ddZ
		- ddV2
		- ddV3
		- ddVrFn
	
-------------------------------------------------------------------------------------------------------------------------*/

//ddV1
template <class T> T ddV1 (const T&, const struct params_for_V&);

//ddZ for ddV2
template <class T> T ddZ (const T&);

//ddV2
template <class T> T ddV2 (const T&, const struct params_for_V&);

//ddV3
template <class T> T ddV3 (const T&, const struct params_for_V&);

//ddVr	
comp ddVrFn (const comp & phi, const double & minimaL, const double & minimaR);

#endif // __POTENTIALS_H_INCLUDED__

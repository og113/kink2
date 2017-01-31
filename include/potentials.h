/*-------------------------------------------------------------------------------------------------------------------------
	declarations for functions giving potentials and their derivatives
-------------------------------------------------------------------------------------------------------------------------*/

#ifndef __POTENTIALS_H_INCLUDED__
#define __POTENTIALS_H_INCLUDED__

#include <complex>
#include "parameters.h"

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
	6 - dddVs
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
	void setFromParameters(const Parameters&);
};

/*-------------------------------------------------------------------------------------------------------------------------
	 2. potential class
		- Potential
-------------------------------------------------------------------------------------------------------------------------*/

template <class T>
class Potential {
public:
	typedef params_for_V ParamsStruct;
	typedef T(*PotentialType)(const T&, const ParamsStruct &);
	Potential();
	Potential (PotentialType, const Parameters&);
	Potential(PotentialType, const ParamsStruct&);
	void operator()(PotentialType, const ParamsStruct&);
	void operator()(PotentialType, const Parameters&);
	~Potential() {}
	T operator()(const T&) const;
	T operator()(const T&, const ParamsStruct&) const;
	void setParams(const ParamsStruct&);
	void setParams(const Parameters&);
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

//VrFn	
comp VrFn (const comp & phi, const double & minimaL, const double & minimaR);

//Vreg
template <class T> T Vreg (const T&, const struct params_for_V&);

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

//dVreg
template <class T> T dVreg (const T&, const struct params_for_V&);
	
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

//ddVreg
template <class T> T ddVreg (const T&, const struct params_for_V&);

/*-------------------------------------------------------------------------------------------------------------------------
	6. third derivatives of potential functions
		- dddV1
		- dddZ
		- dddV2
	
-------------------------------------------------------------------------------------------------------------------------*/

//dddV1
template <class T> T dddV1 (const T&, const struct params_for_V&);

//dddZ for dddV2
template <class T> T dddZ (const T&);

//dddV2
template <class T> T dddV2 (const T&, const struct params_for_V&);

#endif // __POTENTIALS_H_INCLUDED__

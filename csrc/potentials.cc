/*-------------------------------------------------------------------------------------------------------------------------
	definitions of functions giving potentials and their derivatives
-------------------------------------------------------------------------------------------------------------------------*/

#include <cmath>
#include <complex>
#include "potentials.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - params_for_V
	2 - Potential
	3 - Vs
	4 - dVs
	5 - ddVs
	6 - dddVs
	7 - explicit template instantiation
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. params_for_V
		- setFromParameters
-------------------------------------------------------------------------------------------------------------------------*/

void params_for_V::setFromParameters(const Parameters& p) {
	epsi = p.epsilon;
	aa = p.A;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. potential class member function
		- constructors
		- operator()
		- setParams
-------------------------------------------------------------------------------------------------------------------------*/

template <class T>
Potential<T>::Potential(): parameters(), potential(){}

template <class T>
Potential<T>::Potential(PotentialType v, const Parameters& p) {
	parameters.setFromParameters(p);
	potential = v;
}

template <class T>
Potential<T>::Potential(PotentialType v, const ParamsStruct& p) {
	parameters = p;
	potential = v;
}

template <class T>
void Potential<T>::operator()(PotentialType v, const Parameters& p) {
	parameters.setFromParameters(p);
	potential = v;
}

template <class T>
void Potential<T>::operator()(PotentialType v, const params_for_V& p) {
	parameters = p;
	potential = v;
}

template <class T>
T Potential<T>::operator()(const T& p) const {
	return potential(p,parameters);
}

template <class T>
T Potential<T>::operator()(const T& p, const params_for_V& v) const {
	return potential(p,v);
}

template <class T>
void Potential<T>::setParams(const params_for_V& p) {
	parameters = p;
}

template <class T>
void Potential<T>::setParams(const Parameters& p) {
	parameters.setFromParameters(p);
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. potential functions
		- V1
		- Z
		- V2
		- V3
		- VrFn (comp)
	
-------------------------------------------------------------------------------------------------------------------------*/

//V1
template <class T> T V1 (const T& phi, const struct params_for_V& params) {
	return pow(pow(phi,2)-1.0,2.0)/8.0 - params.epsi*(phi-1.0)/2.0;
}

//Z, for V2
template <class T> T Z (const T& x) {
	return exp(-pow(x,2.0))*(x + pow(x,3.0) + pow(x,5.0));
}
	
//V2
template <class T> T V2 (const T& x, const struct params_for_V& params) {
	double epsi = params.epsi;
	double aa = params.aa;
	T y = (x-1.0)/aa;
	return 0.5*pow(x+1.0,2.0)*(1.0-epsi*Z(y));
}

//V3
template <class T> T V3 (const T& phi, const struct params_for_V& params) {
	return pow(phi,2.0)/2.0 - pow(phi,4.0)/4.0/pow(params.epsi,2.0); //epsi here takes the role of r
}

//Vr	
comp VrFn (const comp & phi, const double & minimaL, const double & minimaR) {
	return pow(phi-minimaL,4.0)*pow(phi-minimaR,4.0)/4.0;
}

/*-------------------------------------------------------------------------------------------------------------------------
	4. first derivatives of potential functions
		- dV1
		- dZ
		- dV2
		- dV3
		- dVrFn (comp)

-------------------------------------------------------------------------------------------------------------------------*/

//dV1
template <class T> T dV1 (const T& phi, const struct params_for_V& params) {
	return phi*(pow(phi,2)-1.0)/2.0 - params.epsi/2.0;
}

//dZ for dV2
template <class T> T dZ (const T& phi) {
	return exp(-pow(phi,2.0))*( 1.0 + pow(phi,2.0) + 3.0*pow(phi,4.0) - 2.0*pow(phi,6.0) );
}
	
//dV2
template <class T> T dV2 (const T& x, const struct params_for_V& params) {
	double epsi = params.epsi;
	double aa = params.aa;
	T y = (x-1.0)/aa;
	return (x+1.0)*(1.0-epsi*Z(y)) - (epsi/2.0/aa)*pow(x+1.0,2.0)*dZ(y);
}

//dV3
template <class T> T dV3 (const T& phi, const struct params_for_V& params) {
	double epsi = params.epsi;
	return phi - pow(phi,3.0)/pow(epsi,2.0);
}
	
//dVr
comp dVrFn (const comp & phi, const double & minimaL, const double & minimaR) {
	return pow(phi-minimaL,3.0)*pow(phi-minimaR,4.0) + pow(phi-minimaL,4.0)*pow(phi-minimaR,3.0);
}
	
/*-------------------------------------------------------------------------------------------------------------------------
	5. second derivatives of potential functions
		- dV1
		- dZ
		- dV2
		- dV3
		- dVrFn (comp)
	
-------------------------------------------------------------------------------------------------------------------------*/

//ddV1
template <class T> T ddV1 (const T& phi, const struct params_for_V& params) {
	return (3.0*pow(phi,2)-1.0)/2.0;
}

//ddZ for ddV2
template <class T> T ddZ (const T& x) {
	return exp(-pow(x,2.0))*2.0*pow(x,3.0)*(5.0 - 9.0*pow(x,2.0) + 2.0*pow(x,4.0));
}

//ddV2
template <class T> T ddV2 (const T& x, const struct params_for_V& params) {
	double epsi = params.epsi;
	double aa = params.aa;
	T y = (x-1.0)/aa;
	return 1.0-epsi*Z(y) - (2.0*epsi/aa)*(x+1.0)*dZ(y)\
					- (epsi/2.0/pow(aa,2.0))*pow(x+1.0,2.0)*ddZ(y);
}

//ddV3
template <class T> T ddV3 (const T& phi, const struct params_for_V& params) {
	return 1.0 - 3.0*pow(phi,2.0)/pow(params.epsi,2.0);
}

//ddVr	
comp ddVrFn (const comp & phi, const double & minimaL, const double & minimaR) {
	return 3.0*pow(phi-minimaL,2.0)*pow(phi-minimaR,4.0) + 8.0*pow(phi-minimaL,3.0)*pow(phi-minimaR,3.0)\
				+ 3.0*pow(phi-minimaL,4.0)*pow(phi-minimaR,2.0);
}

/*-------------------------------------------------------------------------------------------------------------------------
	6. third derivatives of potential functions
		- dddV1
		- dddZ
		- dddV2
	
-------------------------------------------------------------------------------------------------------------------------*/

//dddV1
template <class T> T dddV1 (const T& phi, const struct params_for_V& params) {
	return 3.0*phi;
}

//dddZ for dddV2
template <class T> T dddZ (const T& x) {
	return exp(-pow(x,2.0))*2.0*pow(x,2.0)*( 15.0 - 55.0*pow(x,2.0) + 32.0*pow(x,4.0) - 4.0*pow(x,6.0));
}

//dddV2
template <class T> T dddV2 (const T& x, const struct params_for_V& params) {
	double epsi = params.epsi;
	double aa = params.aa;
	T y = (x-1.0)/aa;
	return -3.0*(epsi/aa)*dZ(y) - 3.0*(x+1.0)*(epsi/pow(aa,2.0))*ddZ(y) \
			- 0.5*pow(x+1.0,2.0)*(epsi/pow(aa,3.0))*dddZ(y);
}

/*-------------------------------------------------------------------------------------------------------------------------
	7. explicit template instantiation
		- double
		- complex<double>
-------------------------------------------------------------------------------------------------------------------------*/

template double V1<double>(const double&, const struct params_for_V&);
template double Z<double>(const double&);
template double V2<double>(const double&, const struct params_for_V&);
template double V3<double>(const double&, const struct params_for_V&);
template double dV1<double>(const double&, const struct params_for_V&);
template double dZ<double>(const double&);
template double dV2<double>(const double&, const struct params_for_V&);
template double dV3<double>(const double&, const struct params_for_V&);
template double ddV1<double>(const double&, const struct params_for_V&);
template double ddZ<double>(const double&);
template double ddV2<double>(const double&, const struct params_for_V&);
template double ddV3<double>(const double&, const struct params_for_V&);

template comp V1<comp>(const comp&, const struct params_for_V&);
template comp Z<comp>(const comp&);
template comp V2<comp>(const comp&, const struct params_for_V&);
template comp V3<comp>(const comp&, const struct params_for_V&);
template comp dV1<comp>(const comp&, const struct params_for_V&);
template comp dZ<comp>(const comp&);
template comp dV2<comp>(const comp&, const struct params_for_V&);
template comp dV3<comp>(const comp&, const struct params_for_V&);
template comp ddV1<comp>(const comp&, const struct params_for_V&);
template comp ddZ<comp>(const comp&);
template comp ddV2<comp>(const comp&, const struct params_for_V&);
template comp ddV3<comp>(const comp&, const struct params_for_V&);

template class Potential<double>;
template class Potential<comp>;

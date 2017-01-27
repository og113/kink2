/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for the Angle class
 -------------------------------------------------------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include "simple.h"
#include "angle.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. Angle
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. Angle
		- constructor, void
		- constructor, Angle
		- constructor, double
		- conversion, double
		- copy
		- operator=
		- operator+=
		- operator-=
		- operator*=
		- operator/=
		- isBetween
-------------------------------------------------------------------------------------------------------------------------*/

// constructor, void
Angle::Angle(): Value() {}

// constructor, Angle
Angle::Angle(const Angle& a): Value() {
	copy(a);
}

// constructor, double
Angle::Angle(const double& a): Value() {
	Value = mod(a,Min,Max);
}

// converstion, double
Angle::operator double() const {
	return Value;
}

// copy
void Angle::copy(const Angle& a) {
	Value = (double)a;
}

// operator=
Angle& Angle::operator=(const Angle& rhs) {
	copy(rhs);
	return *this;
}

// operator=
Angle& Angle::operator=(const double& rhs) {
	Value = mod(rhs,Min,Max);
	return *this;
}

// operator+=
Angle& Angle::operator+=(const Angle& rhs) {
	double temp = Value + (double)rhs;
	Value = mod(temp,Min,Max);
	return *this;
}

// operator+=
Angle& Angle::operator+=(const double& rhs) {
	double temp = Value + rhs;
	Value = mod(temp,Min,Max);
	return *this;
}

// operator-=
Angle& Angle::operator-=(const Angle& rhs) {
	double temp = Value - (double)rhs;
	Value = mod(temp,Min,Max);
	return *this;
}

// operator-=
Angle& Angle::operator-=(const double& rhs) {
	double temp = Value - rhs;
	Value = mod(temp,Min,Max);
	return *this;
}

// operator*=
Angle& Angle::operator*=(const double& rhs) {
	double temp = Value*rhs;
	Value = mod(temp,Min,Max);
	return *this;
}

// operator/=
Angle& Angle::operator/=(const double& rhs) {
	double temp = Value/rhs;
	Value = mod(temp,Min,Max);
	return *this;
}

// isBetween
bool Angle::isBetween(const Angle& a, const Angle& b) const {
	double h = ((double)a<(double)b? b: a);
	double l = ((double)a<(double)b? a: b);
	double difference = h-l;
	if (difference<=PI && Value<h && Value>l)
		return true;
	else if (difference<=PI && Value>h && Value<l)
		return true;
	else
		return false;
}

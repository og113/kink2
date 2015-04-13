/*-------------------------------------------------------------------------------------------------------------------------
 	declarations for the Angle class and dependencies
 -------------------------------------------------------------------------------------------------------------------------*/
 
#ifndef __ANGLE_H_INCLUDED__
#define __ANGLE_H_INCLUDED__

#include <cmath>
#include <iostream>
#include "simple.h"


using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. Angle
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. declarations for the Angle class
-------------------------------------------------------------------------------------------------------------------------*/

class Angle {
public:
	Angle();
	Angle(const Angle&);
	Angle(const double&);
	~Angle() {}
	operator double() const;
	Angle& operator=(const Angle&);
	Angle& operator=(const double&);
	Angle& operator+=(const Angle&);
	Angle& operator+=(const double&);
	Angle& operator-=(const Angle&);
	Angle& operator-=(const double&);
	Angle& operator*=(const double&);
	Angle& operator/=(const double&);
	bool isBetween(const Angle&, const Angle&) const;
private:
	double 			Value;
	static double	Max;
	static double	Min;
	void 			copy(const Angle&);
};

double Angle::Max = pi;
double Angle::Min = -pi;

#endif // __ANGLE_H_INCLUDED__

/*-------------------------------------------------------------------------------------------------------------------------
 	declarations for the Stepper class and dependencies
 -------------------------------------------------------------------------------------------------------------------------*/
 
#ifndef __STEPPER_H_INCLUDED__
#define __STEPPER_H_INCLUDED__

#include <string>
#include <vector>
#include <utility> //for pair
#include <iostream>
#include "error.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. Point2d
	2. Stepper
	
n.b. stepper defined in 2d
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. declarations for the simple Point2d class
		- Point2d
		- operator<<
just a collection of two numbers (x,y)
-------------------------------------------------------------------------------------------------------------------------*/

// Point2d
class Point2d {
public:
	Point2d(): X(), Y() {}
	Point2d(const double& XX, const double& YY): X(XX), Y(YY) {}
	Point2d(const Point2d&);
	~Point2d() {}
	Point2d& operator=(const Point2d& rhs);
	void operator()(const double&, const double&);
	void operator()(const Point2d&);
	double X;
	double Y;
private:
	void copy(const Point2d&);
};

// operator<<
ostream& operator<<(ostream&,const Point2d&);

/*-------------------------------------------------------------------------------------------------------------------------
	2. Stepper
		- FxyPair
		- StepperOptions
		- Stepper
-------------------------------------------------------------------------------------------------------------------------*/

// FxyPair
typedef pair<Point2d,double> FxyPair;

// StepperOptions
struct StepperOptions{
	double 			epsi_x;
	double			epsi_y;
	double 			angle0;
	enum			stepTypeList {straight=1, constSimple=2, lagrange=3};
	stepTypeList 	stepType;
	bool			directed;
};

// Stepper
class Stepper {
public:
	Stepper(const StepperOptions& sto, const double& X, const double& Y);
	Stepper(const StepperOptions& sto, const Point2d& P);
	Stepper(const StepperOptions& sto);
	~Stepper() {}
	void 		setStart(const double& X, const double& Y);
	void 		setStart(const Point2d& P);
	void 		step();
	void 		addResult(const double& f);
	void		addResult(const double& f, const double& e, const double& n);
	Point2d 	point() const;
	double		x() const;
	double		y() const;
	uint		offset() const;
	uint		size() const;
	bool		keep() const;
private:
	StepperOptions 	opts;
	vector<FxyPair> f_xy;
	double			angle;
};

#endif // __STEPPER_H_INCLUDED__

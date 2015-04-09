/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for the Stepper class and dependencies
 -------------------------------------------------------------------------------------------------------------------------*/

#include <cmath>
#include <string>
#include <vector>
#include <utility> //for pair
#include <iostream>
#include <iomanip>
#include "error.h"
#include "simple.h"
#include "stepper.h"

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
	1. Point2d
		- copy
		- copy constructor
		- operator=
		- operator()
		- operator<<
-------------------------------------------------------------------------------------------------------------------------*/

// copy
void Point2d::copy(const Point2d& p) {
	X = p.X;
	Y = p.Y;
}

// copy constructor
Point2d::Point2d(const Point2d& p) {
	copy(p);
}

// operator=
Point2d& Point2d::operator=(const Point2d& rhs) {
	copy(rhs);
	return *this;
}

// operator()
void Point2d::operator()(const double& XX, const double& YY) {
	X = XX;
	Y = YY;
}

// operator()
void Point2d::operator()(const Point2d& p) {
	copy(p);
}

// operator<<
ostream& operator<<(ostream& os,const Point2d& p) {
	os << left << setw(15) << p.X << setw(15) << p.Y;
	return os;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. Stepper
		- constructors
-------------------------------------------------------------------------------------------------------------------------*/

// constructor
Stepper::Stepper(const StepperOptions& sto, const double& X, const double& Y): opts(sto), f_xy() {
	Point2d P(X,Y);
	FxyPair toAdd(P,0.0);
	f_xy.push_back(toAdd);
}

// constructor
Stepper::Stepper(const StepperOptions& sto, const Point2d& P): opts(sto), f_xy() {
	FxyPair toAdd(P,0.0);
	f_xy.push_back(toAdd);
}

// constructor
Stepper::Stepper(const StepperOptions& sto): opts(sto), f_xy() {}

// set Start
void Stepper::setStart(const double& X, const double& Y) {
	f_xy.clear();
	Point2d P(X,Y);
	FxyPair toAdd(P,0.0);
	f_xy.push_back(toAdd);
}

// set Start
void Stepper::setStart(const Point2d& P) {
	f_xy.clear();
	FxyPair toAdd(P,0.0);
	f_xy.push_back(toAdd);
}

// x()
double Stepper::x() const {
	return ((f_xy.back()).first).X;
}

// y()
double Stepper::y() const {
	return ((f_xy.back()).first).Y;
}

// offset()
uint Stepper::offset() const {
	uint sizem1 = size()-1;
	if (opts.constant) {
		if (sizem1==3) {
			return 3;
		}
		else if (sizem1==2 || (sizem1%2!=0 && sizem1>3)) {
			return 2;
		}
		else
			return 1;	
	}
	else return 1;
}

// keep()
bool Stepper::keep() const {
	if (opts.constant) {
		if (size()>3 && size()%2==0) {
			return true;
		}
		else
			return false;	
	}
	else
		return true;
}

// point()
Point2d Stepper::point() const{
	return (f_xy.back()).first;
}

// size()
uint Stepper::size() const {
	return f_xy.size();
}

// stepVariables
void Stepper::step() {
	if (size()==0) {
		cerr << "Stepper error: cannot step before giving initial step" << endl;
		return;
	}
	if (abs((f_xy.back()).second)<MIN_NUMBER) {
		cerr << "Stepper error: cannot step before giving result of previous step" << endl;
		return;
	}
	double x_old = x();
	double y_old = y();
	if (opts.constant) {
		if (size()==3) {
			x_old = ((f_xy[0]).first).X;
			y_old = ((f_xy[0]).first).Y;
		}
		else if (size()==2 || (size()%2!=0 && size()>3)) {
			x_old =((f_xy[size()-2]).first).X;
			y_old = ((f_xy[size()-2]).first).Y;
		}		
	}
	double x_new = x_old + opts.epsi_x*cos(opts.angle);
	double y_new = y_old + opts.epsi_y*sin(opts.angle);
	Point2d P(x_new,y_new);
	FxyPair toAdd(P,0.0);
	f_xy.push_back(toAdd);
}

// addResult
void Stepper::addResult(const double& f) {
	(f_xy.back()).second = f;
	if (opts.constant && size()>1) {
		if (size()%2==0) {
			opts.angle += pi/2.0;
			if (opts.angle>2.0*pi) opts.angle -= 2.0*pi;
		}
		else if (size()==3) {
			double angle = opts.angle - pi/2.0;
			if (angle<0) angle += 2.0*pi;
			double dx_n = sqrt(pow(opts.epsi_x*cos(angle),2.0)+pow(opts.epsi_y*sin(angle),2.0));
			double dx_t = sqrt(pow(opts.epsi_x*sin(angle),2.0)+pow(opts.epsi_y*cos(angle),2.0));
			opts.angle = -atan(((f_xy[1]).second-(f_xy[0]).second)*dx_t/((f_xy[2]).second-(f_xy[0]).second)/dx_n);
		}
		else {
			double angle = opts.angle - pi/2.0;
			if (angle<0) angle += 2.0*pi;
			uint j = size()-1;
			double dx_n = sqrt(pow(opts.epsi_x*cos(angle),2.0)+pow(opts.epsi_y*sin(angle),2.0));
			double dx_t = sqrt(pow(opts.epsi_x*sin(angle),2.0)+pow(opts.epsi_y*cos(angle),2.0));
			double angle_n = -atan(((f_xy[j-1]).second-(f_xy[j-3]).second)*dx_t/((f_xy[j]).second-(f_xy[j-1]).second)/dx_n);
			opts.angle = angle + angle_n;
			if (opts.angle>2.0*pi) opts.angle -= 2.0*pi;
			else if (opts.angle<0)opts.angle += 2.0*pi;
		}
	}
}


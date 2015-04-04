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
	2. fns of Point2d
	3. Stepper
	
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

// operator +
Point2d operator+(const Point2d& p1, const Point2d& p2) {
	Point2d p(p1.X+p2.X,p1.Y+p2.Y);
	return p;
}

// operator -
Point2d operator-(const Point2d& p1, const Point2d& p2) {
	Point2d p(p1.X-p2.X,p1.Y-p2.Y);
	return p;
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. functions of Point2d
		- distance
		- calcAngle
			calculates the angle between a line, defined by two points, and the x axis
		- find nth closest
	n.b. all static
-------------------------------------------------------------------------------------------------------------------------*/

// distance
static double distance(const Point2d& p1, const Point2d& p2) {
	return sqrt(pow(p1.X-p2.X,2.0)+pow(p1.Y-p2.Y,2.0));
}


// calcAngle
static double calcAngle(const Point2d& p1, const Point2d& p2) {
	Point2d p3 = p2;
	p3.X += distance(p1,p2);
	return 2.0*asin(distance(p1,p3)/2.0/distance(p1,p2));
}

// find nth closest
uint find_nth_closest(const vector<FxyPair>& fxy, const double& f, const uint& n) {
	if (n>fxy.size()) {
		cerr << "find_nth_closest error: n(" << n << ") chosen larger than f_xy.size() = " << fxy.size() << endl;
		return 1;
	}
	else if (fxy.size()<3) {
		cerr << "find_nth_closest error: f_xy.size() = " << fxy.size() << " smaller than 3" << endl;
		return 1;
	}
	vector<uint> ints(fxy.size());
	for (uint j=0; j<ints.size(); j++)
		ints[j] = j;
	uint loc_smallest = 0;
	uint ints_loc_smallest = 0;
	for (uint j=0; j<n; j++) {
		loc_smallest = ints[0];
		ints_loc_smallest = 0;
		double test_smallest = absDiff((fxy[loc_smallest]).second,f);
		for (uint k=1; k<ints.size(); k++) {
			double testk = absDiff((fxy[ints[k]]).second,f);
			if (testk < test_smallest) {
				loc_smallest = ints[k];
				ints_loc_smallest = k;
				test_smallest = absDiff((fxy[loc_smallest]).second,f);
			}
		}
		ints.erase(ints.begin()+ints_loc_smallest);
	}
	return loc_smallest;
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. Stepper
		- constructors
		- setStart
		- x
		- y
		- stepAngle
		- local
		- keep
		
-------------------------------------------------------------------------------------------------------------------------*/

// constructor
Stepper::Stepper(const StepperOptions& sto, const double& X, const double& Y):\
			 opts(sto), f_xy_local(), f_xy_steps(), angle(sto.angle0){
	Point2d P(X,Y);
	FxyPair toAdd(P,0.0);
	f_xy_local.push_back(toAdd);
	f_xy_steps.push_back(toAdd);
	if (opts.stepType!=StepperOptions::straight && opts.closeness<MIN_NUMBER)
		cerr << "Stepper error: closeness must be larger than 0" << endl;
}

// constructor
Stepper::Stepper(const StepperOptions& sto, const Point2d& P):\
			 opts(sto), f_xy_local(), f_xy_steps(), angle(sto.angle0){
	FxyPair toAdd(P,0.0);
	f_xy_local.push_back(toAdd);
	f_xy_steps.push_back(toAdd);
	if (opts.stepType!=StepperOptions::straight && abs(opts.closeness)<MIN_NUMBER)
		cerr << "Stepper error: closeness must be larger than 0" << endl;
}

// constructor
Stepper::Stepper(const StepperOptions& sto):\
			 opts(sto), f_xy_local(), f_xy_steps(), angle(sto.angle0){
	if (opts.stepType!=StepperOptions::straight && abs(opts.closeness)<MIN_NUMBER)
		cerr << "Stepper error: closeness must be larger than 0" << endl;
}

// set Start
void Stepper::setStart(const double& X, const double& Y) {
	f_xy_local.clear();
	f_xy_steps.clear();
	Point2d P(X,Y);
	FxyPair toAdd(P,0.0);
	f_xy_local.push_back(toAdd);
	f_xy_steps.push_back(toAdd);
}

// set Start
void Stepper::setStart(const Point2d& P) {
	f_xy_local.clear();
	f_xy_steps.clear();
	FxyPair toAdd(P,0.0);
	f_xy_local.push_back(toAdd);
	f_xy_steps.push_back(toAdd);
}

// x()
double Stepper::x() const {
	return ((f_xy_local.back()).first).X;
}

// y()
double Stepper::y() const {
	return ((f_xy_local.back()).first).Y;
}

// stepAngle()
double Stepper::stepAngle() const {
	return angle;
}

// local()
uint Stepper::local() const {
	return f_xy_local.size();
}

// keep()
bool Stepper::keep() const {
	return (absDiff((f_xy_local.back()).second,(f_xy_steps[0]).second)<opts.closeness \
				|| opts.stepType==StepperOptions::straight);
}

// point()
Point2d Stepper::point() const{
	return (f_xy_local.back()).first;
}

// steps()
uint Stepper::steps() const {
	return f_xy_steps.size()-1;
}

// stepVariables
void Stepper::step() {
	if (steps()<0) {
		cerr << "Stepper error: cannot step before giving initial step" << endl;
		return;
	}
	if (abs((f_xy_local.back()).second)<MIN_NUMBER && opts.stepType!=StepperOptions::straight) {
		cerr << "Stepper error: cannot step before giving result of previous step" << endl;
		return;
	}
	double x_old = ((f_xy_steps.back()).first).X;
	double y_old = ((f_xy_steps.back()).first).Y;
	double x_new = x_old + opts.epsi_x*cos(angle);
	double y_new = y_old + opts.epsi_y*sin(angle);
	Point2d P(x_new,y_new);
	FxyPair toAdd(P,0.0);
	f_xy_local.push_back(toAdd);
}

// addResult
void Stepper::addResult(const double& f) {
	(f_xy_local.back()).second = f;
	if (local()==1 || opts.stepType==StepperOptions::straight) {
		(f_xy_steps.back()).second = f;
	}
	else if (opts.stepType==StepperOptions::constSimple) {
		if (keep()) {
			f_xy_steps.push_back(f_xy_local.back());
			f_xy_local.clear();
			f_xy_local.push_back(f_xy_steps[steps()-1]);
			f_xy_local.push_back(f_xy_steps[steps()]);
			if (opts.directed==StepperOptions::local)
				opts.angle0 = angle;
			// angle unchanged if step was successful
		}
		else if (steps()==0) {
			if (local()==2) {
				angle += pi/2.0;	
				// no need to mod angle on first step			
			}
			else if (local()>2) {
				Point2d p_step = (f_xy_steps.back()).first;
				double f_step = (f_xy_steps.back()).second;
				uint l1 = find_nth_closest(f_xy_local,f_step,2);
				uint l2 = find_nth_closest(f_xy_local,f_step,3);
				Point2d p1 = (f_xy_local[l1]).first, p2 = (f_xy_local[l2]).first;
				double f1 = (f_xy_local[l1]).second, f2 = (f_xy_local[l2]).second;
				double norm = absDiff(f1,f_step) + absDiff(f2,f_step);
				angle = (absDiff(f1,f_step)/norm)*calcAngle(p_step,p1) + (absDiff(f2,f_step)/norm)*calcAngle(p_step,p2);
				if (opts.directed!=StepperOptions::undirected) {
					if (angle>pi)
						angle -= 2.0*pi;
					else if (angle<-pi)
						angle +=2.0*pi;
					if (angle<(opts.angle0-pi/2.0))
						angle += pi;
					else if (angle>(opts.angle0+pi/2.0))
						angle -= pi;
				}
			}
			else {
				cerr << "Stepper error: addResult should not have reached this point, 1" << endl;
			}
		}
		else if (local()<2) {
			cerr << "Stepper error: addResult shoudl not have reached this point, 2" << endl;
		}
		else if (local()==3) {
			angle += pi/2.0;
			if (opts.directed!=StepperOptions::undirected) {
				if (angle>pi)
					angle -= 2.0*pi;
				else if (angle<-pi)
					angle += 2.0*pi;
				if (angle<(opts.angle0-pi/2.0-MIN_NUMBER))
					angle += pi;
				else if (angle>(opts.angle0+pi/2.0+MIN_NUMBER))
					angle -= pi;
			}
		}
		else if (local()>3) {
			Point2d p_step = (f_xy_steps.back()).first;
			double f_step = (f_xy_steps.back()).second;
			uint l1 = find_nth_closest(f_xy_local,f_step,2);
			uint l2 = find_nth_closest(f_xy_local,f_step,3);
			Point2d p1 = (f_xy_local[l1]).first, p2 = (f_xy_local[l2]).first;
			double f1 = (f_xy_local[l1]).second, f2 = (f_xy_local[l2]).second;
			double norm = absDiff(f1,f_step) + absDiff(f2,f_step);
			angle = (absDiff(f1,f_step)/norm)*calcAngle(p_step,p1) + (absDiff(f2,f_step)/norm)*calcAngle(p_step,p2);
			if (opts.directed!=StepperOptions::undirected) {
				if (angle>pi)
					angle -= 2.0*pi;
				else if (angle<-pi)
					angle +=2.0*pi;
				if (angle<(opts.angle0-pi/2.0-MIN_NUMBER))
					angle += pi;
				else if (angle>(opts.angle0+pi/2.0+MIN_NUMBER))
					angle -= pi;
			}
		}
	}
	else {
		cerr << "Stepper error: addResult shoudl not have reached this point, 3" << endl;
	}
}

// addResult
void Stepper::addResult(const double& f, const double& e, const double& n) {
	if (opts.stepType!=StepperOptions::lagrange) {
		addResult(f);
	}
	else {
		(f_xy_local.back()).second = f;
		if (steps()>0) {
			cerr << "Stepper error: lagrange stepper not written yet" << endl;
		}
	}
}

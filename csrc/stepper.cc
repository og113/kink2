/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for the Stepper class and dependencies
 -------------------------------------------------------------------------------------------------------------------------*/

#include <ctime> // for time
#include <cstdlib> //for rand, srand
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
static uint find_nth_closest(const vector<FxyPair>& fxy, const double& f, const uint& n) {
	if (n>fxy.size()) {
		//cerr << "find_nth_closest error: n(" << n << ") chosen larger than f_xy.size() = " << fxy.size() << endl;
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
	double test = absDiff((f_xy_local.back()).second,(f_xy_steps[0]).second);
	return ( (test<opts.closeness && ((steps()==0 && local()==1) || (steps()>0 && local()==2)))\
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
		double test = absDiff((f_xy_local.back()).second,(f_xy_steps[0]).second);
		if (test<opts.closeness && local()>3) {
			f_xy_steps.push_back(f_xy_local.back());
			f_xy_local.clear();
			f_xy_local.push_back(f_xy_steps[steps()-1]);
			f_xy_local.push_back(f_xy_steps[steps()]);
			if (opts.directed==StepperOptions::local)
				opts.angle0 = angle;
			angle += pi/2.0;
		}
		else if (local()==2) {
			angle += pi/2.0;		
		}
		else if (local()==3 && steps()==0) {
			double tempAngle = angle - pi/2.0;
			double dx_n = sqrt(pow(opts.epsi_x*cos(tempAngle),2.0)+pow(opts.epsi_y*sin(tempAngle),2.0));
			double dx_t = sqrt(pow(opts.epsi_x*sin(tempAngle),2.0)+pow(opts.epsi_y*cos(tempAngle),2.0));
			angle = -atan(((f_xy_local[1]).second-(f_xy_steps[0]).second)*dx_t/((f_xy_local[2]).second-(f_xy_steps[0]).second)/dx_n);
		}
		else if (local()==3) {
			double tempAngle = angle - pi/2.0;
			double dx_n = sqrt(pow(opts.epsi_x*cos(tempAngle),2.0)+pow(opts.epsi_y*sin(tempAngle),2.0));
			double dx_t = sqrt(pow(opts.epsi_x*sin(tempAngle),2.0)+pow(opts.epsi_y*cos(tempAngle),2.0));
			double angle_n = -atan(((f_xy_steps[steps()]).second-(f_xy_steps[steps()-1]).second)*dx_t \
									/((f_xy_local[local()-1]).second-(f_xy_steps[steps()]).second)/dx_n);
			angle = tempAngle + angle_n;
		}
		else if (local()>3) {
			Point2d p_step = (f_xy_steps.back()).first;
			p_step.X /= opts.epsi_x;
			p_step.Y /= opts.epsi_y;
			double f_0 = (f_xy_steps[0]).second;
			vector <FxyPair> f_low, f_high;
			for (uint k=0; k<local(); k++) {
				if ((f_xy_local[k]).second<(f_0-MIN_NUMBER*1.0e2))
					f_low.push_back(f_xy_local[k]);
				else if ((f_xy_local[k]).second>(f_0+MIN_NUMBER*1.0e2))
					f_high.push_back(f_xy_local[k]);
			}
			if ((f_low.size()==0 || f_high.size()==0) || \
					 absDiff((f_xy_local[local()-1]).second,(f_xy_local[local()-2]).second)<MIN_NUMBER) {
				srand(time(NULL));
				angle += randDouble(-pi,pi);
			}
			else {
				uint ll = find_nth_closest(f_low,f_0,1), lh = find_nth_closest(f_high,f_0,1);				
				Point2d pl = (f_low[ll]).first, ph = (f_high[lh]).first;
				pl.X /= opts.epsi_x;
				pl.Y /= opts.epsi_y;
				ph.X /= opts.epsi_x;
				ph.Y /= opts.epsi_y;
				double fl = (f_low[ll]).second, fh = (f_high[lh]).second;
				double norm = absDiff(fl,f_0) + absDiff(fh,f_0);
				double anglel = calcAngle(p_step,pl), angleh = calcAngle(p_step,ph);
				angle = (absDiff(fl,f_0)/norm)*anglel + (absDiff(fh,f_0)/norm)*angleh;
			}
		}
		
		else
			cerr << "Stepper error: addResult should not have reached this point, 1" << endl;
	}
	else
		cerr << "Stepper error: addResult should not have reached this point, 2" << endl;
		
	if (opts.directed!=StepperOptions::undirected) {	
		if (angle<(opts.angle0-pi/2.0-MIN_NUMBER*1.0e2))
			angle += pi;
		else if (angle>(opts.angle0+pi/2.0+MIN_NUMBER*1.0e2))
			angle -= pi;
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

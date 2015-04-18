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
	2. static functions
	3. Stepper
	
n.b. stepper defined in 2d.
n.b. on 10/04/15 constPlane was the only reliable stepper. considerably better than taylor.
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. Point2d
		- copy
		- copy constructor
		- operator=
		- operator()
		- operator<<
		- operator+
		- operator-
		- operator==
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

// operator==
bool operator==(const Point2d& lhs, const Point2d& rhs) {
	return (abs(lhs.X-rhs.X)<MIN_NUMBER && abs(lhs.Y-rhs.Y)<MIN_NUMBER);
}

/*-------------------------------------------------------------------------------------------------------------------------
	2. static functions
		- staticModder
		- distance
		- calcAngle
			calculates the angle between a line, defined by two points, and the x axis
		- find nth closest
		- isBetween
		- isClockwise
-------------------------------------------------------------------------------------------------------------------------*/

// static modder
static void staticModder(double& a, const StepperOptions& o) {
	a = mod(a,-pi-MIN_NUMBER*1.0e4,pi+MIN_NUMBER*1.0e4);
	if (o.directed!=StepperOptions::undirected) {
		a = mod(a,o.angle0-pi/2.0-MIN_NUMBER*1.0e4,o.angle0+pi/2.0+MIN_NUMBER*1.0e4);
	}
	a = mod(a,-pi-MIN_NUMBER*1.0e4,pi+MIN_NUMBER*1.0e4);
}

// distance
/*static double Distance(const Point2d& p1, const Point2d& p2) {
	return sqrt(pow(p1.X-p2.X,2.0)+pow(p1.Y-p2.Y,2.0));
}

// calcAngle
static double calcAngle(const Point2d& p1, const Point2d& p2) {
	if (Distance(p1,p2)<MIN_NUMBER)
		return 0.0;
	Point2d p3(p2);
	p3.X += Distance(p1,p2);
	return 2.0*asin(Distance(p1,p3)/2.0/Distance(p1,p2));
}*/


// calcAngle
static double calcAngle(const Point2d& p1, const Point2d& p2) {
	if (p1==p2)
		return 0.0;
	double toReturn = atan((p2.Y-p1.Y)/(p2.X-p1.X));
	return toReturn;
}

// calcEpsiAngle
static double calcEpsiAngle(const Point2d& p1, const Point2d& p2, const StepperOptions& o) {
	if (p1==p2)
		return 0.0;
	Point2d P1 = p1, P2 = p2;
	P1.X /= o.epsi_x;
	P1.Y /= o.epsi_y;
	P2.X /= o.epsi_x;
	P2.Y /= o.epsi_y;
	double toReturn = calcAngle(P1,P2);
	staticModder(toReturn,o);
	return toReturn;
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

// isBetween
static bool isBetween(const double& v, const double& a, const double& b) {
	double temph = (a<b? b: a);
	double h = mod(temph,-pi-MIN_NUMBER*1.0e4,pi+MIN_NUMBER*1.0e4);
	double templ = (a<b? a: b);
	double l = mod(templ,-pi-MIN_NUMBER*1.0e4,pi+MIN_NUMBER*1.0e4);
	double difference = h-l;
	double V = mod(v,-pi-MIN_NUMBER*1.0e4,pi+MIN_NUMBER*1.0e4);
	if (difference<=pi && V<h && V>l)
		return true;
	else if (difference<=pi && V>h && V<l)
		return true;
	else
		return false;
}

// isClockwise, is v within pi in a clockwise direction from a?
static bool isClockwise(const double& v, const double& a) {
	double A = mod(a,-pi-MIN_NUMBER*1.0e4,pi+MIN_NUMBER*1.0e4);
	double V = mod(v,-pi-MIN_NUMBER*1.0e4,pi+MIN_NUMBER*1.0e4);
	double difference = abs(A-V);
	if (difference<=pi && V<A)
		return true;
	else if (difference>=pi && V>A)
		return true;
	else
		return false;
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. Stepper
		- constructors
		- setStart
		- x
		- y
		- result
		- closeness
		- stepAngle
		- local
		- keep
		
-------------------------------------------------------------------------------------------------------------------------*/

// constructor
Stepper::Stepper(const StepperOptions& sto, const double& X, const double& Y, const double& f0):\
 		opts(sto), f_xy_local(), f_xy_steps(), angle(sto.angle0){
	Point2d P(X,Y);
	FxyPair toAdd(P,f0);
	f_xy_local.push_back(toAdd);
	f_xy_steps.push_back(toAdd);
	if (opts.stepType!=StepperOptions::straight && opts.closeness<MIN_NUMBER)
		cerr << "Stepper error: closeness must be larger than 0" << endl;
}

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
Stepper::Stepper(const StepperOptions& sto, const Point2d& P, const double& f0):\
			 opts(sto), f_xy_local(), f_xy_steps(), angle(sto.angle0){
	FxyPair toAdd(P,f0);
	f_xy_local.push_back(toAdd);
	f_xy_steps.push_back(toAdd);
	if (opts.stepType!=StepperOptions::straight && abs(opts.closeness)<MIN_NUMBER)
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

// result()
double Stepper::result() const {
	return (f_xy_local.back()).second;
}

// closeness()
double Stepper::closeness() const {
	return opts.closeness;
}

// stepAngle()
double Stepper::stepAngle() const {
	return angle;
}

// local()
uint Stepper::local() const {
	return f_xy_local.size();
}

// keep() - n.b. only works after addResult
bool Stepper::keep() const {
	double test = absDiff((f_xy_local.back()).second,(f_xy_steps[0]).second);
	return ( (test<opts.closeness && (steps()==0 || local()==2)) || opts.stepType==StepperOptions::straight);
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
	Point2d P0 = (f_xy_steps.back()).first;
	double x_old = P0.X;
	double y_old = P0.Y;
	double x_new = x_old + opts.epsi_x*cos(angle);
	double y_new = y_old + opts.epsi_y*sin(angle);
	Point2d P(x_new,y_new);
	if (opts.stepType!=StepperOptions::straight) {
		vector <FxyPair> f_xy_local_minus_loc = f_xy_local;
		steps()==0? f_xy_local_minus_loc.erase(f_xy_local_minus_loc.begin()):\
								 f_xy_local_minus_loc.erase(f_xy_local_minus_loc.begin()+1);
		for (uint k=0; k<(local()-1); k++) {
			if (P==(f_xy_local_minus_loc[k]).first) {
				srand(time(NULL));
				angle += randDouble(-pi,pi)/4.0;
				staticModder(angle,opts);
				x_new = x_old + opts.epsi_x*cos(angle);
				y_new = y_old + opts.epsi_y*sin(angle);
				P.X = x_new;
				P.Y = y_new;
			}
		}
	}
	FxyPair toAdd(P,0.0);
	f_xy_local.push_back(toAdd);
}

// addResult
void Stepper::addResult(const double& f) {
	(f_xy_local.back()).second = f;
	if (local()==1) {
		(f_xy_steps.back()).second = f;
	}
	else if (opts.stepType==StepperOptions::straight) {
		f_xy_steps.push_back(f_xy_local.back());
		f_xy_local.clear();
		f_xy_local.push_back(f_xy_steps[steps()-1]);
		f_xy_local.push_back(f_xy_steps[steps()]);
	}
	else if (opts.stepType==StepperOptions::constPlane) {
		double test = absDiff(f,(f_xy_steps[0]).second);
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
			Point2d p0 = (f_xy_steps[0]).first, p1 = (f_xy_local[1]).first, p2 = (f_xy_local[2]).first;
			double f0 = (f_xy_steps[0]).second, f1 = (f_xy_local[1]).second, f2 = (f_xy_local[2]).second;
			double numerator = (f2-f0)*(p1.Y-p0.Y) - (f1-f0)*(p2.Y-p0.Y);
			double denominator = (f2-f0)*(p1.X-p0.X) - (f1-f0)*(p2.X-p0.X);
			angle = atan(numerator/denominator);
			if (!isBetween(angle,opts.angle0,opts.angle0+pi/2.0))
				angle +=pi;
		}
		else if (local()==3) {
			double tempAngle = angle - pi/2.0;
			Point2d p0 = (f_xy_steps[steps()]).first, p1 = (f_xy_steps[steps()-1]).first, p2 = (f_xy_local[local()-1]).first;
			double f0 = (f_xy_steps[0]).second, f1 = (f_xy_steps[steps()-1]).second, f2 = (f_xy_local[local()-1]).second;
			double numerator = (f2-f0)*(p1.Y-p0.Y) - (f1-f0)*(p2.Y-p0.Y);
			double denominator = (f2-f0)*(p1.X-p0.X) - (f1-f0)*(p2.X-p0.X);
			angle = atan(numerator/denominator);
			if (!isBetween(angle,tempAngle,tempAngle+pi/2.0))
				angle +=pi;
		}
		else if (local()>3) {
			Point2d p0 = (f_xy_steps.back()).first;
			double f0 = (f_xy_steps[0]).second;
			vector <FxyPair> f_low, f_high, f_xy_local_minus_loc = f_xy_local;
			steps()==0? f_xy_local_minus_loc.erase(f_xy_local_minus_loc.begin()):\
								 f_xy_local_minus_loc.erase(f_xy_local_minus_loc.begin()+1);
			for (uint k=0; k<(local()-1); k++) {
				if ((f_xy_local_minus_loc[k]).second<(f0-MIN_NUMBER*1.0e2))
					f_low.push_back(f_xy_local_minus_loc[k]);
				else if ((f_xy_local_minus_loc[k]).second>(f0+MIN_NUMBER*1.0e2))
					f_high.push_back(f_xy_local_minus_loc[k]);
			}
			if (f_low.size()==0) {
				FxyPair f_high_min = f_high[find_nth_closest(f_high,f0,1)];
				FxyPair f_high_second = f_high[find_nth_closest(f_high,f0,2)];
				double angle_min = calcEpsiAngle(p0,f_high_min.first,opts);
				double angle_second = calcEpsiAngle(p0,f_high_second.first,opts);
				staticModder(angle_min,opts);
				staticModder(angle_second,opts);
				int sign = (isClockwise(angle_min,angle_second)? -1.0: 1.0);
				srand(time(NULL));
				angle = angle_min + sign*randDouble(0.0,pi)/2.0;
				if (f_high.size()>5)
					angle += sign*randDouble(0.0,pi)/2.0;
			}
			else if (f_high.size()==0) {
				FxyPair f_low_min = f_low[find_nth_closest(f_low,f0,1)];
				FxyPair f_low_second = f_low[find_nth_closest(f_low,f0,2)];
				double angle_min = calcEpsiAngle(p0,f_low_min.first,opts);
				double angle_second = calcEpsiAngle(p0,f_low_second.first,opts);
				staticModder(angle_min,opts);
				staticModder(angle_second,opts);
				int sign = (isClockwise(angle_min,angle_second)? -1.0: 1.0);
				srand(time(NULL));
				angle = angle_min + sign*randDouble(0.0,pi)/2.0;
				if (f_high.size()>5)
					angle += sign*randDouble(0.0,pi)/2.0;
			}
			else {
				uint ll = find_nth_closest(f_low,f0,1), lh = find_nth_closest(f_high,f0,1);			
				Point2d pl = (f_low[ll]).first, ph = (f_high[lh]).first;
				double fl = (f_low[ll]).second, fh = (f_high[lh]).second;
				double numerator = (fh-f0)*(pl.Y-p0.Y) - (fl-f0)*(ph.Y-p0.Y);
				double denominator = (fh-f0)*(pl.X-p0.X) - (fl-f0)*(ph.X-p0.X);
				angle = atan(numerator/denominator);
				double anglel = calcEpsiAngle(pl,p0,opts);
				double angleh = calcEpsiAngle(ph,p0,opts);
				if (!isBetween(angle,anglel,angleh))
					angle +=pi;
			}
		}
	}
	/*else if (opts.stepType==StepperOptions::constTaylor) {
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
			if (!isBetween(angle,opts.angle0,opts.angle0+pi/2.0))
				angle +=pi;
		}
		else if (local()==3) {
			double tempAngle = angle - pi/2.0;
			double dx_n = sqrt(pow(opts.epsi_x*cos(tempAngle),2.0)+pow(opts.epsi_y*sin(tempAngle),2.0));
			double dx_t = sqrt(pow(opts.epsi_x*sin(tempAngle),2.0)+pow(opts.epsi_y*cos(tempAngle),2.0));
			double angle_n = -atan(((f_xy_steps[steps()]).second-(f_xy_steps[steps()-1]).second)*dx_t \
									/((f_xy_local[local()-1]).second-(f_xy_steps[steps()]).second)/dx_n);
			angle = tempAngle + angle_n;
			if (!isBetween(angle,tempAngle,tempAngle+pi/2.0))
				angle +=pi;
		}
		else if (local()>3) {
			Point2d p_step = (f_xy_steps.back()).first;
			double f_0 = (f_xy_steps[0]).second;
			vector <FxyPair> f_low, f_high;
			for (uint k=0; k<local(); k++) {
				if ((f_xy_local[k]).second<(f_0-MIN_NUMBER*1.0e2))
					f_low.push_back(f_xy_local[k]);
				else if ((f_xy_local[k]).second>(f_0+MIN_NUMBER*1.0e2))
					f_high.push_back(f_xy_local[k]);
			}
			if (f_low.size()==0 || f_high.size()==0) {
				srand(time(NULL));
				angle += randDouble(-pi,pi);
			}
			else {
				uint ll = find_nth_closest(f_low,f_0,1), lh = find_nth_closest(f_high,f_0,1);				
				Point2d pl = (f_low[ll]).first, ph = (f_high[lh]).first;
				double fl = (f_low[ll]).second, fh = (f_high[lh]).second;
				double norm = absDiff(fl,f_0) + absDiff(fh,f_0);
				double anglel = calcEpsiAngle(p_step,pl,opts), angleh = calcEpsiAngle(p_step,ph,opts);
				angle = (absDiff(fl,f_0)/norm)*anglel + (absDiff(fh,f_0)/norm)*angleh;
				if (!isBetween(angle,anglel,angleh))
					angle +=pi; 
			}
		}
	}*/
	else {
		cerr << "Stepper error: addResult should not have reached this point, 2" << endl;
	}
	
	staticModder(angle,opts);
}

/*-------------------------------------------------------------------------------------------------------------------------
	definitions for functions using eigen library functions
-------------------------------------------------------------------------------------------------------------------------*/

#include <Eigen/Dense>
#include <string>
#include <iostream>
#include "eigen_extras.h"
#include "simple.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1 - printErrorInformation
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
 1. printErrorInformation
-------------------------------------------------------------------------------------------------------------------------*/

// printErrorInformation
void printErrorInformation(const vec& v, const string& name) {
	uint minCoeff1 = 0, maxCoeff1 = 0, largeCounter = 0;
	number max = v.maxCoeff(&maxCoeff1);
	number min = v.minCoeff(&minCoeff1);
	number absMax = (-min>max? -min: max);
	for (uint j=0; j<v.size(); j++) {
		if (abs(v(j))>absMax/2.0)
			largeCounter++;
	}
	
	cout << name << " information:" << endl;
	cout << name <<".size()    :         " << v.size() << endl;
	cout << name <<".norm()    :         " << v.norm() << endl;
	cout << name <<".mean()    :         " << v.mean() << endl;
	cout << name <<".minCoeff():         " << min  << endl;
	cout << name <<".maxCoeff():         " << max  << endl;
	cout << name <<" minCoeff  :         " << maxCoeff1 << endl;
	cout << name <<" maxCoeff  :         " << minCoeff1 << endl;
	cout << name <<" absmax/2 counter :  " << largeCounter << endl;
}

// printErrorInformation
void printErrorInformation(const vec& v, const string& name, const uint zm) {
	printErrorInformation(v,name);
	cout << "(" << name <<".tail(zm)).norm()/sqrt(zm): " << (v.tail(zm)).norm()/sqrt(zm) << endl;
	cout << "(" << name <<".head(v.size()-zm)).norm()/sqrt(v.size()-zm):  " << (v.head(v.size()-zm)).norm()/sqrt(v.size()-zm) << endl;

}

// printErrorInformation
void printErrorInformation(const mat& m, const string& name) {
	uint minCoeff1 = 0, maxCoeff1 = 0, minCoeff2 = 0, maxCoeff2 = 0, largeCounter = 0;
	number max = m.maxCoeff(&maxCoeff1,&maxCoeff2);
	number min = m.minCoeff(&maxCoeff1,&maxCoeff2);
	number absMax = (-min>max? -min: max);
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			if (abs(m(j,k))>absMax/2.0)
				largeCounter++;
		}
	}

	cout << name << " information:" << endl;
	cout << name <<".rows()    :         " << m.rows() << endl;
	cout << name <<".cols()    :         " << m.cols() << endl;
	cout << name <<".norm()    :         " << m.norm() << endl;
	cout << name <<".minCoeff():         " << min << endl;
	cout << name <<".maxCoeff():         " << max << endl;
	cout << name <<" minCoeff  :         (" << maxCoeff1 << "," << maxCoeff2 << ")" <<endl;
	cout << name <<" maxCoeff  :         (" << minCoeff1 << "," << minCoeff2 << ")" <<endl;
	cout << name <<" absmax/2 counter :  " << largeCounter << endl;
	cout << name <<".trace():            " << m.trace() << endl;
	cout << name <<".determinant():      " << m.determinant() << endl;
	cout << name <<".mean():             " << m.mean() << endl;
	cout << name <<".sum():              " << m.sum() << endl;
	cout << name <<".prod():             " << m.prod() << endl;

}

// printErrorInformation
void printErrorInformation(const spMat& m, const string& name) {

	cout << name << " information:" << endl;
	cout << name <<".rows()    :         " << m.rows() << endl;
	cout << name <<".cols()    :         " << m.cols() << endl;
	cout << name <<".nonZeros():         " << m.nonZeros() << endl;
	cout << name <<".norm()    :         " << m.norm() << endl;
	cout << name <<".blueNorm():         " << m.blueNorm() << endl;
	cout << name <<".sum()     :         " << m.sum() << endl;
}

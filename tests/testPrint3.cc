/*
	test of binary printing
*/

#include<iostream>
#include<fstream>
#include<string>
#include<complex>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include "folder.h"
#include "parameters.h"
#include "print3.h"
#include "simple.h"

using namespace std;

int main() {

cout << "testPrint3:" << endl << endl;

Filename file = "data/pot_3/field/fMain_N_130_Na_300_Nb_80_Nc_2_LoR_5_Tb_0.79_theta_0_reg_0.data";
Filename file2 = "data/temp/fTemp.data";

vec phi, phi2;

cout << "phi.size() = " << phi.size() << endl;
cout << "phi.norm() = " << phi.norm() << endl << endl;

loadVectorBinary(file,phi);
cout << "phi.size() = " << phi.size() << endl;
cout << "phi.norm() = " << phi.norm() << endl << endl;

cout << "file.exists() = " << file.exists() << endl;
cout << "file2.exists() = " << file2.exists() << endl;
saveVectorBinary(file2,phi);
loadVectorBinary(file2,phi2);
cout << "phi2.size() = " << phi2.size() << endl;
cout << "phi2.norm() = " << phi2.norm() << endl << endl;

cout << "file2.exists() = " << file2.exists() << endl;

return 0;
}

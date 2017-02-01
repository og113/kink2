/*
	test of the function folderToParameters in nr.cc
*/

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "simple.h"
#include "folder.h"
#include "nr.h"

using namespace std;

int main() {
cout << "test FilenameToParameters: " << endl;

Filename file1 = "data/pot_3/field/fMain_N_130_Na_300_Nb_80_Nc_2_LoR_5_DE_0_Tb_0.73_Theta_0_Reg_0.data";
Filename file2 = "data/pot_3/omega/freqsF_N_130_LoR_5_DE_0_Reg_0.data";

Parameters p1 = filenameToParameters(file1);
Parameters p2 = filenameToParameters(file2);

cout << "file1: " << file1 << endl;
cout << "p1   : " << endl << p1 << endl << endl;

cout << "file2: " << file2 << endl;
cout << "p2   : " << endl << p2 << endl << endl;

return 0;
}

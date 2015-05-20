/*
	program to test the classes and functions defined in folder.h and folder.cc
*/

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include "simple.h"
#include "folder.h"

using namespace std;

int main() {
cout << "folderTest.cc:" << endl;
Filename f;
f.Directory = "data";
f.Timenumber = "";
f.ID = "piE";
vector<StringPair> exs(1);
exs[0] = StringPair("N","100");
f.Extras = exs;
f.Suffix = ".dat";
cout << f << " " << f() << endl;
Filename g(f);
g.ID = "minusDSE";
cout << g << endl;
Filename h(f());
cout << h << endl;
h.Timenumber = "0";
cout << h << endl;
Filename p("data/stable/455data1_L_5_Tb_0.8_R_4.2.dat");
cout << " Extras: " << endl;
for (unsigned int loop=0; loop<(p.Extras).size(); loop++) {
	cout << (p.Extras)[loop].first << " " << (p.Extras)[loop].second << endl;
}
cout << p << endl;
FilenameAttributes j;
Filename k = (string)"data/30thingy_L_5_Tb_0.8.dat";
Filename l = (string)"data/blahblah_L_5.txt";
Folder FolderL(l);
cout << "FolderL =  " << FolderL << endl << endl;
k.Timenumber = "40";
cout << "k = " << k << endl;
j.ID = "testError";
Folder folder1(j); // test may die here if not run from the folder tests
cout << "folder1 = " << folder1 << endl;
cout << "folder2 = ";
Folder folder2(k);
cout << folder2 << endl;
folder2.set(k);
//folder2.update();
if (folder2.isPresent(p)) cout << "Present" << endl;
cout << folder1.size() << " " << folder2.size() << endl;
cout << "Folders: " << endl << folder1 << endl << folder2 << endl;
cout << j << endl << k << endl;

removeUnshared(folder1,folder2);
cout << folder1 << endl;

string str = f;
cout << str << endl;

FilenameAttributes fa_low;
fa_low.Directory = "data";
fa_low.Suffix = ".data";
fa_low.ID = "mainp";
fa_low.Timenumber = "0";
FilenameAttributes fa_high(fa_low);
fa_high.Timenumber = "999999999999";
FilenameComparator fc(fa_low,fa_high);
Folder F(fc);
cout << F << endl;
cout << "compared with typing: ls data/*mainp*.data" << endl;

return 0;
}

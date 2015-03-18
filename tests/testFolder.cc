/*
	program to test the classes and functions defined in folder.h and folder.cc
*/

#include <iostream>
#include <string>
#include <vector>
#include <utility>
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
Filename p("data/data1.dat");
cout << p << " " << p.Directory << " " << p.Timenumber << " " << p.ID << " " << p.Suffix << endl;
FilenameAttributes j;
Filename k = (string)"data/30thingy.dat";
k.Timenumber = "40";
Folder folder1(j), folder2(k); // test may die here if not run from the folder tests
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

return 0;
}

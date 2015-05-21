/*----------------------------------------------------------------------------------------------------------------------------
	wrapper
		program to give simple mpi wrapper for main
----------------------------------------------------------------------------------------------------------------------------*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <mpi.h>
#include "folder.h"
#include "main_fn.h"
#include "parameters.h"
#include "simple.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - defining required nodes
		2 - copying files
		3 - finding recent files
		4 - running main
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/*----------------------------------------------------------------------------------------------------------------------------
	1. defining required nodes
		- initializing mpi
----------------------------------------------------------------------------------------------------------------------------*/

int nodes_req = 8;

int nodes, rank;
int returnValue = 0;

MPI::Init(argc, argv);
MPI_Comm_size(MPI_COMM_WORLD, &nodes);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
if (nodes==nodes_req) {
	if (rank==0) {
		cout << "running " << nodes << " nodes" << endl;
	}
}
else {
	if (rank==0) {
		cerr << "have " << nodes << " nodes available" << endl;
		cerr << "require " << nodes_req << " nodes to run" << endl;
	}
	MPI::Finalize();
	return 1;
}

/*----------------------------------------------------------------------------------------------------------------------------
	2. copying files
		- copying optionsM
		- copying negEig, omega etc
----------------------------------------------------------------------------------------------------------------------------*/

bool copyFiles = false;

if (copyFiles && rank==0) {
	for (int k=0; k<nodes_req; k++) {
		string timenumber = "00"+numberToString<int>(k);
		// n.b. in the examples i have seen all quantities are declared before MPI::Init
		// should look this up
	
		// copying optionsM with changed timenumbers
		Options opts;
		opts.load("optionsM");
		Filename newOptsFile = (string)"data/"+timenumber+"optionsM";
		opts.minTimenumberLoad = timenumber;
		opts.maxTimenumberLoad = timenumber;
		opts.save(newOptsFile);

		// copying omega and negEig files
		Filename in = (string)"data/stable/eigVec_pot_3_L_5.dat";
		Filename out = in;
		out.Directory = "data";
		out.Timenumber = timenumber;
		copyFile(in,out);
		in.ID = "freqsExp";
		out.ID = in.ID;
		copyFile(in,out);
		in.ID = "freqs";
		out.ID = in.ID;
		copyFile(in,out);
		in.ID = "modes";
		out.ID = in.ID;
		copyFile(in,out);
		in.ID = "omegaM1";
		out.ID = in.ID;
		copyFile(in,out);
		in.ID = "omega0";
		out.ID = in.ID;
		copyFile(in,out);
		in.ID = "omega1";
		out.ID = in.ID;
		copyFile(in,out);
		in.ID = "omega2";
		out.ID = in.ID;
		copyFile(in,out);
	}
}

/*----------------------------------------------------------------------------------------------------------------------------
	3. finding recent files
		- getting last timenumbers
		- getting last loops
		
	n.b. there is lots of superfluous numTostr and vice versa due to not knowing how to MPI::Send and receive strings.
----------------------------------------------------------------------------------------------------------------------------*/

// getting last timenumbers
bool revertToDefault = false;
vector<string> timenumbers(nodes_req);
vector<string> loops(nodes_req);
string timenumber;
string loop;

FilenameAttributes fa_low;
fa_low.Directory = "data";
fa_low.Suffix = ".data";
fa_low.ID = "mainp";
fa_low.Timenumber = "0";
(fa_low.Extras).push_back(StringPair("step","1"));
FilenameAttributes fa_high(fa_low);
fa_high.Timenumber = "999999999999";
Folder F(fa_low,fa_high);

if (F.size()==0)
	revertToDefault = true;
else if (rank==0 && !revertToDefault) {
	for (int k=0; k<nodes_req; k++) {
		string maxTimenumber = "0";
		if (F.size()==0) {
			revertToDefault = true;
			break;
		}
		for (uint j=0; j<F.size(); j++) {
			if ((F[j]).Timenumber>maxTimenumber) {
				maxTimenumber = (F[j]).Timenumber;
			}
		}
		timenumbers[k] = maxTimenumber;
		unsigned long long l = stringToNumber<unsigned long long>(maxTimenumber);
		l--;
		fa_high.Timenumber = numberToString<unsigned long long>(l);
		F.set(fa_low,fa_high);
	}

	// getting last loops
	if (!revertToDefault) {
		for (int k=0; k<nodes_req; k++) {
			fa_low.Timenumber = timenumbers[k];
			fa_high.Timenumber = timenumbers[k];
			F.set(fa_low,fa_high);
			string maxLoop = "0";
			for (uint j=0; j<F.size(); j++) {
				for (uint l=0; l<(F[j].Extras).size(); l++) {
					if ((((F[j].Extras)[l]).first).compare("loop")==0) {
						if (((F[j].Extras)[l]).second>maxLoop)
							maxLoop = ((F[j].Extras)[l]).second;
					}
				}
			}
			loops[k] = maxLoop;
		}
	}
}

if (revertToDefault && rank==0) {
	for (int k=0; k<nodes_req; k++) {
		timenumbers[k] = "00"+numberToString<int>(k);
		loops[k] = "0";
	}
}

/*if (rank==0) {
	for (int k=1; k<nodes_req; k++) {
		timenumber = stringToNumber<unsigned long long>(timenumbers[k]); 
		loop = stringToNumber<unsigned long long>(loops[k]);
		MPI_Send(&timenumber,1, MPI_UNSIGNED_LONG_LONG, k, 0, MPI_COMM_WORLD);
		MPI_Send(&loop,1, MPI_UNSIGNED_LONG_LONG, k, 1, MPI_COMM_WORLD);
		//cout << "process " << 0 << " sent " << timenumber << " and " << loop << " to process " << k << endl;
	}
	timenumber = stringToNumber<unsigned long long>(timenumbers[0]);
	loop = stringToNumber<unsigned long long>(loops[0]);
}
else {
	MPI_Recv(&timenumber, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(&loop, 1, MPI_UNSIGNED_LONG_LONG, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//cout << "process " << rank << " received " << timenumber << " and " << loop << " from process " << 0 << endl;
}*/


if (rank==0) {
	for (int k=1; k<nodes_req; k++) {
		char temp[timenumbers[k].length()+1];
		timenumbers[k].copy(temp, timenumbers[k].length());
		temp[timenumbers[k].length()] = '\0';
		
		MPI::COMM_WORLD.Send(&temp, timenumbers[k].length(), MPI::CHAR, k, 0);
		
		char temp2[loops[k].length()+1];
		loops[k].copy(temp2, loops[k].length());
		temp2[loops[k].length()] = '\0';
		
		MPI::COMM_WORLD.Send(&temp2, loops[k].length(), MPI::CHAR, k, 1);
		
		//cout << "process " << 0 << " sent " << temp << " and " << temp2 << " to process " << k << endl;
	}
	timenumber = timenumbers[0];
	loop = loops[0];
}
else {
	MPI::Status status;
	MPI::COMM_WORLD.Probe(0, 0, status);
	int l = status.Get_count(MPI::CHAR);
	char *buf = new char[l];
	MPI::COMM_WORLD.Recv(buf, l, MPI::CHAR, 0, 0, status);
	timenumber = string (buf, l);
	delete [] buf;
	
	MPI::COMM_WORLD.Probe(0, 1, status);
	l = status.Get_count(MPI::CHAR);
	char *buf2 = new char[l];
	MPI::COMM_WORLD.Recv(buf2, l, MPI::CHAR, 0, 1, status);
	loop = string(buf2, l);
	delete [] buf2;
	
	//cout << "process " << rank << " recieved " << timenumber << " and " << loop << " from process 0" << endl;
}

/*----------------------------------------------------------------------------------------------------------------------------
	4. running main
		- defining argc, argv for main
		- running main on each node, with short wait between nodes
		- mpi::finalize
----------------------------------------------------------------------------------------------------------------------------*/

int argc_main = 11;
vector <string> argv_main(argc_main);
argv_main[0] = "main";
argv_main[1] = "-mintn";	argv_main[2] = timenumber;
argv_main[3] = "-maxtn";	argv_main[4] = timenumber;
argv_main[5] = "-minll";	argv_main[6] = loop;
argv_main[7] = "-maxll";	argv_main[8] = loop;
argv_main[9] = "-opts";		argv_main[10] = "data/00"+numberToString<int>(rank)+"optionsM";

if (rank==0 && revertToDefault) {
	argc_main += 2;
	argv_main.push_back("-inF");
	argv_main.push_back("p");
}

sleep(rank*2);

cout << "node " << rank << " of " << nodes << " running main with: timenumber " << timenumber;
cout << ", loop " << loop << endl;

returnValue = main_fn(argc_main,argv_main);

if (returnValue!=0) {
	cerr << "return " << returnValue << " for node " << rank << " on running main" << endl;
}

MPI::Finalize();

return 0;
}

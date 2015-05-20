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

#define MASTER 0         // rank of master node

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - defining required nodes
		2 - copying files
		3 - finding recent files
		3 - initializing mpi
		4 - running main
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/*----------------------------------------------------------------------------------------------------------------------------
	1. defining required nodes
		- initializing mpi
----------------------------------------------------------------------------------------------------------------------------*/

int nodes_req = 2;

int nodes, rank;
int returnValue = 0;

MPI::Init(argc, argv);
MPI_Comm_size(MPI_COMM_WORLD, &nodes);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
if (nodes==nodes_req) {
	if (rank==MASTER) {
		cout << "running " << nodes << " nodes" << endl;
	}
}
else {
	if (rank==MASTER) {
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

bool copyFiles = true;

if (copyFiles && rank==MASTER) {
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
----------------------------------------------------------------------------------------------------------------------------*/

// getting last timenumbers
bool revertToDefault = false;
vector<string> timenumbers(nodes_req);
vector<string> loops(nodes_req);
char* timenumber;
char* loop;

FilenameAttributes fa_low;
fa_low.Directory = "data";
fa_low.Suffix = ".data";
fa_low.ID = "mainp";
fa_low.Timenumber = "0";
FilenameAttributes fa_high(fa_low);
fa_high.Timenumber = "999999999999";
Folder F(fa_low,fa_high);

if (F.size()==0)
	revertToDefault = true;
else if (rank==MASTER) {
	for (int k=0; k<nodes_req; k++) {
		string maxTimenumber = "0";
		uint posMax = 0;
		if (F.size()==0) {
			revertToDefault = true;
			break;
		}
		for (uint j=0; j<F.size(); j++) {
			if ((F[j]).Timenumber>maxTimenumber) {
				maxTimenumber = (F[j]).Timenumber;
				posMax = j;
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

if (revertToDefault && rank==MASTER) {
	for (int k=0; k<nodes_req; k++) {
		timenumbers[k] = "00"+numberToString<int>(k);
		loops[k] = "0";
	}
}

/*if (rank==MASTER) {
	timenumber = new char[12];
	loop = new char[12];
	for (int k=0; k<nodes_req; k++) {
		*timenumber = *(timenumbers[k]).c_str();
		*loop = *(loops[k]).c_str();
		MPI_Send(timenumber,(timenumbers[k]).size(), MPI_BYTE, k, 0, MPI_COMM_WORLD);
		MPI_Send(loop,(timenumbers[k]).size(), MPI_BYTE, k, 1, MPI_COMM_WORLD);
		cout << "process " << MASTER << " sent " << timenumber << " and " << loop << " to process " << k << endl;
	}
	delete loop;
	delete timenumber;
}
else {
	MPI_Recv(timenumber, 12, MPI_BYTE, MASTER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Recv(loop, 12, MPI_BYTE, MASTER, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	cout << "process " << rank << " received " << timenumber << " and " << loop << " from process " << MASTER << endl;
}*/


/*----------------------------------------------------------------------------------------------------------------------------
	4. initializing mpi
		- defining argc, argv for main
		- mpi::init
		- checking required number of nodes
----------------------------------------------------------------------------------------------------------------------------*/

int argc_main = 9;
vector <string> argv_main(argc_main);
argv_main[0] = "main";
argv_main[1] = "-mintn";
argv_main[3] = "-maxtn";
argv_main[5] = "-loopMin";
argv_main[7] = "-loopMax";

/*----------------------------------------------------------------------------------------------------------------------------
	5. running main
		- running main on each node, with short wait between nodes
		- mpi::finalize
----------------------------------------------------------------------------------------------------------------------------*/

argv_main[2] = (string)timenumber;
argv_main[4] = (string)timenumber;
argv_main[6] = (string)loop;
argv_main[8] = (string)loop;

sleep(rank*2);

cout << "node " << rank << " of " << nodes << " running main with: timenumber " << timenumbers[rank];
cout << ", loop " << loops[rank] << endl;

returnValue = main_fn(argc_main,argv_main);

if (returnValue!=0) {
	cerr << "return " << returnValue << " for node " << rank << " on running main" << endl;
}

MPI::Finalize();

return 0;
}

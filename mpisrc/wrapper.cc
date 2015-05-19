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
		3 - initializing mpi
		4 - running main
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/*----------------------------------------------------------------------------------------------------------------------------
	1. defining required nodes
----------------------------------------------------------------------------------------------------------------------------*/

int nodes_req = 2;

/*----------------------------------------------------------------------------------------------------------------------------
	2. copying files
		- copying optionsM
		- copying negEig, omega etc
----------------------------------------------------------------------------------------------------------------------------*/

bool copyFiles = false;

if (copyFiles) {
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

FilenameAttributes fa_low;
fa_low.Directory = "data";
fa_low.Suffix = ".data";
fa_low.ID = "mainp";
fa_low.Timenumber = "0";
FilenameAttributes fa_high(fa_low);
fa_high.Timenumber = "999999999999";
FilenameComparator fc(fa_low,fa_high);
Folder F(fc);

if (F.size()==0)
	revertToDefault = true;
else {
	for (int k=0; k<nodes_req; k++) {
		cout << F << endl;
		string maxTimenumber = "0";
		uint posMax = 0;
		if (F.size()==0) {
			revertToDefault = true;
			break;
		}
		for (uint j=0; j<F.size(); j++) {
			if ((F[j]).Timenumber>maxTimenumber) {
				maxTimenumber = (F[j]).Timenumber;
			}
		}
		timenumbers[k] = (F[posMax]).Timenumber;
		long long unsigned int l = stringToNumber<long long unsigned int>((F[posMax]).Timenumber);
		l--;
		fa_high.Timenumber = numberToString<long long unsigned int>(l);
		fc.set(fa_low,fa_high);
		F.set(fc);
	}

	// getting last loops
	if (!revertToDefault) {
		for (int k=0; k<nodes_req; k++) {
			fa_low.Timenumber = timenumbers[k];
			fa_high.Timenumber = timenumbers[k];
			fc.set(fa_low,fa_high);
			F.set(fc);
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

if (revertToDefault) {
	for (int j=0; j<nodes_req; j++) {
		timenumbers[j] = "00"+numberToString<int>(j);
		loops[j] = "0";
	}
}

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
	5. running main
		- running main on each node, with short wait between nodes
		- mpi::finalize
----------------------------------------------------------------------------------------------------------------------------*/

argv_main[2] = timenumbers[rank];
argv_main[4] = timenumbers[rank];
argv_main[6] = loops[rank];
argv_main[8] = loops[rank];

sleep(rank*2);

cout << "node " << rank << " of " << nodes << " running main with input timenumber 00" << rank << endl;

returnValue = main_fn(argc_main,argv_main);

if (returnValue!=0) {
	cerr << "return " << returnValue << " for node " << rank << " on running main" << endl;
}

MPI::Finalize();

return 0;
}

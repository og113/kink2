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
		2 - organizing files
		3 - initializing mpi
		4 - running main
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {
/*----------------------------------------------------------------------------------------------------------------------------
	1. defining required nodes
----------------------------------------------------------------------------------------------------------------------------*/

int nodes_req = 8;

/*----------------------------------------------------------------------------------------------------------------------------
	2. organizing files
		- copying optionsM
		- copying negEig, omega etc
----------------------------------------------------------------------------------------------------------------------------*/

bool copyFiles = false;

if (copyFiles) {
	for (int k=0; k<nodes_req; k++) {
		string timenumber = "00"+numberToString<int>(k);
		int returnValue;
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
	3. initializing mpi
		- defining argc, argv for main
		- mpi::init
		- checking required number of nodes
----------------------------------------------------------------------------------------------------------------------------*/

int argc_main = 3;
vector <string> argv_main(argc_main);
argv_main[0] = "main";
argv_main[1] = "-opts";

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
	4. running main
		- running main on each node, with short wait between nodes
		- mpi::finalize
----------------------------------------------------------------------------------------------------------------------------*/

argv_main[2] = "data/00"+numberToString<int>(rank)+"optionsM";

sleep(rank*2);

cout << "node " << rank << " of " << nodes << " running main with input timenumber 00" << rank << endl;

returnValue = main_fn(argc_main,argv_main);

if (returnValue!=0) {
	cerr << "return " << returnValue << " for node " << rank << " on running main" << endl;
}

MPI::Finalize();

return 0;
}

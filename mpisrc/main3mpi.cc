/*----------------------------------------------------------------------------------------------------------------------------
	main3mpi
		mpi wrapper for main3mpi
----------------------------------------------------------------------------------------------------------------------------*/

#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cstdlib>
#include <string>
#include <vector>
#include "parameters.h"
#include "simple.h"
#include "main_fn3.h"

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
	contents:
		0 - initializing mpi
		1 - sorting inputs
		2 - running nrmain
		3 - clearing up
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

using namespace std;

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	0. initializing mpi
----------------------------------------------------------------------------------------------------------------------------*/

int Nw, rank, root = 0; // Nw, number of workers

MPI_Init(&argc, &argv);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &Nw);

if (rank==root)
	cout << "starting main3mpi with " << Nw << " nodes" << endl;	

/*----------------------------------------------------------------------------------------------------------------------------
	1. sorting inputs
----------------------------------------------------------------------------------------------------------------------------*/

sleep(rank*2);

// argc and argv
uint argcSerial = 1;
string mpiInputsFile = "nrinputs/mpi/inputs0";
string argvList = "";
vector<string> argvSerial(1);
argvSerial[0] = "./nrmpi";

// printing argv
if (rank==root) {
	cout << "input argv to ./main3mpi:" << endl;
	for (int j=0; j<argc; j++) {
		cout << argv[j];
		if (j<(argc-1))
			cout << " ";
		else
			cout << endl;
	}
}

// getting mpi and serial argv
if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("mpiinputs")==0 || id.compare("mpiInputs")==0 || id.compare("inputs")==0)
			mpiInputsFile = (string)argv[2*j+2];
		else if (id.compare("argvlist")==0 || id.compare("argvList")==0)
			argvList = (string)argv[2*j+2];
		else {
			argvSerial.push_back((string)argv[2*j+1]);
			argvSerial.push_back((string)argv[2*j+2]);
			argcSerial += 2;
		}
	}
}
else if (rank==root) {
	cerr << "must provide an even number of arguments after ./main3mpi" << endl;
	MPI_Abort(MPI_COMM_WORLD,1);
}

if (argvList.empty()) {
	
	// loading parameters
	ParametersRange prMpi;
	prMpi.load(mpiInputsFile);
	if (rank==root) {
		Parameters pMpi = prMpi.Min;
		if (pMpi.empty()) {
			cerr << "Parameters empty: nothing in inputs file: " << mpiInputsFile << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
		if (Nw>(int)prMpi.totalSteps()) {
			cerr << "Error: too few steps, " << prMpi.totalSteps() << ", for nodes, " << Nw << endl;
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}

	// distributing parameters
	string inputsSerial = mpiInputsFile+"_rank_"+nts(rank);
	ParametersRange prSerial = prMpi;
	prSerial.Min = prMpi.position(rank);
	uint j=0;
	bool foundToStep = false;
	while (j<Parameters::Size && !foundToStep) {
		if (prMpi.toStep((Parameters::Label)(j+1))) {
			(prSerial.Steps)[j] = 0;
			if (rank==root && (int)(prMpi.Steps)[j]!=Nw) {
				cerr << "Error: number of steps on first parameter, " << (prMpi.Steps)[j];
				cerr << ", not equal to number of nodes, " << Nw << endl;
			}
			foundToStep = true;
		}
		j++;
	}
	prSerial.save(inputsSerial);

	// adding inputs filename to argv
	argvSerial.push_back("-inputs");
	argvSerial.push_back(inputsSerial);
	argcSerial += 2;
}
else {
	ifstream is;
	is.open(argvList.c_str());
	if (is.good()) {
		string line;
		for (int j=0; j<(rank+1); j++) {
			if (is.eof()) {
				cerr << "main3mpi error: argvList file " << argvList << " has too few lines for rank " << rank << endl;
				return 1;
			}
			getline(is,line);
		}
		is.close();
		argvSerial.clear();
		argvSerial = splitString(line," ");
		argcSerial = argvSerial.size();
	}
	else {
		cerr << "main3mpi error: cannot load from " << argvList << endl;
		MPI_Abort(MPI_COMM_WORLD,1);
	}
}

/*----------------------------------------------------------------------------------------------------------------------------
	2. running nrmain
----------------------------------------------------------------------------------------------------------------------------*/

uint returnValue = 0;
returnValue = main_fn3(argcSerial,argvSerial);

if (returnValue!=0) {
	cerr << "----------------------------------------------------------------------" << endl;
	cerr << "return " << returnValue << " for node " << rank << " on running nrmain" << endl;
	cerr << "----------------------------------------------------------------------" << endl;
}

/*----------------------------------------------------------------------------------------------------------------------------
	3. clearing up
----------------------------------------------------------------------------------------------------------------------------*/

MPI_Finalize();

return 0;
}

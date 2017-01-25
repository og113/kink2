/*----------------------------------------------------------------------------------------------------------------------------
	main_fn3
		program to solve boundary value problem on contour ABCD
		main used as a function, not a main process
----------------------------------------------------------------------------------------------------------------------------*/

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include "analysis.h"
#include "eigen_extras.h"
#include "main.h"
#include "main_fn3.h"
#include "print3.h"

//#define NDEBUG //NDEBUG is to remove error and bounds checking on vectors in SparseLU, for speed - only include once everything works

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		0 - enums
		1 - argv inputs, loading options, closenesses
		2 - ParametersRange, Stepper
		3 - beginning parameter loop
		4 - declaring some basic quantities
		5 - assigning potential functions
		6 - omega and negVec
		7 - defining quantitites
		8 - beginning newton-raphson loop
		9 - assigning minusDS, DDS etc
		10 - checks
		11 - printing early 1
		12 - solving for delta
		13 - printing early 2
		14 - convergence
		15 - printing output
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------
	0. enums
----------------------------------------------------------------------------------------------------------------------------*/

struct StepperArgv {
	enum Option { none, w, e, n};
};

/*----------------------------------------------------------------------------------------------------------------------------
	1. loading options, argv inputs
		- defining timenumber and files to load
		- argv inputs
		- loading options
		- loading closenesses
		- beginning cos and ces 
----------------------------------------------------------------------------------------------------------------------------*/

int main_fn3(int argc, vector<string> argv)
{

// defining timenumber
string timenumber = currentDateTime();

// defining rank
string rank = "0";

// options to load
Options opts;
Closenesses tols;
string optionsFile = "nrinputs/options0";
string inputsFile = "nrinputs/inputs0";
string tolsFile = "nrinputs/tols0";
string fIn = "";

// cos and cerr files
string coFile, ceFile;

// stepper stuff
string stepperArgv = "";
string stepperInputsFile = "step0";
Filename stepperOutputsFile = "results/nr/stepper/step.csv";
uint steps = 1;

// getting argv inputs
if (argc==2) timenumber = argv[1];
else if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("opts")==0 || id.compare("options")==0) 						optionsFile = argv[2*j+2];
		else if (id.compare("inputs")==0) 											inputsFile = argv[2*j+2];
		else if (id.compare("tols")==0 || id.compare("tolerances")==0) 				tolsFile = argv[2*j+2];
		else if (id.compare("co")==0) 												coFile = argv[2*j+2];
		else if (id.compare("ce")==0) 												ceFile = argv[2*j+2];
		else if (id.compare("fIn")==0 || id.compare("fin")==0) 						fIn = argv[2*j+2];
		else if (id.compare("stepper")==0) 											stepperArgv = argv[2*j+2];
		else if (id.compare("steps")==0) 											steps = stn<uint>(argv[2*j+2]);
		else if (id.compare("stepperInputs")==0 || id.compare("stepInputs")==0) 	stepperInputsFile = argv[2*j+2];
		else if (id.compare("stepperOutputs")==0 || id.compare("stepResults")==0)	stepperOutputsFile = argv[2*j+2];
	}
}
else if (argc != 1) {
	cerr << "must provide an even number of inputs in format '-name value':" << endl;
	for (int j=1; j<argc; j++) cerr << argv[j] << " ";
	cerr << endl;
	return 1;
}

// loading opts
opts.load(optionsFile);
//opts.print();

// loading closenesses
tols.load(tolsFile);

// other possible inputs
bool trivialChecks = false;
bool verbose = true;
bool redo = true;
bool redoErrors = true;
string baseFolder = "";

// getting argv inputs
if (argc==2) timenumber = argv[1];
else if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("opts")==0 || id.compare("options")==0);
		else if (id.compare("inputs")==0);
		else if (id.compare("tols")==0 || id.compare("tolerances")==0);
		else if (id.compare("co")==0);
		else if (id.compare("ce")==0);
		else if (id.compare("fIn")==0 || id.compare("fin")==0);
		else if (id.compare("stepper")==0);
		else if (id.compare("steps")==0);
		else if (id.compare("stepperInputs")==0 || id.compare("stepInputs")==0);
		else if (id.compare("stepperOutputs")==0 || id.compare("stepResults")==0);
		else if (id.compare("tn")==0 || id.compare("timenumber")==0) 			timenumber = argv[2*j+2];
		else if (id.compare("zmx")==0) 											opts.zmx = argv[2*j+2];
		else if (id.compare("zmt")==0) 											opts.zmt = argv[2*j+2];
		else if (id.compare("bds")==0) 											opts.bds = argv[2*j+2];
		else if (id.compare("inF")==0) 											opts.inF = argv[2*j+2];
		else if (id.compare("epsiTb")==0) 										opts.epsiTb = stn<double>(argv[2*j+2]);
		else if (id.compare("epsiTheta")==0) 									opts.epsiTheta = stn<double>(argv[2*j+2]);
		else if (id.compare("loops")==0) 										opts.loops = stn<uint>(argv[2*j+2]);
		else if (id.compare("printChoice")==0) 									opts.printChoice = argv[2*j+2];
		else if (id.compare("rank")==0) 										rank = argv[2*j+2];
		else if (id.compare("verbose")==0) 										verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("trivialChecks")==0) 								trivialChecks = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("redo")==0) 										redo = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("redoErrors")==0) 									redoErrors = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("baseFolder")==0) 									baseFolder = argv[2*j+2];
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

// filling empty co and ce filenames
if (coFile.empty()) coFile = "data/coe/"+timenumber+"co.txt";
if (ceFile.empty()) ceFile = "data/coe/"+timenumber+"ce.txt";

// beginning cos and ces streams
fstream cos;
cos.open(coFile.c_str(),fstream::app);

fstream ces;
ces.open(ceFile.c_str(),fstream::app);

// stepper
StepperArgv::Option stepargv = StepperArgv::none;
if (!stepperArgv.empty()) {
	if (stepperArgv.compare("W")==0) {
		stepargv = StepperArgv::w;
		(stepperOutputsFile.Extras).push_back(StringPair("const","W"));
		cos << "stepping with constant W" << endl;
		if (verbose)
			cout << "stepping with constant W" << endl;
	}
	else if (stepperArgv.compare("E")==0) {
		stepargv = StepperArgv::e;
		(stepperOutputsFile.Extras).push_back(StringPair("const","E"));
		cos << "stepping with constant E" << endl;
		if (verbose)
			cout << "stepping with constant E" << endl;
	}
	else if (stepperArgv.compare("N")==0) {
		stepargv = StepperArgv::n;
		(stepperOutputsFile.Extras).push_back(StringPair("const","N"));
		cos << "stepping with constant N" << endl;
		if (verbose)
			cout << "stepping with constant N" << endl;
	}
	else if (stepperArgv.compare("none")==0)
		stepargv = StepperArgv::none;
	else {
		ces << endl << "loopChoice and stepper options not possible:" << endl;
		ces << opts << endl;
		cerr << endl << "loopChoice and stepper options not possible:" << endl;
		cerr << opts << endl;
		return 1;
	}
}

/*----------------------------------------------------------------------------------------------------------------------------
	2. ParametersRange, Stepper
----------------------------------------------------------------------------------------------------------------------------*/

// ParametersRange
ParametersRange pr;
pr.load(inputsFile);
Parameters p = pr.Min, pold = pr.Min;
if (p.empty()) {
	ces << "Parameters empty: nothing in inputs file: " << inputsFile << endl;
	cerr << "Parameters empty: nothing in inputs file: " << inputsFile << endl;
	return 1;
}

// initializing stepper
StepperOptions stepOpts;
stepOpts.tol = 1.0;
Point2d point;
if (stepargv!=StepperArgv::none) {
	ifstream is;
	is.open(stepperInputsFile.c_str());
	if (is.good()) {
		is >> stepOpts;
		is.close();
		(stepperOutputsFile.Extras).push_back(StringPair("tol",nts(stepOpts.tol)));
		(stepperOutputsFile.Extras).push_back(StringPair("aim",nts(stepOpts.aim)));
	}
	else {
		ces << "Error: cannot open stepper inputs file, " << stepperInputsFile << endl;
		cerr << "Error: cannot open stepper inputs file, " << stepperInputsFile << endl;
		return 1;
	}
	point(p.Tb,p.Theta);
}
Stepper stepper(stepOpts,point);
if (stepperOutputsFile.exists() && stepargv!=StepperArgv::none) {
	stepper.load(stepperOutputsFile);
	p.Tb = stepper.x();
	p.Theta = stepper.y();
	pold.Tb = (stepper.lastStep()).X;
	pold.Theta = (stepper.lastStep()).Y;
}

//stepOpts.closeness = closenesses.Step; // NOT USING THIS

// results
string resultsFile = (pass? "results/nr/nr0pass.csv":"results/nr/nr0.csv");
uint idSizeResults = 3, datumSizeResults = 17; // THIS WILL NEED CHANGING
vector<string> idCheck(idSizeResults);
idCheck[idSizeResults-1] = potExtras.second;
NewtonRaphsonData results(resultsFile,idSizeResults,datumSizeResults);

// errors
string errorsFile = "results/nr/nr0error.csv";
uint idSizeErrors = 4, datumSizeErrors = 13;  // THIS WILL NEED CHANGING
vector<string> idCheckErrors(idSizeErrors);
idCheckErrors[idSizeErrors-1] = potExtras.second;
NewtonRaphsonData errors(errorsFile,idSizeErrors,datumSizeErrors);

// filename suffix
string suffix = ((opts.inF).back()=='n'? ".data": ".dat");

// printing timenumber and pot
cos << "timenumber: " << timenumber << ", pot: " << potExtras.second << endl;
cout << "timenumber: " << timenumber << ", pot: " << potExtras.second << endl;

/*----------------------------------------------------------------------------------------------------------------------------
	3. beginning parameter loop
		- defining a time
		- changing parameters (if required)
		- copying a verson of parameters with timenumber
		- printing timenumber
		- declaring checks
----------------------------------------------------------------------------------------------------------------------------*/

// number of loops
uint Npl = (stepargv==StepperArgv::none? pr.totalSteps(): steps);
cos << "looping over " << Npl << " steps" << endl;
if (verbose)
	cout << "looping over " << Npl << " steps" << endl;

// starting loop
for (uint pl=0; pl<Npl; pl++) {
//defining a time and starting the clock
	clock_t time;
	time = clock();
	
	// doing things before parameter step
	Filename loadFile, stepFile;
	
	// stepping parameters
	if (pl>0) {
		if (stepargv==StepperArgv::none) {
			p = pr.position(pl);
			pold = pr.neigh(pl);
		}
		else {
			// getting step base
			p.Tb = stepper.x();
			p.Theta = stepper.y();
			pold.Tb = (stepper.lastStep()).X;
			pold.Theta = (stepper.lastStep()).Y;
		}
	}
	
	// checking if have results already
	if (!redo && results.find(idCheck,p)) {
		if (verbose) {
			cout << "result found in " << resultsFile << " for pl = " << pl << ", ";
			cout << "continuing to next step" << endl;
		}
		continue;
	}
	if (!redoErrors && errors.find(idCheckErrors,p)) {
		if (verbose) {
			cout << "result found in " << errorsFile << " for pl = " << pl << ", ";
			cout << "continuing to next step" << endl;
		}
		continue;
	}
	
	// getting step file
	if (pl>0) {
		if (opts.InF[0]=='m') {
			stepFile = filenameMain(pold,baseFolder,"field","fmain",suffix);
		}
		else if (opts.InF[0]=='p')
			stepFile = filenameMain(pold,baseFolder,"field","fpi",suffix);
		else {
			ces << "opts.inF not understood: " << opts.inF << endl;
			cerr << "opts.inF not understood: " << opts.inF << endl;
			return 1;
		}
	}
			
	// opening file to print c-style
	cos.close();
	FILE* cof;
	cof = fopen(coFile.c_str(),"a");
	
	//printing angle
	if (stepargv!=StepperArgv::none) {
		double angleModTwoPi = mod(stepper.stepAngle(),-pi,pi);		
		fprintf(cof,"%12s%12.3g\n","step angle: ",angleModTwoPi);
		if (verbose)
			printf("%12s%12.3g\n","step angle: ",angleModTwoPi);
	}

/*----------------------------------------------------------------------------------------------------------------------------
4. declaring some basic quantities
	- checks
	- some rerived quantities
----------------------------------------------------------------------------------------------------------------------------*/

	// declaring Checks
	Check checkAction("action",closenesses.Action);
	Check checkSoln("solution",closenesses.Soln);
	Check checkSolnMax("solution max",closenesses.SolnMax);
	Check checkDelta("delta",closenesses.Delta);
	Check checkInv("matrix inversion",closenesses.Inv*p.N*p.NT);
	Check checkCon("energy conservation",closenesses.Con);
	Check checkLin("linear energy flat",closenesses.Lin);
	Check checkTrue("linear energy equal true energy",closenesses.True);
	Check checkLatt("lattice small enough for energy",closenesses.Latt);
	Check checkReg("regularisation term",closenesses.Reg);
	Check checkIE("imaginary part of energy",closenesses.IE);
	Check checkContm("linear energy equal continuum expression",closenesses.Contm);
	Check checkOS("linear energy equal on shell expression",closenesses.OS);
	Check checkAB("a_k = Gamma*b_k",closenesses.AB);
	Check checkABNE("N = Sum(a_k*b_k), E = Sum(w_k*a_k*b_k)",closenesses.ABNE);
	Check checkLR("linear representation of phi",closenesses.LR);
	Check checkBoundRe("initial real boundary condition",closenesses.Soln);
	Check checkBoundIm("initial imaginary boundary condition",closenesses.Soln);
	Check checkChiT("time zero mode condition",closenesses.Soln);
	
	// some derived quantities
	uint zm = (p.Pot<3? 2: 1); // number of zero modes

/*----------------------------------------------------------------------------------------------------------------------------
5. assigning potential functions
	- assigning potential functions
	- assigning preliminary parameter structs
----------------------------------------------------------------------------------------------------------------------------*/

	// assigning potential functions
	Potential<comp> V, dV, ddV;
	if (p.Pot==1) {
		V((Potential<comp>::PotentialType)&V1<comp>,ps);
		dV((Potential<comp>::PotentialType)&dV1<comp>,ps);
		ddV((Potential<comp>::PotentialType)&ddV1<comp>,ps);
	}
	else if (p.Pot==2) {
		V((Potential<comp>::PotentialType)&V2<comp>,ps);
		dV((Potential<comp>::PotentialType)&dV2<comp>,ps);
		ddV((Potential<comp>::PotentialType)&ddV2<comp>,ps);
	}
	else if (p.Pot==3) {
		V((Potential<comp>::PotentialType)&V3<comp>,ps);
		dV((Potential<comp>::PotentialType)&dV3<comp>,ps);
		ddV((Potential<comp>::PotentialType)&ddV3<comp>,ps);
		}
	else {
		ces << "pot option not available, pot = " << p.Pot << endl;
		if ((opts.printChoice).compare("gui")==0)
			cerr << "pot option not available, pot = " << p.Pot << endl;
		return 1;
	}

	// assigning preliminary parameter structs
	params_for_V paramsV  = {p.epsilon, p.A};

	//lambda functions for pot_r
	auto Vr = [&] (const comp& phi) {
		return -ii*p.Reg*VrFn(phi,p.minima[0],p.minima[1]);
	};
	auto dVr = [&] (const comp& phi) {
		return -ii*p.Reg*dVrFn(phi,p.minima[0],p.minima[1]);
	};
	auto ddVr = [&] (const comp& phi) {
		return -ii*p.Reg*ddVrFn(phi,p.minima[0],p.minima[1]);
	};

/*----------------------------------------------------------------------------------------------------------------------------
6. omega and negVec
	- loading or constructing omega
	- loading negVec
----------------------------------------------------------------------------------------------------------------------------*/

	//deterimining omega matrices for fourier transforms in spatial direction
	vec freqs(p.N), freqs_exp(p.N);
	mat modes(p.N,p.N);
	mat omega_m1(p.N,p.N), omega_0(p.N,p.N), omega_1(p.N,p.N), omega_2(p.N,p.N);
	
	omegaM1F = filenameMain(p,baseFolder,"omega","omegaM1",suffix);
	omega0F = filenameMain(p,baseFolder,"omega","omega0F",suffix);
	omega1F = filenameMain(p,baseFolder,"omega","omega1F",suffix);
	omega2F = filenameMain(p,baseFolder,"omega","omega2F",suffix);
	modesF = filenameMain(p,baseFolder,"omega","modesF",suffix);
	freqsF = filenameMain(p,baseFolder,"omega","freqsF",suffix);
	freqsExpF = filenameMain(p,baseFolder,"omega","freqsExpF",suffix);
	
	if (omegaM1F.exists() && omega0F.exists() && omega1F.exists() && omega2F.exists() \
			&& omegaM1F.modesF() && freqsF.exists() && freqsExpF.exists()) {
			loadMatrixBinary(omegaM1F,omega_m1);
			loadMatrixBinary(omega0F,omega_0);
			loadMatrixBinary(omega1F,omega_1);
			loadMatrixBinary(omega2F,omega_2);
			loadMatrixBinary(modesF,modes);
			loadVectorBinary(freqsF,freqs);
			loadVectorBinary(freqsExpF,freqs_exp);
		}
		else {
			bool approxOmega = true;
			if (!approxOmega) {
				numericalModes(modes,freqs,freqs_exp,ps);
			}
			else {
				analyticModes(modes,freqs,freqs_exp,ps);
			}
			omegasFn(approxOmega,modes,freqs,omega_m1,omega_0,omega_1,omega_2,ps);
			saveMatrixBinary(omegaM1F,omega_m1);
			saveMatrixBinary(omega0F,omega_0);
			saveMatrixBinary(omega1F,omega_1);
			saveMatrixBinary(omega2F,omega_2);
			saveMatrixBinary(modesF,modes);
			saveVectorBinary(freqsF,freqs);
			saveVectorBinary(freqsExpF,freqs_exp);
		}

//############################################################ CURRENT STATE ############################################################################

	// loading negVec
	vec negVec;
	if (opts.zmt[0]=='n' || opts.zmx[0]=='n') {
		if (p.Pot==3) {
			Filename eigVecFile = "data/00"+rank+"eigVec_pot_3_L_" + numberToString<double>(p.L) + ".dat";	
			so_simple.printType = SaveOptions::ascii;
			load(eigVecFile,so_simple,negVec); // should automatically interpolate
			so_simple.printType = SaveOptions::binary;
		}
		else {
			Filename eigVecFile;
			string N_load;
			string Nb_load;
			Filename lower = "data/00"+rank+"eigVec_pot_" + numberToString<uint>(p.Pot)\
								 + "_N_100_Nb_100_L_" + numberToString<double>(p.L) + ".dat";
			Filename upper = lower;
			vector<StringPair> upperExtras = lower.Extras;
			upper.Extras[1] = StringPair("N","1000");
			upper.Extras[2] = StringPair("Nb","1000");
			Folder eigVecFolder(lower,upper);
			if (eigVecFolder.size()==0) {
				ces << "no negative eigenvector files found between:" << endl << lower << endl << upper << endl;
				if ((opts.printChoice).compare("gui")==0)
					cerr << "no negative eigenvector files found between:" << endl << lower << endl << upper << endl;
				return 1;
			}
			else {
				eigVecFile = eigVecFolder[0];
				for (uint j=1; j<eigVecFolder.size(); j++) {
					// picking the file with the largest N, for some reason
					if (((eigVecFolder[j].Extras)[1]).second>((eigVecFile.Extras)[1]).second) {
						eigVecFile = eigVecFolder[j];
						N_load = ((eigVecFile.Extras)[1]).second;
						Nb_load = ((eigVecFile.Extras)[2]).second;
					}
				}
			}
			SaveOptions eigVecOpts;
			eigVecOpts.printType = SaveOptions::ascii;
			eigVecOpts.vectorType = SaveOptions::realB;
			eigVecOpts.extras = SaveOptions::coords;
			eigVecOpts.paramsOut = psu; // USE OF PSU
			eigVecOpts.paramsOut = ps;
			eigVecOpts.printMessage = false;
			Parameters pIn = ps;
			pIn.N = stn<uint>(N_load);
			pIn.Nb = stn<uint>(Nb_load);
			pIn.NT = 1000; // a fudge so that interpolate realises the vector is only on BC
			load(eigVecFile,eigVecOpts,negVec);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------------
7. defining quantities
	- erg, linErg etc
	- f, minusDS, DDS
	- printing parameters
	- print input phi
----------------------------------------------------------------------------------------------------------------------------*/	

	//defining energy and number vectors
	cVec erg(p.NT);
	cVec linErg(p.NT);
	cVec linNum(p.NT);
	cVec linErgOffShell(p.NT);
	cVec linNumOffShell(p.NT);
	cVec derivErg(p.NT), potErg(p.NT);
	comp linErgContm, linNumContm;
	comp linNumAB, linErgAB;
	
	//defining the action and bound and W and zero of energy
	double ergZero = (p.Pot==3? 0.0: p.N*p.a*real(V(p.minima[0])) );
	comp action = ii*p.action0;
	double bound = 0.0;
	double W = 0.0;
	double E = 0.0;
	double E_exact = 0.0;
	double Num = 0.0;
	
	//defining some quantities used to stop the Newton-Raphson loop when action stops varying
	comp action_last = action;
	uint runs_count = 0;
	uint min_runs = 1;
	uint max_runs = 100;

	//initializing phi (=f)
	vec f;
	SaveOptions so_tp;
	so_tp.printType = SaveOptions::binary;
	so_tp.paramsIn = psu;
	so_tp.paramsOut = ps;
	so_tp.vectorType = SaveOptions::complex;
	so_tp.extras = SaveOptions::coords;
	so_tp.zeroModes = 2;
	so_tp.printMessage = false;
	if (loop==0) {
		if (fileLoop==0 && ((pFolder[0]).Suffix).compare(".dat")==0) so_tp.printType = SaveOptions::ascii;
		load(pFolder[fileLoop],so_tp,f); //n.b there may be some problems with zero modes for binary printing
		if (f.size()==(2*p.N*p.NT+1)) { // if came from pi.cc and in binary
			f.conservativeResize(2*p.N*p.NT+2);
			f(2*p.N*p.NT+1) = 0.5;
		}
		so_tp.paramsIn = ps;
		so_tp.printType = SaveOptions::binary;
		fprintf(cof,"%12s%30s\n","input: ",(pFolder[fileLoop]()).c_str());
		if ((opts.printChoice).compare("gui")==0)
			printf("%12s%30s\n","input: ",(pFolder[fileLoop]()).c_str());
	}
	else {
		Filename lastPhi;
		uint sigma = (stepper.steps()>0? 1:0);
		if (((opts.loopChoice).substr(0,5)).compare("const")!=0) {
			lastPhi = (string)("data/"+timenumber+"mainp_fLoop_"+numberToString<uint>(fileLoop)\
					+"_loop_"+numberToString<uint>(loop-1)+".data");
		}
		else {
			lastPhi = (string)("data/"+timenumber+"mainp_fLoop_"+numberToString<uint>(fileLoop)\
					+"_loop_"+numberToString<uint>(loop-stepper.local()+1+sigma)+"_step_1.data");
		}
		load(lastPhi,so_tp,f);
		fprintf(cof,"%12s%30s\n","input: ",(lastPhi()).c_str());
		if ((opts.printChoice).compare("gui")==0)
			printf("%12s%30s\n","input: ",(lastPhi()).c_str());
	}
	
	// printing parameters
	p.print(cof);
	if ((opts.printChoice).compare("gui")==0)
		p.print();
	
	//defining complexified vector Cf
	cVec Cf;
	Cf = vecComplex(f,ps);

	//defining DDS and minusDS
	spMat DDS(2*p.N*p.NT+2,2*p.N*p.NT+2);
	vec minusDS(2*p.N*p.NT+2);
	
	// defining a couple of vectors to test whether the initial boundary conditions are satisfied
	vec boundRe(p.N);
	vec boundIm(p.N);
		
	//very early vector print
	so_tp.paramsIn = ps;
	if ((opts.printChoice).compare("n")!=0) {
		Filename earlyPrintFile = (string)("data/"+timenumber+"mainpE_fLoop_"+numberToString<uint>(fileLoop)\
				 +"_loop_"+numberToString<uint>(loop)+"_run_" + "0.data");
		save(earlyPrintFile,so_tp,f);
	}
/*----------------------------------------------------------------------------------------------------------------------------
8. beginning newton-raphson loop
	- beginning newton-raphson loop	
	- chiT, chiX
	- reserving space in DDS
	- zeroing erg etc
	- crude test of potential (to remove once done)
----------------------------------------------------------------------------------------------------------------------------*/

	//beginning newton-raphson loop	
	while (!checkSoln.good() || !checkSolnMax.good() || runs_count<min_runs) {
		runs_count++;
		if (runs_count>max_runs) {
			ces << "main error: max_runs(" << max_runs << ") exceeded for:" << endl;
			ces << "timenumber: " << timenumber << "; fileLoop: " << fileLoop << "; loop: " << loop << endl;
			if ((opts.printChoice).compare("gui")==0) {
				cerr << "main error: max_runs(" << max_runs << ") exceeded for:" << endl;
				cerr << "timenumber: " << timenumber << "; fileLoop: " << fileLoop << "; loop: " << loop << endl;
			}
			return 1;				
		}
		
		// zero modes - fixed with chiX and chiT
		vec chiX(p.NT*p.N);	chiX = Eigen::VectorXd::Zero(p.N*p.NT); //to fix spatial zero mode
		vec chiT(p.NT*p.N);	chiT = Eigen::VectorXd::Zero(p.N*p.NT); //to fix real time zero mode
		for (uint j=0; j<p.N; j++) {
			uint posX, posT, posCe;
			uint slicesX, slicesT;
			if (getLastInt(opts.zmx)<0) {
				ces << "getLastInt error with zmx = " << opts.zmx << endl;
				if ((opts.printChoice).compare("gui")==0)
					cerr << "getLastInt error with zmx = " << opts.zmx << endl;
				return 1;
			}
			if (getLastInt(opts.zmt)<0) {
				ces << "getLastInt error with zmt = " << opts.zmt << endl;
				if ((opts.printChoice).compare("gui")==0)
					cerr << "getLastInt error with zmt = " << opts.zmt << endl;
				return 1;
			}
			slicesX = getLastInt(opts.zmx);
			slicesT = getLastInt(opts.zmt);
			posCe = j*p.Nb+p.Nb-slicesT; //position C for Euclidean vector, i.e. for negVec
			map<char,unsigned int> posMap;
			posMap['A'] = j*p.NT;
			posMap['B'] = j*p.NT + p.Na-1;
			posMap['C'] = j*p.NT+p.Na+p.Nb-slicesT;
			posMap['D'] = j*p.NT+p.NT-slicesT-1; // extra -1 (so actually D-1) as chiT orthogonal to forward time derivative
			if ((opts.zmt).size()<3) {
				ces << "zmt lacks info, zmt = " << opts.zmt << endl;
				if ((opts.printChoice).compare("gui")==0)
					cerr << "zmt lacks info, zmt = " << opts.zmt << endl;
			}
			for (uint l=0;l<((opts.zmt).size()-2);l++) {
				if (posMap.find(opts.zmt[1+l])!=posMap.end()) {
					posT = posMap.at(opts.zmt[1+l]);
					for (uint k=0;k<slicesT;k++) {
						if (opts.zmt[0]=='n' && p.Pot!=3) 	chiT(posT+k) = negVec(2*(posCe+k));
						else if (opts.zmt[0]=='n' && p.Pot==3) {
												double r = p.r0 + j*p.a;
												chiT(posT+k) = negVec(j)*r;
						}
						else if (opts.zmt[0]=='d')	chiT(posT+k) = f(2*(posT+k+1))-f(2*(posT+k));
						else {
							ces << "choice of zmt(" << opts.zmt << ") not allowed" << endl;
							if ((opts.printChoice).compare("gui")==0)
								cerr << "choice of zmt(" << opts.zmt << ") not allowed" << endl;
							return 1;
						}
					}
				}
			}
			posMap.erase('C');
			posMap.erase('D');
			posMap['C'] = j*p.NT+p.Na+p.Nb-slicesX;
			posMap['D'] = j*p.NT+p.NT-slicesX;
			posCe = j*p.Nb+p.Nb-slicesX;
			if ((opts.zmx).size()<3) {
				ces << "zmx lacks info, zmx = " << opts.zmx << endl;
				if ((opts.printChoice).compare("gui")==0)
					cerr << "zmx lacks info, zmx = " << opts.zmx << endl;
			}
			for (uint l=0;l<((opts.zmx).size()-2);l++) {
				if (posMap.find(opts.zmx[1+l])!=posMap.end()) {
					posX = posMap.at(opts.zmx[1+l]);
					for (uint k=0;k<slicesX;k++) {
						if (opts.zmx[0]=='n' && p.Pot!=3)		chiX(posX+k) = negVec(2*(posCe+k));
						else if (opts.zmx[0]=='n' && p.Pot==3) {
															double r = p.r0 + j*p.a;
															chiX(posX+k) = negVec(j)*r;
						}
						else if (opts.zmx[0]=='d' && p.Pot!=3)
							chiX(posX+k) = f(2*neigh(posX+k,1,1,ps))-f(2*neigh(posX+k,1,-1,ps));
						else {
							ces << "choice of zmx(" << opts.zmx << ") not allowed" << endl;
							if ((opts.printChoice).compare("gui")==0)
								cerr << "choice of zmx(" << opts.zmx << ") not allowed" << endl;
							return 1;
						}
					}
				}
			}
		}
		double normX = chiX.norm();
		double normT = chiT.norm();
		normT = pow(normT,0.5);
		if (abs(normX)<MIN_NUMBER || abs(normT)<MIN_NUMBER) {
			ces << "norm of chiX = " << normX << ", norm of chiT = " << normT << endl;
			if ((opts.printChoice).compare("gui")==0)
				cerr << "norm of chiX = " << normX << ", norm of chiT = " << normT << endl;
		}
		chiX = chiX/normX;
		chiT = chiT/normT;
		
		// allocating memory for DS, DDS
		minusDS = Eigen::VectorXd::Zero(2*p.N*p.NT+2); //initializing to zero
		DDS.setZero(); //just making sure
		Eigen::VectorXi DDS_to_reserve(2*p.N*p.NT+2);//number of non-zero elements per column
		DDS_to_reserve = Eigen::VectorXi::Constant(2*p.N*p.NT+2,13);
		if (abs(p.Theta)<2.0e-16) {
			DDS_to_reserve(0) = p.N+3;
			DDS_to_reserve(1) = 11;
		}
		else {
			DDS_to_reserve(0) = p.N+11;
			DDS_to_reserve(1) = p.N+11;
		}
		DDS_to_reserve(2*p.N*p.NT-2) = 4;
		DDS_to_reserve(2*p.N*p.NT-1) = 4;
		DDS_to_reserve(2*p.N*p.NT) = p.N*((opts.zmx).size()-2);
		DDS_to_reserve(2*p.N*p.NT+1) = 2*p.N*((opts.zmt).size()-2);
		DDS.reserve(DDS_to_reserve);
		
		//initializing to zero
		comp kineticS = 0.0;
		comp kineticT = 0.0;
		comp potV = 0.0;
		comp pot_r = 0.0;
		bound = 0.0;
		erg = Eigen::VectorXcd::Constant(p.NT,-ergZero);
		linErg = Eigen::VectorXcd::Zero(p.NT);
		linErgOffShell = Eigen::VectorXcd::Zero(p.NT);
		linNum = Eigen::VectorXcd::Zero(p.NT);
		linNumOffShell = Eigen::VectorXcd::Zero(p.NT);
		derivErg = Eigen::VectorXcd::Zero(p.NT);
		potErg = Eigen::VectorXcd::Constant(p.NT,-ergZero);
		linErgContm = 0.0;
		linNumContm = 0.0;
		boundRe = Eigen::VectorXd::Zero(p.N);
		boundIm = Eigen::VectorXd::Zero(p.N);
		
		//testing that the potential term is working for pot3
		if (p.Pot==3 && trivialChecks) {
			comp Vtrial = 0.0, Vcontrol = 0.0;
			for (unsigned int j=0; j<p.N; j++) {
				double r = p.r0 + j*p.a;
				paramsV.epsi = r;
				V.setParams(paramsV);
				Vcontrol += pow(f(2*j*p.Nb),2.0)/2.0 - pow(f(2*j*p.Nb),4.0)/4.0/pow(r,2.0);
				Vtrial += V(f(2*j*p.Nb));
			}
			double potTest = pow(pow(real(Vcontrol-Vtrial),2.0) + pow(imag(Vcontrol-Vtrial),2.0),0.5);
			fprintf(cof,"potTest = %8.4g\n",potTest);
			if ((opts.printChoice).compare("gui")==0)
				printf("potTest = %8.4g\n",potTest);
		}

/*----------------------------------------------------------------------------------------------------------------------------
9. assigning minusDS, DDS etc
	- beginning loop over lattice points
	- fixing zero modes
	- linErg, linNum
	- t=(NT-1)
	- t=0
	- bulk
	- extras (*4.0*pi, etc)
	- E, N, bound, W
----------------------------------------------------------------------------------------------------------------------------*/
		
		// beginning loop over lattice points
		for (lint j = 0; j < p.N*p.NT; j++) {		
			uint t 					= intCoord(j,0,ps); //coordinates
			uint x 					= intCoord(j,1,ps);
			int neighPosX 			= neigh(j,1,1,ps);
			int neighNegX  			= neigh(j,1,-1,ps);
			
			comp Dt 		= DtFn(t,ps);
			comp dt 		= dtFn(t,ps);
			comp dtm 		= (t>0? dtFn(t-1,ps): p.b);
			double Dx 		= DxFn(x,ps);
			double dx 		= dxFn(x,ps);
			double dxm 		= (x>0? dxFn(x-1,ps): p.a);
			
			if (p.Pot==3) {
				paramsV.epsi = p.r0+x*p.a;
				V.setParams(paramsV);
				dV.setParams(paramsV);
				ddV.setParams(paramsV);
			}
		
			if (abs(chiX(j))>MIN_NUMBER && p.Pot!=3) { //spatial zero mode lagrange constraint
				DDS.insert(2*j,2*p.N*p.NT) 		= Dx*chiX(j); 
				DDS.insert(2*p.N*p.NT,2*j) 		= Dx*chiX(j);
				minusDS(2*j) 						+= -Dx*chiX(j)*f(2*p.N*p.NT);
				minusDS(2*p.N*p.NT) 				+= -Dx*chiX(j)*f(2*j);
			}
				
			if (abs(chiT(j))>MIN_NUMBER && t<(p.NT-1) && (p.Pot!=3 || (x>0 && x<(p.N-1)))) {
				DDS.coeffRef(2*(j+1),2*p.N*p.NT+1) 	+= Dx*chiT(j); //chiT should be 0 at t=(p.NT-1) or this line will go wrong
				DDS.coeffRef(2*p.N*p.NT+1,2*(j+1)) 	+= Dx*chiT(j);
				DDS.coeffRef(2*j,2*p.N*p.NT+1) 		+= -Dx*chiT(j);
				DDS.coeffRef(2*p.N*p.NT+1,2*j) 		+= -Dx*chiT(j);
	            minusDS(2*(j+1)) 						+= -Dx*chiT(j)*f(2*p.N*p.NT+1);
	            minusDS(2*j) 							+= Dx*chiT(j)*f(2*p.N*p.NT+1);
	            minusDS(2*p.N*p.NT+1) 				+= -Dx*chiT(j)*(f(2*(j+1))-f(2*j));
			}
				
			if (t<(p.NT-1)) {
				for (uint k=0;k<p.N;k++) {
					lint l = k*p.NT+t;
					linErgOffShell(t) += 0.5*( omega_2(x,k)*( Cf(l)-p.minima[0] )*( Cf(j)-p.minima[0] ) \
									+ omega_0(x,k)*( Cf(l+1)-Cf(l) )*( Cf(j+1)-Cf(j) )/pow(dt,2.0));
					linNumOffShell(t) += 0.5*(omega_1(x,k)*( Cf(l)-p.minima[0] )*( Cf(j)-p.minima[0] ) \
									+ omega_m1(x,k)*( Cf(l+1)-Cf(l) )*( Cf(j+1)-Cf(j) )/pow(dt,2.0));
				}
			}
				
			if (abs(p.Theta)>MIN_NUMBER) {
				for (uint k=0;k<p.N;k++) {
					lint l = k*p.NT+t;
					linNum(t) += 2.0*p.Gamma*omega_1(x,k)*( (f(2*l)-p.minima[0])*(f(2*j)-p.minima[0])/pow(1.0+p.Gamma,2.0)\
							 + f(2*j+1)*f(2*l+1)/pow(1.0-p.Gamma,2.0) );
					linErg(t) += 2.0*p.Gamma*omega_2(x,k)*( (f(2*l)-p.minima[0])*(f(2*j)-p.minima[0])/pow(1.0+p.Gamma,2.0)\
							+ f(2*j+1)*f(2*l+1)/pow(1.0-p.Gamma,2.0) );
				}
			}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//boundaries
			if (p.Pot==3 && x==(p.N-1)) {
				DDS.insert(2*j,2*j) 	= 1.0; // f=0 at r=R
				DDS.insert(2*j+1,2*j+1) = 1.0;
			}
			else if (p.Pot==3 && x==0) {
				kineticS 				+= 	Dt*pow(Cf(neighPosX),2.0)/dx/2.0;
				derivErg(t) 			+= 	pow(Cf(neighPosX),2.0)/dx/2.0; //n.b. the Dt/dt difference is ignored for erg(t)
				erg(t) 					+= 	pow(Cf(neighPosX),2.0)/dx/2.0;
				
				DDS.insert(2*j,2*j) 	= 1.0; // f=0 at r=0
				DDS.insert(2*j+1,2*j+1) = 1.0;
			}			
			else if (t==(p.NT-1)) {
				if (neighPosX!=-1) {
					kineticS			+= Dt*pow(Cf(neighPosX)-Cf(j),2.0)/dx/2.0;
					derivErg(t) 		+= pow(Cf(neighPosX)-Cf(j),2.0)/dx/2.0;
					erg(t) 				+= pow(Cf(neighPosX)-Cf(j),2.0)/dx/2.0;
				}			
				potV 					+= Dt*Dx*V(Cf(j));
				pot_r 					+= Dt*Dx*Vr(Cf(j));
				erg(t) 					+= Dx*V(Cf(j)) + Dx*Vr(Cf(j));
				potErg(t) 				+= Dx*V(Cf(j)) + Dx*Vr(Cf(j));
			
				DDS.insert(2*j,2*(j-1)+1) = 1.0; //zero imaginary part of time derivative
				DDS.insert(2*j+1,2*j+1)   = 1.0; //zero imaginary part
			}
			else if (t==0) {
				kineticT 	+= Dx*pow(Cf(j+1)-Cf(j),2.0)/dt/2.0;
				derivErg(t) += Dx*pow(Cf(j+1)-Cf(j),2.0)/pow(dt,2.0)/2.0;
				erg(t) 		+= Dx*pow(Cf(j+1)-Cf(j),2.0)/pow(dt,2.0)/2.0;
				
				///////////////////////////////// including other terms in action at t=0 ///////////////////////////
				if (neighPosX!=-1) {
					kineticS 	+= Dt*pow(Cf(neighPosX)-Cf(j),2.0)/dx/2.0;
					derivErg(t) += pow(Cf(neighPosX)-Cf(j),2.0)/dx/2.0;
					erg(t) 		+= pow(Cf(neighPosX)-Cf(j),2.0)/dx/2.0;
				}
				potV 		+= Dt*Dx*V(Cf(j));
				pot_r 		+= Dt*Dx*Vr(Cf(j));
				potErg(t) 	+= Dx*V(Cf(j)) + Dx*Vr(Cf(j));
				erg(t) 		+= Dx*V(Cf(j)) + Dx*Vr(Cf(j));
				double tnt = 1.0;
				if ((opts.bds).compare("uc")!=0)
					tnt *= p.Theta;
				//////////////////////////////////////equation I - both///////////////////////////////////
				for (uint k=1; k<2*2; k++) {
			        int sign = pow(-1,k+1);
			        uint direc = (uint)(k/2.0);
			        int neighb = neigh(j,direc,sign,ps);
			        double dxd = (sign==1? dx: dxm);
			        if (direc == 0) {
			            minusDS(2*j+1) 					+= Dx*imag(Cf(j+sign)/dt);
			            DDS.coeffRef(2*j+1,2*(j+sign)) 	+= -imag(Dx/dt);
			            DDS.coeffRef(2*j+1,2*(j+sign)+1)+= -real(Dx/dt);
			           }
			        else if (neighb!=-1) {
			            minusDS(2*j+1) 					+= -imag(Dt*Cf(neighb))/dxd;
			            DDS.coeffRef(2*j+1,2*neighb) 	+= imag(Dt)/dxd;
			            DDS.coeffRef(2*j+1,2*neighb+1) 	+= real(Dt)/dxd;
			        }
			    }
			    comp temp0 = Dx/dt - Dt*(1.0/dx+1.0/dxm);
			    if (neighPosX==-1) 		temp0 += Dt/dx;
			    else if (neighNegX==-1) temp0 += Dt/dxm;
			    comp temp1 = Dt*Dx*( dV(Cf(j)) + dVr(Cf(j)) );//dV terms should be small
			    comp temp2 = Dt*Dx*(ddV(Cf(j)) + ddVr(Cf(j)));
			    
			    minusDS(2*j+1) 				+= imag(-temp0*Cf(j) + temp1 );
			    DDS.coeffRef(2*j+1,2*j) 	+= imag(temp0 - temp2 );
			    DDS.coeffRef(2*j+1,2*j+1) 	+= real(temp0 - temp2 );
				/////////////////////////////////////////////////////////////////////////////////////////
				if (abs(p.Theta)<MIN_NUMBER) {
					//simplest boundary conditions replaced by ones continuously connected to theta!=0 ones
					//DDS.insert(2*j+1,2*(j+1)+1) = 1.0; //zero imaginary part of time derivative
					//DDS.insert(2*j,2*j+1) = 1.0; //zero imaginary part
				
					/////////////////////////////////////equation R - theta=0//////////////////////////////////////
					for (uint k=0;k<p.N;k++) {
						if (abs(omega_1(x,k))>MIN_NUMBER) {
							lint m=k*p.NT;
							DDS.coeffRef(2*j,2*m+1) += -2.0*omega_1(x,k);
							minusDS(2*j) 			+= 2.0*omega_1(x,k)*f(2*m+1);
						}
					}
					////////////////////////////////////////////////////////////////////////////////////////
				}
				else {
					for (uint k=0;k<p.N;k++) {
						if (abs(omega_1(x,k))>MIN_NUMBER) {
							/////////////////////equation I - theta!=0//////////////
							lint m=k*p.NT;
							DDS.coeffRef(2*j+1,2*m) += (1.0-p.Gamma)*omega_1(x,k)/(1.0+p.Gamma);
							minusDS(2*j+1) 			+= -(1.0-p.Gamma)*omega_1(x,k)*(f(2*m)-p.minima[0])/(1.0+p.Gamma);
							/////////////////////equation R - theta!=0//////////////
							minusDS(2*j) 			+= tnt*f(2*m+1)*omega_1(x,k)*(1+p.Gamma)/(1-p.Gamma);
							DDS.coeffRef(2*j,2*m+1)	+= -tnt*omega_1(x,k)*(1.0+p.Gamma)/(1.0-p.Gamma);
							bound 				+= -(1.0-p.Gamma)*omega_1(x,k)*(f(2*j)-p.minima[0])\
														*(f(2*m)-p.minima[0])/(1.0+p.Gamma)
														 + (1.0+p.Gamma)*omega_1(x,k)*f(2*j+1)*f(2*m+1)/(1.0-p.Gamma);
						}
					}
					//////////////////////////////////////equation R - theta!=0//////////////////////////////
					for (uint k=1; k<2*2; k++){
				        int sign = pow(-1,k+1);
				        uint direc = (uint)(k/2.0);
				        int neighb = neigh(j,direc,sign,ps);
				        double dxd = (sign==1? dx: dxm);
				        if (direc == 0) {
				            minusDS(2*j) 					+= tnt*Dx*real(Cf(j+sign)/dt);
			            	DDS.coeffRef(2*j,2*(j+sign)) 	+= -tnt*real(Dx/dt);
			            	DDS.coeffRef(2*j,2*(j+sign)+1) 	+= tnt*imag(Dx/dt);
				        }
				        else if (neighb!=-1) {
				            minusDS(2*j) 					+= -tnt*real(Dt*Cf(neighb))/dxd;
				            DDS.coeffRef(2*j,2*neighb) 		+= tnt*real(Dt)/dxd;
			            	DDS.coeffRef(2*j,2*neighb+1) 	+= -tnt*imag(Dt)/dxd;
				        }
				    }
				    minusDS(2*j) 			+= tnt*real(-temp0*Cf(j) + temp1 );
			    	DDS.coeffRef(2*j,2*j) 	+= tnt*real(temp0 - temp2 );
			    	DDS.coeffRef(2*j,2*j+1) += tnt*imag(-temp0 + temp2 );
				}
			}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//bulk
			else {
				if (neighPosX!=-1) {
					kineticS 	+= Dt*pow(Cf(neighPosX)-Cf(j),2.0)/dx/2.0;
					erg(t) 	 	+= pow(Cf(neighPosX)-Cf(j),2.0)/dx/2.0;
					derivErg(t) += pow(Cf(neighPosX)-Cf(j),2.0)/dx/2.0;
				}
				
				kineticT 	+= Dx*pow(Cf(j+1)-Cf(j),2.0)/dt/2.0;
				potV 		+= Dt*Dx*V(Cf(j));
				pot_r 		+= Dt*Dx*Vr(Cf(j));
				erg(t) 		+= Dx*pow(Cf(j+1)-Cf(j),2.0)/pow(dt,2.0)/2.0 + Dx*V(Cf(j)) + Dx*Vr(Cf(j));
				derivErg(t) += Dx*pow(Cf(j+1)-Cf(j),2.0)/pow(dt,2.0)/2.0;
				potErg(t) 	+= Dx*V(Cf(j)) + Dx*Vr(Cf(j));
			
	            for (uint k=0; k<2*2; k++) {
	                int sign = pow(-1,k);
	                uint direc = (uint)(k/2.0);
	                int neighb = neigh(j,direc,sign,ps);
	                comp dtd = (sign==1? dt: dtm);
	                double dxd = (sign==1? dx: dxm);
	                if (direc == 0) {
	                    minusDS(2*j) 					+= real(Dx*Cf(j+sign)/dtd);
	                    minusDS(2*j+1) 					+= imag(Dx*Cf(j+sign)/dtd);
	                    DDS.insert(2*j,2*(j+sign)) 		= -real(Dx/dtd);
	                    DDS.insert(2*j,2*(j+sign)+1) 	= imag(Dx/dtd);
	                    DDS.insert(2*j+1,2*(j+sign)) 	= -imag(Dx/dtd);
	                    DDS.insert(2*j+1,2*(j+sign)+1) 	= -real(Dx/dtd);
	                }
	                else if (neighb!=-1) {                        
	                    minusDS(2*j) 					+= -real(Dt*Cf(neighb)/dxd);
	                    minusDS(2*j+1) 					+= -imag(Dt*Cf(neighb)/dxd);
	                    DDS.insert(2*j,2*neighb) 		= real(Dt/dxd);
	                    DDS.insert(2*j,2*neighb+1) 		= -imag(Dt/dxd);
	                    DDS.insert(2*j+1,2*neighb) 		= imag(Dt/dxd);
	                    DDS.insert(2*j+1,2*neighb+1) 	= real(Dt/dxd);
	                }
            	}
	            comp temp0 = Dx*(1.0/dt + 1.0/dtm) - Dt*(1.0/dx + 1.0/dxm);
	            if (neighPosX==-1) 		temp0 += Dt/dx;
	            else if (neighNegX==-1) temp0 += Dt/dxm;
	            comp temp1 = Dt*Dx*(dV(Cf(j)) + dVr(Cf(j)));
	            comp temp2 = Dt*Dx*(ddV(Cf(j)) + ddVr(Cf(j)));
	                
	            minusDS(2*j) 			+= real(temp1 - temp0*Cf(j));
	            minusDS(2*j+1) 			+= imag(temp1 - temp0*Cf(j));
	            DDS.insert(2*j,2*j) 	= real(-temp2 + temp0);
	            DDS.insert(2*j,2*j+1) 	= imag(temp2 - temp0);
	            DDS.insert(2*j+1,2*j) 	= imag(-temp2 + temp0);
	            DDS.insert(2*j+1,2*j+1) = real(-temp2 + temp0);
	        }
	    } // end of loop over j
	    
	    if (p.Pot==3) 		DDS.insert(2*p.N*p.NT,2*p.N*p.NT) = 1.0;
	    action = kineticT - kineticS - potV - pot_r;
	    linErgOffShell(p.NT-1) = linErgOffShell(p.NT-2);
    	linNumOffShell(p.NT-1) = linNumOffShell(p.NT-2);
    	
	    if (p.Pot==3) {
	    	action 			*= 4.0*pi;
	    	derivErg 		*= 4.0*pi;
	    	potErg 			*= 4.0*pi;
	    	erg				*= 4.0*pi;
	    	linErg			*= 4.0*pi;
	    	linNum			*= 4.0*pi;
	    	linNumOffShell 	*= 4.0*pi;
	    	linErgOffShell 	*= 4.0*pi;
	    }
	   
	    if (abs(p.Theta)<MIN_NUMBER) {
	    	linErg = linErgOffShell;
	    	linNum = linNumOffShell;
	    }
	    
	    for (uint x=1; x<(p.N-1); x++) {
	    	lint m = x*p.NT;
	    	boundRe(x) = minusDS(2*m);
	    	boundIm(x) = minusDS(2*m+1);
	    }
	    
	    //defining E, Num and W
		E = real(linErg(0));
		Num = real(linNum(0));
		W = - E*2.0*p.Tb - p.Theta*Num - bound + 2.0*imag(action);
		
/*----------------------------------------------------------------------------------------------------------------------------
10. checks
	- erg = potErg + derivErg (trivial)
	- checkContm
	- checkAB
	- checkABNE
	- checkLR (linear representation of phi)
	- checkReg
	- checkLin
	- checkOS
	- checkCon
	- checkIE
	- checkTrue
	- checkLatt
	- checkBound
	- checkChiT
	
----------------------------------------------------------------------------------------------------------------------------*/
		
		//trivial test that erg=potErg+derivErg
		if (trivialChecks) {
			for (uint j=0; j<p.NT; j++) {
				double diff = absDiff(erg(j),potErg(j)+derivErg(j));
				if (diff>1.0e-14) {
					ces << "erg(" << j << ") != potErg + derivErg. absDiff = " << diff << endl;
					if ((opts.printChoice).compare("gui")==0)
						cerr << "erg(" << j << ") != potErg + derivErg. absDiff = " << diff << endl;
				}
			}
		}
		
		//calculating continuum approx to linErg and linNum on initial time slice - redundant
		if (p.Pot==3 && abs(p.Theta)<MIN_NUMBER) {
			for (uint k=1; k<p.N; k++) {
				double momtm = k*pi/(p.L-p.r0);
				double freqSqrd = 1.0+pow(momtm,2.0);
				double Asqrd, integral1 = 0.0, integral2 = 0.0;
				for (unsigned int l=0; l<p.N; l++) {
					double r = p.r0 + l*p.a;
					lint m = l*p.NT;
					integral1 += p.a*f(2*m)*pow(2.0/p.L,0.5)*sin(momtm*r);
					integral2 += p.a*(f(2*(m+1))-f(2*m))*pow(2.0/p.L,0.5)*sin(momtm*r)/p.b;
				}
				Asqrd = pow(integral1,2.0) + pow(integral2,2.0)/freqSqrd;
				linErgContm += 2.0*pi*Asqrd*freqSqrd;
				linNumContm += 2.0*pi*Asqrd*pow(freqSqrd,0.5);
			}
		}
		double contmErgTest = absDiff(E,linErgContm);
		double contmNumTest = absDiff(Num,linNumContm);
		contmErgTest>contmNumTest? checkContm.add(contmErgTest): checkContm.add(contmNumTest);
		
		//calculating a_k and b*_k
		cVec a_k(p.N), b_k(p.N); //p.N.b. b_k means b*_k but can't put * in the name
		a_k = Eigen::VectorXcd::Zero(p.N);
		b_k = Eigen::VectorXcd::Zero(p.N);
		double T0 = 0.0;
		comp dt0 = dtFn(0,ps);
		for (uint n=0; n<p.N; n++) {
			double w_n = freqs(n);
			double w_n_e = freqs_exp(n);
			for (uint j=0; j<p.N; j++) {
				lint m=j*p.NT;
				double sqrtDj = sqrt(4.0*pi*DxFn(j,ps));
				if (abs(w_n)>1.0e-16 && abs(w_n_e)>1.0e-16) {
					a_k(n) += exp(ii*w_n_e*T0)*sqrt(2.0*w_n)*modes(j,n)* \
								sqrtDj*((Cf(m+1)-p.minima[0])-(Cf(m)-p.minima[0])*exp(ii*w_n_e*dt0)) \
									/(exp(-ii*w_n_e*dt0)-exp(ii*w_n_e*dt0));
					b_k(n) += exp(-ii*w_n_e*T0)*sqrt(2.0*w_n)*modes(j,n)* \
								sqrtDj*((Cf(m+1)-p.minima[0])-(Cf(m)-p.minima[0])*exp(-ii*w_n_e*dt0)) \
									/(exp(ii*w_n_e*dt0)-exp(-ii*w_n_e*dt0));
				}
			}
		}			
		
		//using a_k and b*_k to check that a_k=Gamma*b*_k and f->f_lin as t->0 and that Sum_k(a_k*b*_k)=linNum
		// and that Sum_k(w_k*a_k*b*_k)=linErg
		double ABtest = 0.0, linRepTest;
		linNumAB = 0.0, linErgAB = 0.0;
		cVec linRep(p.N), p0(p.N);
		linRep = Eigen::VectorXcd::Constant(p.N,p.minima[0]);
		for (uint j=0; j<p.N; j++) {
			ABtest += absDiff(a_k(j),conj(b_k(j))*p.Gamma);
			p0(j) = Cf(j*p.NT);
			linNumAB += a_k(j)*b_k(j);
			linErgAB += freqs(j)*a_k(j)*b_k(j);
			if (trivialChecks) {
				for (uint n=0; n<p.N; n++) {
					double w_n = freqs(n);
					double w_n_e = freqs_exp(n);
					double sqrtDj = sqrt(4.0*pi*DxFn(j,ps));
					if (abs(w_n)>1.0e-16) {
						linRep(j) += modes(j,n)*(a_k(n)*exp(-ii*w_n_e*T0)+b_k(n)*exp(ii*w_n_e*T0)) \
										/sqrt(2.0*w_n)/sqrtDj;
					}
				}
			}
		}
		if (trivialChecks) {
			linRepTest = absDiff(p0,linRep);
			checkLR.add(linRepTest);
		}
		ABtest /= (double)p.N;
		checkAB.add(ABtest);
		double ABNtest = absDiff(linNumAB,linNum(0));
		double ABEtest = absDiff(linErgAB,linErg(0));
		ABNtest>ABEtest? checkABNE.add(ABNtest): checkABNE.add(ABEtest);
		
		//checking pot_r is much smaller than the other potential terms
		checkReg.add(abs(pot_r/potV));
		checkReg.checkMessage(ces);
		if ((opts.printChoice).compare("gui")==0)
			checkReg.checkMessage(cerr);
					
		//checking linearisation of linErg and linNum
		double linTestE;		double linTestN;
		double linEMax = 0.0;	double linEMin = 5.0e15; //surely it's going to be less than this
		double linNMax = 0.0;	double linNMin = 5.0e15;
		uint linearInt = (uint)(p.Na/10);
		for (uint j=1;j<(linearInt+1);j++) {
			if (abs(linErgOffShell(j))>linEMax) linEMax = abs(linErgOffShell(j));
			if (abs(linErgOffShell(j))<linEMin) linEMin = abs(linErgOffShell(j));
			if (abs(linNumOffShell(j))>linNMax) linNMax = abs(linNumOffShell(j));
			if (abs(linNumOffShell(j))<linNMin) linNMin = abs(linNumOffShell(j));
		}
		linTestE = absDiff(linEMax,linEMin);
		linTestN = absDiff(linNMax,linNMin);
		linTestE>linTestN? checkLin.add(linTestE): checkLin.add(linTestN);
		
		//checking agreement of on and off shell linear energy at initial time
		double onShellTest = absDiff(linErg(0),linErgOffShell(0));
		checkOS.add(onShellTest);
		
		//checking conservation of E
		double conservTest = absDiff(erg(1),erg(p.NT-2));
		checkCon.add(conservTest);
		
		//testing imaginary part of energy
		double imErgTest = 0.0;
		for (uint j=0; j<p.NT; j++) if (abs(erg(j))>MIN_NUMBER) imErgTest += imag(erg(j))/abs(erg(j));
		imErgTest /= (double)p.NT;
		checkIE.add(imErgTest);
		
		//checking agreement between erg and linErg
		E_exact = 0.0;
		for (uint j=1; j<(linearInt+1); j++) E_exact += real(erg(j));
		E_exact /= (double)linearInt;
		double trueTest = absDiff(E,E_exact);
		checkTrue.add(trueTest);
		if (!isfinite(trueTest)) {
			ces << "E = " << E << ", E_exact = " << E_exact << ", linearInt = " << linearInt << endl;
			if ((opts.printChoice).compare("gui")==0)
				cerr << "E = " << E << ", E_exact = " << E_exact << ", linearInt = " << linearInt << endl;
		}
		
		//checking lattice small enough for E, should have parameter for this
		double momTest = E*p.b/Num/pi; //perhaps should have a not b here
		checkLatt.add(momTest);
		
		//checking initial boundary conditions satisfied
		double normP = f.norm();
		double boundReTest = boundRe.norm()*(2*p.N*p.NT+2)/normP/(p.N-2.0);
		double boundImTest = boundIm.norm()*(2*p.N*p.NT+2)/normP/(p.N-2.0);
		checkBoundRe.add(boundReTest);
		checkBoundIm.add(boundImTest);
		
		// checking chiT orthogonality satisfied
		double chiTTest = minusDS(2*p.N*p.NT+1)*(2*p.N*p.NT+2)/normP;
		checkChiT.add(chiTTest);

/*----------------------------------------------------------------------------------------------------------------------------
11. printing early 1
	- filenames and saveoptions
	- phi
	- minusDS
	- DDS
	- chiX, chiT
	- energies
	- a_k, b_k
	
----------------------------------------------------------------------------------------------------------------------------*/

	//printing early if desired
	if ((opts.printChoice).compare("n")!=0) {
		so_tp.printType = SaveOptions::ascii;
		so_simple.printType = SaveOptions::ascii;
		Filename basic = (string)("data/"+timenumber+"basic_fLoop_"+numberToString<uint>(fileLoop)\
							+"_loop_"+numberToString<uint>(loop)+"_run_"+numberToString<uint>(runs_count)+".dat");
		if ((opts.printChoice).compare("v")==0 || (opts.printChoice).compare("e")==0) {
			Filename vEFile = basic;
			vEFile.ID = "mainminusDSE";
			save(vEFile,so_tp,minusDS);
		}
		if ((opts.printChoice).compare("m")==0 || (opts.printChoice).compare("e")==0) {
			Filename mEFile = basic;
			mEFile.ID = "mainDDSE";
			save(mEFile,so_simple,DDS);
		}
		if ((opts.printChoice).compare("z")==0 || (opts.printChoice).compare("e")==0) {
			so_tp.vectorType = SaveOptions::real;			
			Filename zEFile = basic;
			zEFile.ID = "mainchiTE";
			save(zEFile,so_tp,chiT);
			zEFile.ID = "mainchiXE";
			save(zEFile,so_tp,chiX);
			so_tp.vectorType = SaveOptions::complex;
		}
		if ((opts.printChoice).compare("l")==0 || (opts.printChoice).compare("e")==0) {
			Filename lEFile = basic;
			lEFile.ID = "mainlinErgE";
			save(lEFile,so_simple,linErg);
			lEFile.ID = "mainergE";
			save(lEFile,so_simple,erg);
		}
		if ((opts.printChoice).compare("ab")==0 || (opts.printChoice).compare("e")==0) {
			Filename abEFile = basic;
			abEFile.ID = "mainakE";
			save(abEFile,so_simple,a_k);
			abEFile.ID = "mainbkE";
			save(abEFile,so_simple,b_k);
		}
		so_tp.printType = SaveOptions::binary;
		so_simple.printType = SaveOptions::binary;
	}
	
	
/*----------------------------------------------------------------------------------------------------------------------------
12. solving for delta	
	- defining delta, solver etc
	- analysing pattern
	- factorising
	- solving
	- checking matrix inversion
	- f' = f + delta
	- Cf'
----------------------------------------------------------------------------------------------------------------------------*/
		
		//solving for delta in DDS*delta=minusDS, where f' = f + delta		
		vec delta(2*p.N*p.NT+2);
		delta = Eigen::VectorXd::Zero(2*p.N*p.NT+2);
		DDS.prune(MIN_NUMBER);
		DDS.makeCompressed();
		Eigen::SparseLU<spMat> solver;
	
		/*solver.analyzePattern(DDS);
		if(solver.info()!=Eigen::Success) {
			ces << "DDS pattern analysis failed, solver.info() = "<< solver.info() << endl;
			if ((opts.printChoice).compare("gui")==0) {
				cerr << "DDS pattern analysis failed, solver.info() = "<< solver.info() << endl;
				printErrorInformation(f,"f",2);
				cout << endl;
				printErrorInformation(minusDS,"mds",2);
				cout << endl;
				printErrorInformation(DDS,"DDS");
				cout << endl;
			}
			return 1;
		}		
		solver.factorize(DDS);
		if(solver.info()!=Eigen::Success) {
			ces << "Factorization failed, solver.info() = "<< solver.info() << endl;
			if ((opts.printChoice).compare("gui")==0) {
				cerr << "Factorization failed, solver.info() = "<< solver.info() << endl;
				printErrorInformation(f,"f",2);
				cout << endl;
				printErrorInformation(minusDS,"mds",2);
				cout << endl;
				printErrorInformation(DDS,"DDS");
				cout << endl;
			}
			return 1;
		}*/
		solver.compute(DDS);
		if(solver.info()!=Eigen::Success) {
			ces << "Compute failed, solver.info() = "<< solver.info() << endl;
			if ((opts.printChoice).compare("gui")==0) {
				cerr << "Compute failed, solver.info() = "<< solver.info() << endl;
				printErrorInformation(f,"f",2);
				cout << endl;
				printErrorInformation(minusDS,"mds",2);
				cout << endl;
				printErrorInformation(DDS,"DDS");
				cout << endl;
			}
			return 1;
		}
		delta = solver.solve(minusDS);// use the factorization to solve for the given right hand side
		if(solver.info()!=Eigen::Success) {
			ces << "Solving failed, solver.info() = "<< solver.info() << endl;
			ces << "log(abs(det(DDS))) = " << solver.logAbsDeterminant() << endl;
			ces << "sign(det(DDS)) = " << solver.signDeterminant() << endl;
			if ((opts.printChoice).compare("gui")==0) {
				cerr << "Solving failed, solver.info() = "<< solver.info() << endl;
				cerr << "log(abs(det(DDS))) = " << solver.logAbsDeterminant() << endl;
				cerr << "sign(det(DDS)) = " << solver.signDeterminant() << endl;
				printErrorInformation(f,"f",2);
				cout << endl;
				printErrorInformation(minusDS,"mds",2);
				cout << endl;
				printErrorInformation(delta,"delta",2);
				cout << endl;
				printErrorInformation(DDS,"DDS");
				cout << endl;
			}
			return 1;
		}
	
		//independent check on whether calculation worked
		vec diff(2*p.N*p.NT+2);
		diff = DDS*delta-minusDS;
		double maxDiff = diff.maxCoeff();
		maxDiff = abs(maxDiff);
		checkInv.add(maxDiff);
		checkInv.checkMessage(ces);
		if ((opts.printChoice).compare("gui")==0)
			checkInv.checkMessage(cerr);
		if (!checkInv.good()) return 1;

		//assigning values to phi
		f += delta;
	
		//passing changes on to complex vector
		Cf = vecComplex(f,p.N*p.NT);
		
/*----------------------------------------------------------------------------------------------------------------------------
13. printing early 2
	- filenames and saveoptions
	- delta
	
----------------------------------------------------------------------------------------------------------------------------*/

	//printing early if desired
	if ((opts.printChoice).compare("n")!=0) {
		so_tp.printType = SaveOptions::ascii;
		Filename basic = (string)("data/"+timenumber+"basic_fLoop_"+numberToString<uint>(fileLoop)\
							+"_loop_"+numberToString<uint>(loop)+"_run_"+numberToString<uint>(runs_count)+".dat");
		if ((opts.printChoice).compare("f")==0 || (opts.printChoice).compare("e")==0) {
			Filename pEFile = basic;
			pEFile.ID = "mainpE";
			save(pEFile,so_tp,f);
		}
		if ((opts.printChoice).compare("d")==0 || (opts.printChoice).compare("e")==0) {
			Filename dEFile = basic;
			dEFile.ID = "maindeltaE";
			save(dEFile,so_tp,delta);
		}
		so_tp.printType = SaveOptions::binary;
	}
	
/*----------------------------------------------------------------------------------------------------------------------------
14. convergence
	- evaluating norms
	- adding convergence checks
	- printing convergence checks
	
----------------------------------------------------------------------------------------------------------------------------*/
		//convergence issues
		//evaluating norms
		double normDS = minusDS.norm();
		double maxDS = minusDS.maxCoeff();
		double minDS = minusDS.minCoeff();
		if (-minDS>maxDS) maxDS = -minDS;
		double maxP = f.maxCoeff();
		double minP = f.minCoeff();
		if (-minP>maxP) maxP = -minP;
		double normDelta = delta.norm();
	
		// adding convergence checks
		checkAction.add(absDiff(action,action_last));
		action_last = action;
		checkSoln.add(normDS/normP);
		checkSolnMax.add(maxDS/maxP);
		checkDelta.add(normDelta/normP);

		// rate of convergence
		double convergRate = (checkDelta.size()>1? \
			log((checkDelta.tests())[checkDelta.size()-1])/log((checkDelta.tests())[checkDelta.size()-2]):0.0);
		
		//printing tests to see convergence
		if (runs_count==1) {
			fprintf(cof,"%5s%5s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s\n","loop","run","sol","solM","delta","converg"\
			,"linear","true erg","on shell","AB","ABNE","conserv","latt");
			if ((opts.printChoice).compare("gui")==0)
				printf("%5s%5s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s\n","loop","run","sol","solM","delta","converg"\
				,"linear","true erg","on shell","AB","ABNE","conserv","latt");
		}
		fprintf(cof,"%5i%5i%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g\n",loop,runs_count\
		,checkSoln.back(),checkSolnMax.back(),checkDelta.back(),convergRate,checkLin.back(),checkTrue.back()\
		,checkOS.back(),checkAB.back(),checkABNE.back(),\
			checkCon.back(),checkLatt.back());
		if ((opts.printChoice).compare("gui")==0)
			printf("%5i%5i%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g\n",loop,runs_count\
			,checkSoln.back(),checkSolnMax.back(),checkDelta.back(),convergRate,checkLin.back(),checkTrue.back()\
			,checkOS.back(),checkAB.back(),checkABNE.back(),checkCon.back(),checkLatt.back());
		
		if (!checkDelta.good()) {
			checkDelta.checkMessage(ces);
			if ((opts.printChoice).compare("gui")==0)
				checkDelta.checkMessage(cerr);
			break;
		}
		
	} //ending while loop
/*----------------------------------------------------------------------------------------------------------------------------
15. printing output
	- check messages
	- stopping clock
	- stepping stepper
	- printing results to terminal
	- printing results to file
	- printing (and plotting if vectors):
		- f
		- linErg
		- erg
		- if printEverything:
			- minusDS
			- DDS
			- linErgOffShell
			- linNum
			- derivErg
			- potErg
	- return
	
----------------------------------------------------------------------------------------------------------------------------*/

	// check messages
	checkTrue.checkMessage(ces);
	checkOS.checkMessage(ces);
	checkABNE.checkMessage(ces);
	if (p.Pot==3 && abs(p.Theta)<MIN_NUMBER) checkContm.checkMessage(ces);
	checkCon.checkMessage(ces);
	checkIE.checkMessage(ces);
	checkLatt.checkMessage(ces);
	checkBoundRe.checkMessage(ces);
	checkBoundIm.checkMessage(ces);
	if ((opts.printChoice).compare("gui")==0) {
		checkTrue.checkMessage(cerr);
		checkOS.checkMessage(cerr);
		checkABNE.checkMessage(cerr);
		if (p.Pot==3 && abs(p.Theta)<MIN_NUMBER) checkContm.checkMessage(cerr);
		checkCon.checkMessage(cerr);
		checkIE.checkMessage(cerr);
		checkLatt.checkMessage(cerr);
		checkBoundRe.checkMessage(cerr);
		checkBoundIm.checkMessage(cerr);
	}
	
	//stopping clock
	time = clock() - time;
	double realtime = time/1000000.0;
	
	// stepping stepper
	if(((opts.loopChoice).substr(0,5)).compare("const")==0) {
		double F = 0.0;
		if ((opts.loopChoice)[(opts.loopChoice).size()-1]=='W') 
			F = W;
		else if ((opts.loopChoice)[(opts.loopChoice).size()-1]=='E')
			F = E;
		else if ((opts.loopChoice)[(opts.loopChoice).size()-1]=='N')
			F = Num;
		else {
			ces << "Stepper error: option " << opts.loopChoice << " not possible" << endl;
			if ((opts.printChoice).compare("gui")==0)
				cerr << "Stepper error: option " << opts.loopChoice << " not possible" << endl;
			return 1;
		}
		double angleToPrint = (loop==0? 0.0: stepper.stepAngle());
		if (absDiff(opts.loopMin,F)<stepper.closeness() && loop==0)
			stepper.addResult(opts.loopMin);
		else if (loop==0) {
			opts.loopMin = F;
			opts.inF = "mn";
			opts.save(optionsFile);
			stepper.addResult(F);
		}
		else
			stepper.addResult(F);
		string keep = (stepper.keep()? "y": "n");
		FILE * stepOs;
		string stepFile = "data/"+timenumber+"mainStep_fLoop_"+numberToString<uint>(fileLoop)+".dat";
		stepOs = fopen(stepFile.c_str(),"a");
		fprintf(stepOs,"%12s%5i%5i%5i%5i%6g%13.5g%13.5g%13.5g%13.5g%13.5g%8s\n",\
					timenumber.c_str(),fileLoop,loop,p.N,p.NT,p.L,p.DE,p.Tb,p.Theta,angleToPrint,F,keep.c_str());
		fclose(stepOs);
		//fprintf(cof,"%12s%30s\n","steps:",stepFile.c_str());
		//if ((opts.printChoice).compare("gui")==0)
		//printf("%12s%30s\n","steps:",stepFile.c_str());
	}
	else
		stepper.addResult(1.0); // choice irrelevant but a value is require to make step

	// printing results to terminal
	fprintf(cof,"\n");
	fprintf(cof,"%8s%8s%8s%8s%8s%8s%8s%8s%14s%14s%14s%14s\n","runs","time","p.N","NT","L","Tb","dE","theta","Num","E","im(action)","W");
	fprintf(cof,"%8i%8g%8i%8i%8g%8g%8g%8g%14.4g%14.4g%14.4g%14.4g\n",\
			runs_count,realtime,p.N,p.NT,p.L,p.Tb,p.DE,p.Theta,Num,E,imag(action),W);
	fprintf(cof,"\n");
	if ((opts.printChoice).compare("gui")==0) {
		printf("\n");
		printf("%8s%8s%8s%8s%8s%8s%8s%8s%14s%14s%14s%14s\n","runs","time","p.N","NT","L","Tb","dE",\
			"theta","Num","E","im(action)","W");
		printf("%8i%8g%8i%8i%8g%8g%8g%8g%14.4g%14.4g%14.4g%14.4g\n",\
			runs_count,realtime,p.N,p.NT,p.L,p.Tb,p.DE,p.Theta,Num,E,imag(action),W);
		printf("\n");
	}

	// printing results to file
	stepped = stepper.keep();
	if (stepped && checkDelta.good()) {
		FILE * actionfile;
		string resultsFile = "results/"+timenumber+"mainResults.dat";
		actionfile = fopen(resultsFile.c_str(),"a");
		fprintf(actionfile,"%12s%5i%5i%5i%5i%6g%13.5g%13.5g%13.5g%13.5g%13.5g%13.5g%13.5g%8.2g%8.2g%8.2g\n",\
					timenumber.c_str(),fileLoop,loop,p.N,p.NT,p.L,p.Tb,p.DE,p.Theta,E,Num,(2.0*imag(action)-bound)\
					,W,checkSoln.back(),checkLin.back(),checkTrue.back());
		fclose(actionfile);
		fprintf(cof,"%12s%30s\n","results:",resultsFile.c_str());
		if ((opts.printChoice).compare("gui")==0)
			printf("%12s%30s\n","results:",resultsFile.c_str());
	}
	stepper.step();
	
	// print everything?, plot too
	bool printEverything = ( ((opts.printChoice).compare("E")==0 || (opts.printChoice).compare("P")==0)? true: false);
	bool plotEverything = ( (opts.printChoice).compare("P")==0? true: false);

	// printing messages for saved files
	so_tp.printMessage = true;
	so_simple.printMessage = true;
	
	// plot options
	PlotOptions po_tp;
	po_tp.gp = "gp/repi.gp";
	po_tp.style = "points";
	Filename plotFile = (string)("data/"+timenumber+"mainp_fLoop_"+numberToString<uint>(fileLoop)\
				+"_loop_"+numberToString<uint>(loop)+".png");
	po_tp.output = plotFile;
	po_tp.printMessage = true;
	
	PlotOptions po_simple;
	po_simple.column = 1;
	po_simple.style = "linespoints";
	po_simple.printMessage = true;
	
	//copying a version of ps with timenumber
	Filename paramsRunFile;
	if (((opts.loopChoice).substr(0,5)).compare("const")!=0) {
		paramsRunFile = "data/"+timenumber+"inputsM_fLoop_"+numberToString<uint>(fileLoop)\
				+"_loop_"+numberToString<uint>(loop);
	}
	else if (stepped) {
		paramsRunFile = "data/"+timenumber+"inputsM_fLoop_"+numberToString<uint>(fileLoop)\
				+"_loop_"+numberToString<uint>(loop)+"_step_1";
	}
	else {
		paramsRunFile = "data/"+timenumber+"inputsM_fLoop_"+numberToString<uint>(fileLoop)\
				+"_loop_"+numberToString<uint>(loop)+"_step_0";
	}
	p.save(paramsRunFile);

	//printing output phi
	Filename tpFile;
	if (((opts.loopChoice).substr(0,5)).compare("const")!=0) {
		tpFile = "data/"+timenumber+"mainp_fLoop_"+numberToString<uint>(fileLoop)\
				+"_loop_"+numberToString<uint>(loop)+".data";
	}
	else if (stepped) {
		tpFile = "data/"+timenumber+"mainp_fLoop_"+numberToString<uint>(fileLoop)\
				+"_loop_"+numberToString<uint>(loop)+"_step_1.data";
	}
	else {
		tpFile = "data/"+timenumber+"mainp_fLoop_"+numberToString<uint>(fileLoop)\
				+"_loop_"+numberToString<uint>(loop)+"_step_0.data";
	}
	save(tpFile,so_tp,f);
	if (plotEverything)
		plot(tpFile,po_tp);
	
	if ((opts.printChoice).compare("n")!=0) {
		//printing linErg
		tpFile.ID = "mainlinErg";
		linErg.conservativeResize(p.Na);
		save(tpFile,so_simple,linErg);
		plotFile = tpFile;
		plotFile.Suffix = ".png";
		po_simple.output = plotFile;
		if (plotEverything)
			plot(tpFile,po_simple);

		//printing erg
		tpFile.ID = "mainerg";
		//erg.conservativeResize(p.Na);
		save(tpFile,so_simple,erg);
		plotFile = tpFile;
		plotFile.Suffix = ".png";
		po_simple.output = plotFile;
		if (plotEverything)
			plot(tpFile,po_simple);
	}

	if (printEverything) {
		//printing output minusDS
		tpFile.ID = "mainminusDS";
		save(tpFile,so_tp,minusDS);
		plotFile = tpFile;
		plotFile.Suffix = ".png";
		po_tp.output = plotFile;
		if (plotEverything)
			plot(tpFile,po_tp);	
			
		//printing output DDS
		tpFile.ID = "mainDDS";
		save(tpFile,so_simple,DDS);
	
		//printing linErgOffShell
		tpFile.ID = "mainlinErgOffShell";
		linErgOffShell.conservativeResize(p.Na);
		save(tpFile,so_simple,linErgOffShell);
		
		//printing linNum
		tpFile.ID = "mainlinNum";
		tpFile.Suffix = ".dat";
		linNum.conservativeResize(p.Na);
		save(tpFile,so_simple,linNum);
	
		//printing derivErg
		tpFile.ID = "mainderivErg";
		derivErg.conservativeResize(p.Na);
		save(tpFile,so_simple,derivErg);
	
		//printing potErg
		tpFile.ID = "mainpotErg";
		potErg.conservativeResize(p.Na);
		save(tpFile,so_simple,potErg);
	}
	
	if (!checkDelta.good()) {
			return 1;
		}
	fprintf(cof,"\n----------------------------------------------------------------------------------------------------------------------------\n\n");
	fclose(cof);
	psu = ps;
	} //ending parameter loop
ces.close();

return 0;
}


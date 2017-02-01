/*----------------------------------------------------------------------------------------------------------------------------
	main_fn3old
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
#include "nr.h"
#include "print3.h"

//#define NDEBUG //NDEBUG is to remove error and bounds checking on vectors in SparseLU, for speed - only include once everything works

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		0 - enums
		1 - argv inputs, loading options, tols
		2 - ParametersRange, Stepper
		3 - beginning parameter loop
		4 - declaring some basic quantities
		5 - assigning potential functions
		6 - omega and negVec
		7 - defining quantitites
		8 - beginning newton-raphson loop
		9 - chiT, chiX
		10 - assigning mds, dds etc
		11 - checks
		12 - printing early 1
		13 - solving for delta
		14 - printing early 2
		15 - convergence
		16 - printing output
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
		- loading tols
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
string fIn;

// cos and cerr files
string coFile, ceFile;

// stepper stuff
string stepperArgv = "";
string stepperInputsFile = "step0";
Filename stepperOutputsFile = "results/nr/step1.csv";
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
cout << "main3" << endl << "====================================================================================================================" << endl;
cout << endl << "options   : " << optionsFile << endl;

// loading tols
tols.load(tolsFile);
cout << "tols      : " << tolsFile << endl;

// other possible inputs
bool extraChecks = false;
bool verbose = true;
bool redo = true;
bool redoErrors = true;
bool step = true;
bool pass = false;
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
		else if (id.compare("tn")==0 || id.compare("timenumber")==0)	 			timenumber = argv[2*j+2];
		else if (id.compare("zmx")==0) 												opts.zmx = argv[2*j+2];
		else if (id.compare("zmt")==0) 												opts.zmt = argv[2*j+2];
		else if (id.compare("bds")==0) 												opts.bds = argv[2*j+2];
		else if (id.compare("inF")==0) 												opts.inF = argv[2*j+2];
		else if (id.compare("epsiTb")==0) 											opts.epsiTb = stn<double>(argv[2*j+2]);
		else if (id.compare("epsiTheta")==0) 										opts.epsiTheta = stn<double>(argv[2*j+2]);
		else if (id.compare("loops")==0) 											opts.loops = stn<uint>(argv[2*j+2]);
		else if (id.compare("printChoice")==0) 										opts.printChoice = argv[2*j+2];
		else if (id.compare("rank")==0) 											rank = argv[2*j+2];
		else if (id.compare("verbose")==0) 											verbose = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("extraChecks")==0 || id.compare("extraTests")==0) 		extraChecks = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("redo")==0) 											redo = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("redoErrors")==0) 										redoErrors = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("step")==0) 											step = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("pass")==0) 											pass = (stn<uint>(argv[2*j+2])!=0);
		else if (id.compare("baseFolder")==0) 										baseFolder = argv[2*j+2];
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
ofstream cos;
cos.open(coFile.c_str(),fstream::app);

ofstream ces;
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
cout << "inputsFile: " << inputsFile << endl;
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

//stepOpts.closeness = tols.Step; // NOT USING THIS

// results
string resultsFile = (pass? "results/nr/nr1pass.csv":"results/nr/nr1.csv");
uint idSizeResults = 2, datumSizeResults = 20;
vector<string> idCheck(idSizeResults);
NewtonRaphsonData results(resultsFile,idSizeResults,datumSizeResults);

// errors
string errorsFile = "results/nr/nr1error.csv";
uint idSizeErrors = 3, datumSizeErrors = 29;
vector<string> idCheckErrors(idSizeErrors);
NewtonRaphsonData errors(errorsFile,idSizeErrors,datumSizeErrors);

// filename suffix
string suffix = ((opts.inF)[(opts.inF).size()-1]=='n'? ".data": ".dat");

// printing timenumber and pot
cos << "timenumber: " << timenumber << endl;
cos << "pot       : " << p.Pot << endl;
cout << "timenumber: " << timenumber << endl;
cout << "pot       : " << p.Pot << endl;

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
cos << "steps     : " << Npl << " steps" << endl << endl;
if (verbose)
	cout << "steps     : " << Npl << " steps" << endl << endl;

// starting loop
for (uint pl=0; pl<Npl; pl++) {
	cos << "--------------------------------------------------------------------------------------------------------------------" << endl;
	if (verbose)
		cout << "--------------------------------------------------------------------------------------------------------------------" << endl;

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
	if (!redo && results.find(p)) {
		if (verbose) {
			cout << "result found in " << resultsFile << " for pl = " << pl << ", ";
			cout << "continuing to next step" << endl;
		}
		continue;
	}
	if (!redoErrors && errors.find(p)) {
		if (verbose) {
			cout << "result found in " << errorsFile << " for pl = " << pl << ", ";
			cout << "continuing to next step" << endl;
		}
		continue;
	}
	
	// getting step file
	if (pl>0) {
		if (opts.inF[0]=='m') {
			stepFile = filenameMain(pold,baseFolder,"field","fMain",suffix);
		}
		else if (opts.inF[0]=='p')
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
		double angleModTwoPi = mod(stepper.stepAngle(),-PI,PI);		
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
	Check checkAction("action",tols.Action);
	Check checkSol("solution",tols.Soln);
	Check checkSolMax("solution max",tols.SolnMax);
	Check checkDelta("delta",tols.Delta);
	Check checkInv("matrix inversion",tols.Inv*p.N*p.NT);
	Check checkCon("energy conservation",tols.Con);
	Check checkLin("linear energy flat",tols.Lin);
	Check checkTrue("linear energy equal true energy",tols.True);
	Check checkLatt("lattice small enough for energy",tols.Latt);
	Check checkReg("regularisation term",tols.Reg);
	Check checkIE("imaginary part of energy",tols.IE);
	Check checkContm("linear energy equal continuum expression",tols.Contm);
	Check checkOS("linear energy equal on shell expression",tols.OS);
	Check checkAB("a_k = Gamma*b_k",tols.AB);
	Check checkABNE("N = Sum(a_k*b_k), E = Sum(w_k*a_k*b_k)",tols.ABNE);
	Check checkLR("linear representation of phi",tols.LR);
	Check checkBoundRe("initial real boundary condition",tols.Soln);
	Check checkBoundIm("initial imaginary boundary condition",tols.Soln);
	Check checkChiT("time zero mode condition",tols.Soln);
	
	// some derived quantities
	uint zm = (p.Pot<3? 2: 1); // number of zero modes
	uint Len = 2*(p.N*p.NT+zm); // not sure about this. do i impose a real or a complex equation to fix the zero mode?

/*----------------------------------------------------------------------------------------------------------------------------
5. assigning potential functions
	- assigning potential functions
	- assigning preliminary parameter structs
----------------------------------------------------------------------------------------------------------------------------*/

	// assigning potential functions
	Potential<comp> V, dV, ddV;
	if (p.Pot==1) {
		V((Potential<comp>::PotentialType)&V1<comp>,p);
		dV((Potential<comp>::PotentialType)&dV1<comp>,p);
		ddV((Potential<comp>::PotentialType)&ddV1<comp>,p);
	}
	else if (p.Pot==2) {
		V((Potential<comp>::PotentialType)&V2<comp>,p);
		dV((Potential<comp>::PotentialType)&dV2<comp>,p);
		ddV((Potential<comp>::PotentialType)&ddV2<comp>,p);
	}
	else if (p.Pot==3) {
		V((Potential<comp>::PotentialType)&V3<comp>,p);
		dV((Potential<comp>::PotentialType)&dV3<comp>,p);
		ddV((Potential<comp>::PotentialType)&ddV3<comp>,p);
	}
	else {
		ces << "pot option not available, pot = " << p.Pot << endl;
		if (verbose)
			cerr << "pot option not available, pot = " << p.Pot << endl;
		return 1;
	}

	//Vr
	Potential<comp> Vr, dVr, ddVr;
	Vr((Potential<comp>::PotentialType)&Vreg<comp>,p);
	dVr((Potential<comp>::PotentialType)&dVreg<comp>,p);
	ddVr((Potential<comp>::PotentialType)&ddVreg<comp>,p);
	
	// assigning preliminary parameter structs
	params_for_V paramsV  = {p.epsilon, p.A};
	params_for_V paramsVr  = {p.minima[0], p.minima[1]};
	V.setParams(paramsV);
	Vr.setParams(paramsVr);

/*----------------------------------------------------------------------------------------------------------------------------
6. omega and negVec
	- loading or constructing omega
	- loading negVec
----------------------------------------------------------------------------------------------------------------------------*/

	//deterimining omega matrices for fourier transforms in spatial direction
	vec freqs(p.N), freqs_exp(p.N);
	mat modes(p.N,p.N);
	mat omega_m1(p.N,p.N), omega_0(p.N,p.N), omega_1(p.N,p.N), omega_2(p.N,p.N);
	
	Filename omegaM1F = filenameSpatial(p,baseFolder,"omega","omegaM1",suffix);
	Filename omega0F = filenameSpatial(p,baseFolder,"omega","omega0F",suffix);
	Filename omega1F = filenameSpatial(p,baseFolder,"omega","omega1F",suffix);
	Filename omega2F = filenameSpatial(p,baseFolder,"omega","omega2F",suffix);
	Filename modesF = filenameSpatial(p,baseFolder,"omega","modesF",suffix);
	Filename freqsF = filenameSpatial(p,baseFolder,"omega","freqsF",suffix);
	Filename freqsExpF = filenameSpatial(p,baseFolder,"omega","freqsExpF",suffix);
	
	if (omegaM1F.exists() && omega0F.exists() && omega1F.exists() && omega2F.exists() \
			&& modesF.exists() && freqsF.exists() && freqsExpF.exists()) {
			loadMatrixBinary(omegaM1F,omega_m1);
			loadMatrixBinary(omega0F,omega_0);
			loadMatrixBinary(omega1F,omega_1);
			loadMatrixBinary(omega2F,omega_2);
			loadMatrixBinary(modesF,modes);
			loadVectorBinary(freqsF,freqs);
			loadVectorBinary(freqsExpF,freqs_exp);
	}
	else {
		bool approxOmega = false;
		if (!approxOmega) {
			numericalModes(modes,freqs,freqs_exp,p);
		}
		else {
			analyticModes(modes,freqs,freqs_exp,p);
		}
		omegasFn(approxOmega,modes,freqs,omega_m1,omega_0,omega_1,omega_2,p);
		saveMatrixBinary(omegaM1F,omega_m1);
		saveMatrixBinary(omega0F,omega_0);
		saveMatrixBinary(omega1F,omega_1);
		saveMatrixBinary(omega2F,omega_2);
		saveMatrixBinary(modesF,modes);
		saveVectorBinary(freqsF,freqs);
		saveVectorBinary(freqsExpF,freqs_exp);
	}
		
	// loading negVec
	vec negVec;
	if (opts.zmt[0]=='n' || opts.zmx[0]=='n') {
		Filename negVecFile = baseFolder+"data/00"+rank+"eigVec_pot_3_L_"+nts(p.L)+".dat";
		//Filename negVecFile = filenameSpatial(p,baseFolder,"eigenvector","negVec",suffix);
		fprintf(cof,"%12s%28s\n","negVecFile: ",((string)negVecFile).c_str());
		if (verbose)
			printf("%12s%28s\n","negVecFile: ",((string)negVecFile).c_str());
		if (negVecFile.exists()) {
			//loadVectorBinary(negVecFile,negVec);
			loadVectorAscii(negVecFile,negVec);
		}
		else {
			ces << "Error: negVecFile doesn't exist." << endl;
			if (verbose)
				cerr << "Error: negVecFile doesn't exist." << endl;
		}
		if (negVec.size()!=p.N)
			negVec = interpolate1d(negVec,p.N);
		// THERE WAS A BUNCH MORE STUFF HERE FOR POT!=3 WHICH I HAVE REMOVED, SEE MAIN_FN.CC
	}

/*----------------------------------------------------------------------------------------------------------------------------
7. defining quantities
	- erg, linErg etc
	- f, mds, dds
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
	comp kineticS = 0.0;
	comp kineticT = 0.0;
	comp potV = 0.0;
	comp potr = 0.0;
	comp bound = 0.0;
	double W = 0.0;
	double E = 0.0;
	double E_exact = 0.0;
	double Num = 0.0;
	
	//defining some quantities used to stop the Newton-Raphson loop when action stops varying
	comp action_last = action;
	uint runsCount = 0;
	uint min_runs = 1;
	uint max_runs = 100;
	
	if (pl==0 || !step || pass) {
		if (!fIn.empty())
			loadFile = fIn;
		else if (opts.inF[0]=='m')
			loadFile = filenameMain(p,baseFolder,"field","fMain",suffix);
		else if (opts.inF[0]=='p')
			loadFile = filenameMain(p,baseFolder,"field","fpi",suffix);
	}
	else 
		loadFile = stepFile;
	
	if (!loadFile.exists()) {
		ces << "nrmain error: " << loadFile << " doesn't exist on pl = " << pl << ", moving to next parameter loop" << endl;
		cerr << "nrmain error: " << loadFile << " doesn't exist on pl = " << pl << ", moving to next parameter loop" << endl;
		continue; ///////// CONTINUE STATEMENT IF FILE DOESN'T EXIST
	}
	
	fprintf(cof,"%12s%30s\n\n","loadFile: ",((string)loadFile).c_str());
	if (verbose)
		printf("%12s%30s\n\n","loadFile: ",((string)loadFile).c_str());
		
	//initializing phi (=f)
	vec f;
		
	// loading f
	loadVectorBinary(loadFile,f);
	if ((pold.NT!=p.NT || pold.N!=p.N) && loadFile==stepFile) {
		f = interpolate(f, pold, p);
		fprintf(cof,"interpolating\n");
		if (verbose)
			printf("interpolating\n");
	}
	if (f.size()<Len) {
		uint oldLen = f.size();
		f.conservativeResize(Len);
		for (uint mu=0; mu<2*zm; mu++)
			f[2*p.N*p.NT+mu] = 1.0e-4;
		fprintf(cof,"resizing from %u to %u\n",oldLen,Len);
		if (verbose)
			printf("resizing from %u to %u\n",oldLen,Len);
	}
	else if (f.size()>Len)
		f.conservativeResize(Len);
	// MAY BE WORTH INTRODUCING A STRETCHING HERE IF Tb OR LoR IS BEING CHANGED
	
	// printing parameters
	p.print(cof);
	if (verbose)
		p.print();
	
	//defining complexified vector Cf
	cVec Cf;
	Cf = vecComplex(f,p);
	

	//defining dds and mds
	spMat dds(Len,Len);
	vec mds(Len);
	
	// defining a couple of vectors to test whether the initial boundary conditions are satisfied
	vec boundRe(p.N);
	vec boundIm(p.N);
		
	//very early vector print
	if ((opts.printChoice).compare("n")!=0) {
		Filename earlyPrintFile = filenameMain(p,baseFolder,"ascii","fMainEarly",".dat");
		(earlyPrintFile.Extras).push_back(StringPair("run","0"));
		saveVectorAscii(earlyPrintFile,f);
	}
	
/*----------------------------------------------------------------------------------------------------------------------------
8. beginning newton-raphson loop
	- beginning newton-raphson loop	
	- reserving space in dds
	- zeroing erg etc
	- crude test of potential (to remove once done)
----------------------------------------------------------------------------------------------------------------------------*/

	bool passThrough = false;
	//beginning newton-raphson loop	
	while ((!checkSol.good() || !checkSolMax.good() || runsCount<min_runs) && !passThrough) {
		passThrough = pass;
		runsCount++;
		if (runsCount>max_runs) {
			ces << "main error: max_runs(" << max_runs << ") exceeded for:" << endl;
			ces << "timenumber: " << timenumber << "; loop: " << pl << endl;
			if (verbose) {
				cerr << "main error: max_runs(" << max_runs << ") exceeded for:" << endl;
				cerr << "timenumber: " << timenumber << "; loop: " << pl << endl;
			}
			continue;				
		}
		
		// allocating memory for DS, dds
		mds = Eigen::VectorXd::Zero(Len); //initializing to zero
		dds.setZero(); //just making sure
		Eigen::VectorXi dds_to_reserve(Len);//number of non-zero elements per column
		dds_to_reserve = Eigen::VectorXi::Constant(Len,13);
		if (abs(p.Theta)<2.0e-16) {
			dds_to_reserve(0) = p.N+3;
			dds_to_reserve(1) = 11;
		}
		else {
			dds_to_reserve(0) = p.N+11;
			dds_to_reserve(1) = p.N+11;
		}
		dds_to_reserve(2*p.N*p.NT-2) = 4;
		dds_to_reserve(2*p.N*p.NT-1) = 4;
		dds_to_reserve(2*p.N*p.NT) = p.N*((opts.zmx).size()-2);
		dds_to_reserve(2*p.N*p.NT+1) = 2*p.N*((opts.zmt).size()-2);
		dds.reserve(dds_to_reserve);
		
		//initializing to zero
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
		kineticS = 0.0;
		kineticT = 0.0;
		potV = 0.0;
		potr = 0.0;
		boundRe = Eigen::VectorXd::Zero(p.N);
		boundIm = Eigen::VectorXd::Zero(p.N);
		
/*----------------------------------------------------------------------------------------------------------------------------
9. chiT, chiX
----------------------------------------------------------------------------------------------------------------------------*/
		
		// zero modes - fixed with chiX and chiT
		vec chiX(p.NT*p.N);	chiX = Eigen::VectorXd::Zero(p.N*p.NT); //to fix spatial zero mode
		vec chiT(p.NT*p.N);	chiT = Eigen::VectorXd::Zero(p.N*p.NT); //to fix real time zero mode
		for (uint j=0; j<p.N; j++) {
			uint posX, posT, posCe;
			uint slicesX, slicesT;
			if (getLastInt(opts.zmx)<0) {
				ces << "getLastInt error with zmx = " << opts.zmx << endl;
				if (verbose)
					cerr << "getLastInt error with zmx = " << opts.zmx << endl;
				return 1;
			}
			if (getLastInt(opts.zmt)<0) {
				ces << "getLastInt error with zmt = " << opts.zmt << endl;
				if (verbose)
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
				if (verbose)
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
							if (verbose)
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
				if (verbose)
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
							chiX(posX+k) = f(2*neigh(posX+k,1,1,p))-f(2*neigh(posX+k,1,-1,p));
						else {
							ces << "choice of zmx(" << opts.zmx << ") not allowed" << endl;
							if (verbose)
								cerr << "choice of zmx(" << opts.zmx << ") not allowed" << endl;
							return 1;
						}
					}
				}
			}
		}
		double normX = chiX.norm();
		double normT = chiT.norm();
		if (abs(normX)<MIN_NUMBER || abs(normT)<MIN_NUMBER) {
			ces << "norm of chiX = " << normX << ", norm of chiT = " << normT << endl;
			if (verbose)
				cerr << "norm of chiX = " << normX << ", norm of chiT = " << normT << endl;
		}
		chiX = chiX/normX;
		chiT = chiT/normT;
		
/*----------------------------------------------------------------------------------------------------------------------------
9.1 some extra checks, mostly trivial
----------------------------------------------------------------------------------------------------------------------------*/
		
		if (extraChecks) {
			//testing that the potential term is working for pot3
			if (p.Pot==3) {
				comp Vtrial = 0.0, Vcontrol = 0.0;
				for (unsigned int j=0; j<p.N; j++) {
					double r = p.r0 + j*p.a;
					paramsV.epsi = r;
					V.setParams(paramsV);
					Vcontrol += pow(f(2*j*p.Nb),2.0)/2.0 - pow(f(2*j*p.Nb),4.0)/4.0/pow(r,2.0);
					Vtrial += V(f(2*j*p.Nb));
				}
				double potTest = pow(pow(real(Vcontrol-Vtrial),2.0) + pow(imag(Vcontrol-Vtrial),2.0),0.5);
				fprintf(cof,"potTest         = %40.20g\n",potTest);
				if (verbose)
					printf("potTest         = %40.20g\n",potTest);
			}
			fprintf(cof,"f.norm()        = %40.20g\n",f.norm());
			fprintf(cof,"chiT.norm()    = %40.20g\n",normX);
			fprintf(cof,"chiX.norm()     = %40.20g\n",normT);
			if (verbose) {
				printf("f.norm()        = %40.20g\n",f.norm());
				printf("chiT.norm()     = %40.20g\n",normX);
				printf("chiX.norm()     = %40.20g\n",normT);
			}
		}

/*----------------------------------------------------------------------------------------------------------------------------
10.0 coefficients for assigning mds, dds etc
----------------------------------------------------------------------------------------------------------------------------*/

		// coefficient vectors
		cVec coeff_kineticT(p.N*p.NT);
		cVec coeff_kineticS(p.N*p.NT);
		cVec coeff_pot(p.N*p.NT);
		cVec coeff_potr(p.N*p.NT);
		cVec coeff_kineticT_tnt(p.N*p.NT);
		cVec coeff_kineticS_tnt(p.N*p.NT);
		cVec coeff_pot_tnt(p.N*p.NT);
		cVec coeff_potr_tnt(p.N*p.NT);
		cVec coeff_kineticT_erg(p.N*p.NT);
		cVec coeff_kineticS_erg(p.N*p.NT);
		cVec coeff_pot_erg(p.N*p.NT);
		cVec coeff_potr_erg(p.N*p.NT);
		double tnt = 1.0;
		//if ((opts.bds).compare("uc")!=0) // removed this requirement because multiplication by Theta is built into the boundary functions in nr.cc
		tnt *= p.Theta;
				
		// preparatory loop over lattice points
		for (lint j = 0; j < p.N*p.NT; j++) {
			uint t 					= intCoord(j,0,p); //coordinates
			uint x 					= intCoord(j,1,p);
			
			// coefficients in action
			coeff_kineticT[j] 		= DxFn(x,p)/dtFn(t,p);
			coeff_kineticS[j] 		= -DtFn(t,p)/dxFn(x,p);
			coeff_pot[j]	 		= -DtFn(t,p)*DxFn(x,p);
			coeff_potr[j]	 		= (ii*p.Reg)*DtFn(t,p)*DxFn(x,p);
			
			// with Theta multiplying the real equation
			coeff_kineticT_tnt[j] 	= comp(tnt*real(coeff_kineticT[j]),imag(coeff_kineticT[j]));
			coeff_kineticS_tnt[j] 	= comp(tnt*real(coeff_kineticS[j]),imag(coeff_kineticS[j]));
			coeff_pot_tnt[j]	 	= comp(tnt*real(coeff_pot[j]),imag(coeff_pot[j]));
			coeff_potr_tnt[j]	 	= comp(tnt*real(coeff_potr[j]),imag(coeff_potr[j]));
			
			// for the energy
			// the following 4 expressions are discretisations of continuum expressions (can we do better than this?)
			// We could just derive the energy as a Noether charge.
			coeff_kineticT_erg[j] 	= DxFn(x,p)/pow(dtFn(t,p),2.0);
			coeff_kineticS_erg[j] 	= 1.0/dxFn(x,p);
			coeff_pot_erg[j]	 	= DxFn(x,p);
			coeff_potr_erg[j]	 	= -ii*p.Reg*DxFn(x,p);
		}

/*----------------------------------------------------------------------------------------------------------------------------
10. assigning mds, dds etc
	- beginning loop over lattice points
	- fixing zero modes
	- linErg, linNum
	- t=(NT-1)
	- t=0
	- bulk
	- extras (*4.0*PI, etc)
	- E, N, bound, W
----------------------------------------------------------------------------------------------------------------------------*/
		
		// beginning loop over lattice points
		for (lint j = 0; j < p.N*p.NT; j++) {
			//coordinates
			uint t 					= intCoord(j,0,p); 
			uint x 					= intCoord(j,1,p);
			comp dt 		= dtFn(t,p);
			double Dx 		= DxFn(x,p);
			
			// assigning parameters
			if (p.Pot==3) {
				paramsV.epsi = p.r0+x*p.a;
				V.setParams(paramsV);
				dV.setParams(paramsV);
				ddV.setParams(paramsV);
			}
			
			// lagrange multiplier terms
			if (abs(chiX(j))>MIN_NUMBER && p.Pot!=3) { //spatial zero mode lagrange constraint
				dds.insert(2*j,2*p.N*p.NT) 		= Dx*chiX(j); 
				dds.insert(2*p.N*p.NT,2*j) 		= Dx*chiX(j);
				mds(2*j) 						+= -Dx*chiX(j)*f(2*p.N*p.NT);
				mds(2*p.N*p.NT) 				+= -Dx*chiX(j)*f(2*j);
			}
			if (abs(chiT(j))>MIN_NUMBER && t<(p.NT-1) && (p.Pot!=3 || (x>0 && x<(p.N-1)))) {
				dds.coeffRef(2*(j+1),2*p.N*p.NT+1) 	+= Dx*chiT(j); //chiT should be 0 at t=(p.NT-1) or this line will go wrong
				dds.coeffRef(2*p.N*p.NT+1,2*(j+1)) 	+= Dx*chiT(j);
				dds.coeffRef(2*j,2*p.N*p.NT+1) 		+= -Dx*chiT(j);
				dds.coeffRef(2*p.N*p.NT+1,2*j) 		+= -Dx*chiT(j);
	            mds(2*(j+1)) 						+= -Dx*chiT(j)*f(2*p.N*p.NT+1);
	            mds(2*j) 							+= Dx*chiT(j)*f(2*p.N*p.NT+1);
	            mds(2*p.N*p.NT+1) 					+= -Dx*chiT(j)*(f(2*(j+1))-f(2*j));
			}
				
			// linear energy and particle number
			if (t<(p.NT-1)) {
				for (uint k=0;k<p.N;k++) {
					lint l = k*p.NT+t;
					linErgOffShell(t) += 0.5*( omega_2(x,k)*( Cf(l)-p.minima[0] )*( Cf(j)-p.minima[0] ) \
									+ omega_0(x,k)*( Cf(l+1)-Cf(l) )*( Cf(j+1)-Cf(j) )/pow(dt,2.0)); // Are we sure about the 1/dt^2? Are there any factors of Dt?
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

			//boundaries
			if (p.Pot==3 && x==(p.N-1)) {
				dds.insert(2*j,2*j) 	= 1.0; // f=0 at r=R
				dds.insert(2*j+1,2*j+1) = 1.0;
			}
			else if (p.Pot==3 && x==0) {
				Kinetic_nr (j, 1, f, p, coeff_kineticS_erg, derivErg, t); //#### The hashes signify that this term is new and untested.
				
				Kinetic_nr (j, 1, f, p, coeff_kineticS, kineticS);
				dds.insert(2*j,2*j) 	= 1.0; // f=0 at r=0
				dds.insert(2*j+1,2*j+1) = 1.0;
			}			
			else if (t==(p.NT-1)) {		
				Kinetic_nr (j, 1, f, p, coeff_kineticS_erg, derivErg, t);	//####
				Potential_nr (j, f, p, V, coeff_pot_erg, potErg, t);	//####
				Potential_nr (j, f, p, Vr, coeff_potr_erg, potErg, t);	//####
				
				Kinetic_nr (j, 1, f, p, coeff_kineticS, kineticS);
				Potential_nr (j, f, p, V, coeff_pot, potV);
				Potential_nr (j, f, p, Vr, coeff_potr, potr);
			
				dds.insert(2*j,2*(j-1)+1) = 1.0; //zero imaginary part of time derivative
				dds.insert(2*j+1,2*j+1)   = 1.0; //zero imaginary part
			}
			else if (t==0) {
				Kinetic_nr (j, 0, f, p, coeff_kineticT_erg, derivErg, t);	//####
				Kinetic_nr (j, 1, f, p, coeff_kineticS_erg, derivErg, t);	//####
				Potential_nr (j, f, p, V, coeff_pot_erg, potErg, t);		//####
				Potential_nr (j, f, p, Vr, coeff_potr_erg, potErg, t);	//####
				
				Kinetic_nr (j, 0, f, p, coeff_kineticT, kineticT);
				Kinetic_nr (j, 1, f, p, coeff_kineticS, kineticS);
				Potential_nr (j, f, p, V, coeff_pot, potV);
				Potential_nr (j, f, p, Vr, coeff_potr, potr);
				Boundary_nr (j, f, p, omega_1, bound);
				
				mdKinetic_nr (j, 0, f, p, coeff_kineticT, mds, Complex_nr::imaginary);
				mdKinetic_nr (j, 1, f, p, coeff_kineticS, mds, Complex_nr::imaginary);
				mdPotential_nr (j, f, p, dV, coeff_pot, mds, Complex_nr::imaginary);
				mdPotential_nr (j, f, p, dVr, coeff_potr, mds, Complex_nr::imaginary);
				mdBoundary_nr (j, f, p, omega_1, mds, Complex_nr::both);
				
				ddKinetic_nr (j, 0, f, p, coeff_kineticT, dds, Complex_nr::imaginary);
				ddKinetic_nr (j, 1, f, p, coeff_kineticS, dds, Complex_nr::imaginary);
				ddPotential_nr (j, f, p, ddV, coeff_pot, dds, Complex_nr::imaginary);
				ddPotential_nr (j, f, p, ddVr, coeff_potr, dds, Complex_nr::imaginary);
				ddBoundary_nr (j, f, p, omega_1, dds, Complex_nr::real);
				
				if (abs(p.Theta)>MIN_NUMBER) {
					mdKinetic_nr (j, 0, f, p, coeff_kineticT_tnt, mds, Complex_nr::real);
					mdKinetic_nr (j, 1, f, p, coeff_kineticS_tnt, mds, Complex_nr::real);
					mdPotential_nr (j, f, p, dV, coeff_pot_tnt, mds, Complex_nr::real);
					mdPotential_nr (j, f, p, dVr, coeff_potr_tnt, mds, Complex_nr::real);
					
					ddKinetic_nr (j, 0, f, p, coeff_kineticT_tnt, dds, Complex_nr::real);
					ddKinetic_nr (j, 1, f, p, coeff_kineticS_tnt, dds, Complex_nr::real);
					ddPotential_nr (j, f, p, ddV, coeff_pot_tnt, dds, Complex_nr::real);
					ddPotential_nr (j, f, p, ddVr, coeff_potr_tnt, dds, Complex_nr::real);
				}
			}
			//bulk
			else {
				Kinetic_nr (j, 0, f, p, coeff_kineticT_erg, derivErg, t);	//####
				Kinetic_nr (j, 1, f, p, coeff_kineticS_erg, derivErg, t);	//####
				Potential_nr (j, f, p, V, coeff_pot_erg, potErg, t);	//####
				Potential_nr (j, f, p, Vr, coeff_potr_erg, potErg, t);	//####
            	
				Kinetic_nr (j, 0, f, p, coeff_kineticT, kineticT);
				Kinetic_nr (j, 1, f, p, coeff_kineticS, kineticS);
				Potential_nr (j, f, p, V, coeff_pot, potV);
				Potential_nr (j, f, p, Vr, coeff_potr, potr);
				
				mdKinetic_nr (j, 0, f, p, coeff_kineticT, mds);
				mdKinetic_nr (j, 1, f, p, coeff_kineticS, mds);
				mdPotential_nr (j, f, p, dV, coeff_pot, mds);
				mdPotential_nr (j, f, p, dVr, coeff_potr, mds);
				
				ddKinetic_nr (j, 0, f, p, coeff_kineticT, dds);
				ddKinetic_nr (j, 1, f, p, coeff_kineticS, dds);
				ddPotential_nr (j, f, p, ddV, coeff_pot, dds);
				ddPotential_nr (j, f, p, ddVr, coeff_potr, dds);
	        }
	    } // end of loop over j
	    
	    if (p.Pot==3)
	    	dds.insert(2*p.N*p.NT,2*p.N*p.NT) = 1.0;
	    
	    action = kineticT + kineticS + potV + potr;
	    linErgOffShell(p.NT-1) = linErgOffShell(p.NT-2);
    	linNumOffShell(p.NT-1) = linNumOffShell(p.NT-2);
    	
    	// assigning energy
		erg = derivErg + potErg;
    	
	    if (p.Pot==3) {
	    	action 			*= 4.0*PI;
	    	derivErg 		*= 4.0*PI;
	    	potErg 			*= 4.0*PI;
	    	erg				*= 4.0*PI;
	    	linErg			*= 4.0*PI;
	    	linNum			*= 4.0*PI;
	    	linNumOffShell 	*= 4.0*PI;
	    	linErgOffShell 	*= 4.0*PI;
	    }
	   
	    if (abs(p.Theta)<MIN_NUMBER) {
	    	linErg = linErgOffShell;
	    	linNum = linNumOffShell;
	    }
	    
	    for (uint x=1; x<(p.N-1); x++) {
	    	lint m = x*p.NT;
	    	boundRe(x) = mds(2*m);
	    	boundIm(x) = mds(2*m+1);
	    }
	    
	    //defining E, Num and W
		E = real(linErg(0));
		Num = real(linNum(0));
		W = - E*2.0*p.Tb - p.Theta*Num - imag(bound) + 2.0*imag(action);
		
/*----------------------------------------------------------------------------------------------------------------------------
11. checks
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
		if (extraChecks) {
			for (uint j=0; j<p.NT; j++) {
				double diff = absDiff(erg(j),potErg(j)+derivErg(j));
				if (diff>1.0e-14) {
					ces << "erg(" << j << ") != potErg + derivErg. absDiff = " << diff << endl;
					if (verbose)
						cerr << "erg(" << j << ") != potErg + derivErg. absDiff = " << diff << endl;
				}
			}
		}
		
		//calculating continuum approx to linErg and linNum on initial time slice - redundant
		if (p.Pot==3 && abs(p.Theta)<MIN_NUMBER) {
			for (uint k=1; k<p.N; k++) {
				double momtm = k*PI/(p.L-p.r0);
				double freqSqrd = 1.0+pow(momtm,2.0);
				double Asqrd, integral1 = 0.0, integral2 = 0.0;
				for (unsigned int l=0; l<p.N; l++) {
					double r = p.r0 + l*p.a;
					lint m = l*p.NT;
					integral1 += p.a*f(2*m)*pow(2.0/p.L,0.5)*sin(momtm*r);
					integral2 += p.a*(f(2*(m+1))-f(2*m))*pow(2.0/p.L,0.5)*sin(momtm*r)/p.b;
				}
				Asqrd = pow(integral1,2.0) + pow(integral2,2.0)/freqSqrd;
				linErgContm += 2.0*PI*Asqrd*freqSqrd;
				linNumContm += 2.0*PI*Asqrd*pow(freqSqrd,0.5);
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
		comp dt0 = dtFn(0,p);
		for (uint n=0; n<p.N; n++) {
			double w_n = freqs(n);
			double w_n_e = freqs_exp(n);
			for (uint j=0; j<p.N; j++) {
				lint m=j*p.NT;
				double sqrtDj = sqrt(4.0*PI*DxFn(j,p));
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
			if (extraChecks) {
				for (uint n=0; n<p.N; n++) {
					double w_n = freqs(n);
					double w_n_e = freqs_exp(n);
					double sqrtDj = sqrt(4.0*PI*DxFn(j,p));
					if (abs(w_n)>1.0e-16) {
						linRep(j) += modes(j,n)*(a_k(n)*exp(-ii*w_n_e*T0)+b_k(n)*exp(ii*w_n_e*T0)) \
										/sqrt(2.0*w_n)/sqrtDj;
					}
				}
			}
		}
		if (extraChecks) {
			linRepTest = absDiff(p0,linRep);
			checkLR.add(linRepTest);
		}
		ABtest /= (double)p.N;
		checkAB.add(ABtest);
		double ABNtest = absDiff(linNumAB,linNum(0));
		double ABEtest = absDiff(linErgAB,linErg(0));
		ABNtest>ABEtest? checkABNE.add(ABNtest): checkABNE.add(ABEtest);
		
		//checking potr is much smaller than the other potential terms
		checkReg.add(abs(potr/potV));
		checkReg.checkMessage(ces);
		if (verbose)
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
			if (verbose)
				cerr << "E = " << E << ", E_exact = " << E_exact << ", linearInt = " << linearInt << endl;
		}
		
		//checking lattice small enough for E, should have parameter for this
		double momTest = E*p.b/Num/PI; //perhaps should have a not b here
		checkLatt.add(momTest);
		
		//checking initial boundary conditions satisfied
		double normP = f.norm();
		double boundReTest = boundRe.norm()*(Len)/normP/(p.N-2.0);
		double boundImTest = boundIm.norm()*(Len)/normP/(p.N-2.0);
		checkBoundRe.add(boundReTest);
		checkBoundIm.add(boundImTest);
		
		// checking chiT orthogonality satisfied
		double chiTTest = mds(2*p.N*p.NT+1)*(Len)/normP;
		checkChiT.add(chiTTest);
		
/*----------------------------------------------------------------------------------------------------------------------------
11.1 some extra checks, mostly trivial
----------------------------------------------------------------------------------------------------------------------------*/		
		
		// some extra checks, mostly trivial
		if (extraChecks) {
	    	fprintf(cof,"imag(action)    = %40.20g\n",imag(action));
			fprintf(cof,"imag(kineticT)  = %40.20g\n",imag(kineticT));
			fprintf(cof,"imag(kineticS)  = %40.20g\n",imag(kineticS));
			fprintf(cof,"imag(potV)      = %40.20g\n",imag(potV));
			fprintf(cof,"imag(potr)      = %40.20g\n",imag(potr));
			fprintf(cof,"E               = %40.20g\n",E);
			fprintf(cof,"Num             = %40.20g\n",Num);
			fprintf(cof,"derivErg.norm() = %40.20g\n",derivErg.norm());
			fprintf(cof,"potErg.norm()   = %40.20g\n",potErg.norm());
			fprintf(cof,"f.norm()        = %40.20g\n",f.norm());
			fprintf(cof,"mds.norm()      = %40.20g\n",mds.norm());
			fprintf(cof,"dds.norm()      = %40.20g\n",dds.norm());
			if (verbose) {
				printf("imag(action)    = %40.20g\n",imag(action));
				printf("imag(kineticT)  = %40.20g\n",imag(kineticT));
				printf("imag(kineticS)  = %40.20g\n",imag(kineticS));
				printf("imag(potV)      = %40.20g\n",imag(potV));
				printf("imag(potr)      = %40.20g\n",imag(potr));
				printf("E               = %40.20g\n",E);
				printf("Num             = %40.20g\n",Num);
				printf("derivErg.norm() = %40.20g\n",derivErg.norm());
				printf("potErg.norm()   = %40.20g\n",potErg.norm());	
				printf("f.norm()        = %40.20g\n",f.norm());
				printf("mds.norm()      = %40.20g\n",mds.norm());
				printf("dds.norm()      = %40.20g\n",dds.norm());
			}
		}

/*----------------------------------------------------------------------------------------------------------------------------
12. printing early 1
	- filenames and saveoptions
	- phi
	- mds
	- dds
	- chiX, chiT
	- energies
	- a_k, b_k
	
----------------------------------------------------------------------------------------------------------------------------*/

	//printing early if desired
	if ((opts.printChoice).compare("n")!=0) {
		if ((opts.printChoice).compare("v")==0 || (opts.printChoice).compare("e")==0) {
			Filename earlyPrintFile = filenameMain(p,baseFolder,"ascii","mdsMainEarly",".dat");
			(earlyPrintFile.Extras).push_back(StringPair("run",nts(runsCount)));
			saveVectorAscii(earlyPrintFile,mds);
		}
		if ((opts.printChoice).compare("m")==0 || (opts.printChoice).compare("e")==0) {
			Filename earlyPrintFile = filenameMain(p,baseFolder,"ascii","ddsMainEarly",".dat");
			(earlyPrintFile.Extras).push_back(StringPair("run",nts(runsCount)));
			saveSparseMatrixAscii(earlyPrintFile,dds);
		}
		if ((opts.printChoice).compare("z")==0 || (opts.printChoice).compare("e")==0) {
			Filename earlyPrintFile = filenameMain(p,baseFolder,"ascii","chiTMainEarly",".dat");
			(earlyPrintFile.Extras).push_back(StringPair("run",nts(runsCount)));
			saveVectorAscii(earlyPrintFile,chiT);
			earlyPrintFile = filenameMain(p,baseFolder,"ascii","chiXMainEarly",".dat");
			(earlyPrintFile.Extras).push_back(StringPair("run",nts(runsCount)));
			saveVectorAscii(earlyPrintFile,chiX);
		}
		if ((opts.printChoice).compare("l")==0 || (opts.printChoice).compare("e")==0) {
			Filename earlyPrintFile = filenameMain(p,baseFolder,"ascii","linErgMainEarly",".dat");
			(earlyPrintFile.Extras).push_back(StringPair("run",nts(runsCount)));
			saveComplexVectorAscii(earlyPrintFile,linErg);
			earlyPrintFile = filenameMain(p,baseFolder,"ascii","ergMainEarly",".dat");
			(earlyPrintFile.Extras).push_back(StringPair("run",nts(runsCount)));
			saveComplexVectorAscii(earlyPrintFile,erg);
		}
		if ((opts.printChoice).compare("ab")==0 || (opts.printChoice).compare("e")==0) {
			Filename earlyPrintFile = filenameMain(p,baseFolder,"ascii","bkMainEarly",".dat");
			(earlyPrintFile.Extras).push_back(StringPair("run",nts(runsCount)));
			saveComplexVectorAscii(earlyPrintFile,a_k);
			earlyPrintFile = filenameMain(p,baseFolder,"ascii","akMainEarly",".dat");
			(earlyPrintFile.Extras).push_back(StringPair("run",nts(runsCount)));
			saveComplexVectorAscii(earlyPrintFile,b_k);
		}
	}
	
	
/*----------------------------------------------------------------------------------------------------------------------------
13. solving for delta	
	- defining delta, solver etc
	- analysing pattern
	- factorising
	- solving
	- checking matrix inversion
	- f' = f + delta
	- Cf'
----------------------------------------------------------------------------------------------------------------------------*/
		
		//solving for delta in dds*delta=mds, where f' = f + delta		
		vec delta(Len);
		delta = Eigen::VectorXd::Zero(Len);
		dds.prune(MIN_NUMBER);
		dds.makeCompressed();
		
		if (!pass) {
			Eigen::SparseLU<spMat> solver;
	
			solver.compute(dds);
			if(solver.info()!=Eigen::Success) {
				ces << "Compute failed, solver.info() = "<< solver.info() << endl;
				if (verbose) {
					cerr << "Compute failed, solver.info() = "<< solver.info() << endl;
					printErrorInformation(f,"f",2);
					cout << endl;
					printErrorInformation(mds,"mds",2);
					cout << endl;
					printErrorInformation(dds,"dds");
					cout << endl;
				}
				return 1;
			}
			delta = solver.solve(mds);// use the factorization to solve for the given right hand side
			if(solver.info()!=Eigen::Success) {
				ces << "Solving failed, solver.info() = "<< solver.info() << endl;
				ces << "log(abs(det(dds))) = " << solver.logAbsDeterminant() << endl;
				ces << "sign(det(dds)) = " << solver.signDeterminant() << endl;
				if (verbose) {
					cerr << "Solving failed, solver.info() = "<< solver.info() << endl;
					cerr << "log(abs(det(dds))) = " << solver.logAbsDeterminant() << endl;
					cerr << "sign(det(dds)) = " << solver.signDeterminant() << endl;
					printErrorInformation(f,"f",2);
					cout << endl;
					printErrorInformation(mds,"mds",2);
					cout << endl;
					printErrorInformation(delta,"delta",2);
					cout << endl;
					printErrorInformation(dds,"dds");
					cout << endl;
				}
				return 1;
			}
	
			//independent check on whether calculation worked
			vec diff(Len);
			diff = dds*delta-mds;
			double maxDiff = diff.maxCoeff();
			maxDiff = abs(maxDiff);
			checkInv.add(maxDiff);
			checkInv.checkMessage(ces);
			if (verbose)
				checkInv.checkMessage(cerr);
			if (!checkInv.good()) return 1;

			//assigning values to phi
			f += delta;
	
			//passing changes on to complex vector
			Cf = vecComplex(f,p.N*p.NT);
		}
		
/*----------------------------------------------------------------------------------------------------------------------------
13.1 some extra checks, mostly trivial
----------------------------------------------------------------------------------------------------------------------------*/		
		
		// some extra checks, mostly trivial
		if (extraChecks) {
			fprintf(cof,"delta.norm()    = %40.20g\n",delta.norm());
			if (verbose) {
				printf("delta.norm()    = %40.20g\n",delta.norm());
			}
		}
		
/*----------------------------------------------------------------------------------------------------------------------------
14. printing early 2
	- filenames and saveoptions
	- delta
	
----------------------------------------------------------------------------------------------------------------------------*/

	//printing early if desired
	if ((opts.printChoice).compare("n")!=0) {
		if ((opts.printChoice).compare("f")==0 || (opts.printChoice).compare("e")==0) {
			Filename earlyPrintFile = filenameMain(p,baseFolder,"ascii","fMainEarly",".dat");
			(earlyPrintFile.Extras).push_back(StringPair("run",nts(runsCount)));
			saveVectorAscii(earlyPrintFile,f);
		}
		if ((opts.printChoice).compare("d")==0 || (opts.printChoice).compare("e")==0) {
			Filename earlyPrintFile = filenameMain(p,baseFolder,"ascii","deltaMainEarly",".dat");
			(earlyPrintFile.Extras).push_back(StringPair("run",nts(runsCount)));
			saveVectorAscii(earlyPrintFile,delta);
		}
	}
	
/*----------------------------------------------------------------------------------------------------------------------------
15. convergence
	- evaluating norms
	- adding convergence checks
	- printing convergence checks
	
----------------------------------------------------------------------------------------------------------------------------*/
		//convergence issues
		//evaluating norms
		double normDS = mds.norm();
		double maxDS = mds.maxCoeff();
		double minDS = mds.minCoeff();
		if (-minDS>maxDS) maxDS = -minDS;
		double maxP = f.maxCoeff();
		double minP = f.minCoeff();
		if (-minP>maxP) maxP = -minP;
	
		// adding convergence checks
		checkAction.add(absDiff(action,action_last));
		action_last = action;
		checkSol.add(normDS/normP);
		checkSolMax.add(maxDS/maxP);
		
		double normDelta = delta.norm();
		checkDelta.add(normDelta/normP);

		// rate of convergence
		double convergRate = (checkDelta.size()>1? \
			log((checkDelta.tests())[checkDelta.size()-1])/log((checkDelta.tests())[checkDelta.size()-2]):0.0);
		
		//printing tests to see convergence
		if (runsCount==1) {
			fprintf(cof,"%5s%5s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s\n","loop","run","sol","solM","delta","converg"\
			,"linear","true erg","on shell","AB","ABNE","conserv","latt");
			if (verbose)
				printf("%5s%5s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s\n","loop","run","sol","solM","delta","converg"\
				,"linear","true erg","on shell","AB","ABNE","conserv","latt");
		}
		fprintf(cof,"%5i%5i%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g\n",pl,runsCount\
		,checkSol.back(),checkSolMax.back(),checkDelta.back(),convergRate,checkLin.back(),checkTrue.back()\
		,checkOS.back(),checkAB.back(),checkABNE.back(),\
			checkCon.back(),checkLatt.back());
		if (verbose)
			printf("%5i%5i%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g\n",pl,runsCount\
			,checkSol.back(),checkSolMax.back(),checkDelta.back(),convergRate,checkLin.back(),checkTrue.back()\
			,checkOS.back(),checkAB.back(),checkABNE.back(),checkCon.back(),checkLatt.back());
		
		if (!checkDelta.good() && !pass) {
			checkDelta.checkMessage(ces);
			if (verbose)
				checkDelta.checkMessage(cerr);
			break;
		}
		
	} //ending while loop
/*----------------------------------------------------------------------------------------------------------------------------
16. printing output
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
			- mds
			- dds
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
	if (p.Pot==3 && abs(p.Theta)<MIN_NUMBER)
		checkContm.checkMessage(ces);
	checkCon.checkMessage(ces);
	checkIE.checkMessage(ces);
	checkLatt.checkMessage(ces);
	checkBoundRe.checkMessage(ces);
	checkBoundIm.checkMessage(ces);
	if (verbose) {
		checkTrue.checkMessage(cerr);
		checkOS.checkMessage(cerr);
		checkABNE.checkMessage(cerr);
		if (p.Pot==3 && abs(p.Theta)<MIN_NUMBER)
			checkContm.checkMessage(cerr);
		checkCon.checkMessage(cerr);
		checkIE.checkMessage(cerr);
		checkLatt.checkMessage(cerr);
		checkBoundRe.checkMessage(cerr);
		checkBoundIm.checkMessage(cerr);
	}
	
	//stopping clock
	time = clock() - time;
	double realtime = time/1000000.0;
	
	// printing stepper results
	if (stepargv!=StepperArgv::none) {
		// adding result to stepper
		number F = 0.0;
		if(stepargv==StepperArgv::w)
			F = W;
		else if(stepargv==StepperArgv::e)
			F = E;
		else if(stepargv==StepperArgv::n)
			F = Num;
		else {
			cerr << "main_fn3 error: stepargv, " << stepargv << ", not recognized" << endl;
			return 1;
		}
		stepper.addResult(F);
		
		// printing step
		stepper.save(stepperOutputsFile);
		fprintf(cof,"%12s%24s\n","stepper outputs:",((string)stepperOutputsFile).c_str());
		if (verbose)
			printf("%12s%24s\n","stepper outputs:",((string)stepperOutputsFile).c_str());
		
		// stepping
		stepper.step();
	
	}

	// printing results to terminal
	fprintf(cof,"\n");
	fprintf(cof,"%8s%8s%8s%8s%8s%8s%8s%8s%14s%14s%14s%14s\n","runs","time","p.N","NT","L","Tb","dE","theta","Num","E","im(action)","W");
	fprintf(cof,"%8i%8g%8i%8i%8g%8g%8g%8g%14.4g%14.4g%14.4g%14.4g\n",\
			runsCount,realtime,p.N,p.NT,p.L,p.Tb,p.DE,p.Theta,Num,E,imag(action),W);
	fprintf(cof,"\n");
	if (verbose) {
		printf("\n");
		printf("%8s%8s%8s%8s%8s%8s%8s%8s%14s%14s%14s%14s\n","runs","time","p.N","NT","L","Tb","dE",\
			"theta","Num","E","im(action)","W");
		printf("%8i%8g%8i%8i%8g%8g%8g%8g%14.4g%14.4g%14.4g%14.4g\n",\
			runsCount,realtime,p.N,p.NT,p.L,p.Tb,p.DE,p.Theta,Num,E,imag(action),W);
		printf("\n");
	}
	
	// printing results to file
	if ((checkDelta.good() && checkSol.good() && checkSolMax.good()) || pass) {
		// printing good results to file
		
		// id
		vector<string> idResult(idSizeResults);
		idResult[0] = timenumber;
		idResult[1] = nts(pl);
		
		// actual results
		vector<number> datumResult(datumSizeResults);
		datumResult[0] = Num;
		datumResult[1] = E;
		datumResult[2] = W;
		datumResult[3] = imag(action);
		datumResult[4] = imag(kineticT);
		datumResult[5] = imag(kineticS);
		datumResult[6] = imag(potV);
		datumResult[7] = imag(potr);
		datumResult[8] = imag(bound);
		datumResult[9] = checkSol.back();
		datumResult[10] = checkSolMax.back();
		datumResult[11] = checkLatt.back();
		datumResult[12] = checkLin.back();
		datumResult[13] = checkTrue.back();
		datumResult[14] = checkOS.back();
		datumResult[15] = checkAB.back();
		datumResult[16] = checkABNE.back();
		datumResult[17] = checkCon.back();
		datumResult[18] = checkIE.back();
		datumResult[19] = checkReg.back();
		
		// saving
		NewtonRaphsonDatum result(idResult,p,datumResult);
		result.save(resultsFile);
		fprintf(cof,"%12s%30s\n","results:",resultsFile.c_str());
		if (verbose)
			printf("%12s%24s\n","results:",resultsFile.c_str());
	}
	else {
		// printing error results to file	
		
		// id
		vector<string> idError(idSizeErrors);
		idError[0] = timenumber;
		idError[1] = nts(pl);
		idError[2] = nts(runsCount);
		
		// actual results
		vector<number> datumError(datumSizeErrors);
			datumError[0] = Num;
		datumError[1] = E;
		datumError[2] = W;
		datumError[3] = imag(action);
		datumError[4] = imag(kineticT);
		datumError[5] = imag(kineticS);
		datumError[6] = imag(potV);
		datumError[7] = imag(potr);
		datumError[8] = imag(bound);
		datumError[9] = checkSol.back();
		datumError[10] = checkSolMax.back();
		datumError[11] = checkLatt.back();
		datumError[12] = checkLin.back();
		datumError[13] = checkTrue.back();
		datumError[14] = checkOS.back();
		datumError[15] = checkAB.back();
		datumError[16] = checkABNE.back();
		datumError[17] = checkCon.back();
		datumError[18] = checkIE.back();
		datumError[19] = checkReg.back();
		datumError[20] = checkAction.back();
		datumError[21] = checkDelta.back();
		datumError[22] = checkInv.back();
		datumError[23] = checkReg.back();
		datumError[24] = checkContm.back();
		datumError[25] = checkLR.back();
		datumError[26] = checkBoundRe.back();
		datumError[27] = checkBoundIm.back();
		datumError[28] = checkChiT.back();
		
		// saving
		NewtonRaphsonDatum error(idError,p,datumError);
		error.save(errorsFile);
		fprintf(cof,"%12s%30s\n","results:",resultsFile.c_str());
		if (verbose)
			printf("%12s%24s\n","results:",errorsFile.c_str());
	}
	
	// printing field and its energy
	if (checkDelta.good() && checkSol.good() && checkSolMax.good()) {		
		// printing f to file
		Filename fRes = filenameMain(p,baseFolder,"field","fMain",suffix);
		saveVectorBinary(fRes,f);
		Filename ergRes = filenameMain(p,baseFolder,"energy","ergMain",suffix);
		saveComplexVectorBinary(ergRes,erg);
		Filename linErgRes = filenameMain(p,baseFolder,"energy","linErgMain",suffix);
		saveComplexVectorBinary(linErgRes,linErg);
		
		fprintf(cof,"%12s%50s\n","f     :",((string)fRes).c_str());
		fprintf(cof,"%12s%50s\n","erg   :",((string)ergRes).c_str());
		fprintf(cof,"%12s%50s\n","linErg:",((string)linErgRes).c_str());
		
		if (verbose) {
			printf("%12s%50s\n","f     :",((string)fRes).c_str());
			printf("%12s%50s\n","erg   :",((string)ergRes).c_str());
			printf("%12s%50s\n","linErg:",((string)linErgRes).c_str());
		}
		
	}
	
	if (!checkDelta.good() && !pass) {
			return 1;
		}
	fprintf(cof,"\n----------------------------------------------------------------------------------------------------------------------------\n\n");
	fclose(cof);
	} //ending parameter loop
ces.close();
cout << endl << "====================================================================================================================" << endl;

return 0;
}



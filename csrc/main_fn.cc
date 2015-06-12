/*----------------------------------------------------------------------------------------------------------------------------
	main_fn
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
#include "main_fn.h"
#include "main.h"

//#define NDEBUG //NDEBUG is to remove error and bounds checking on vectors in SparseLU, for speed - only include once everything works

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - argv inputs, loading options, closenesses
		2 - Folders
		3 - beginning file loop
		4 - beginning parameter loop
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

int main_fn(int argc, vector<string> argv)
{
/*----------------------------------------------------------------------------------------------------------------------------
	1. loading options, argv inputs
		- defining timenumber and files to load
		- argv inputs
		- loading options
		- loading closenesses
		- beginning cos and ces 
----------------------------------------------------------------------------------------------------------------------------*/

// defining timenumber
string timenumber = currentDateTime();

// defining rank
string rank = "0";

// options to load
Options opts;
Closenesses closenesses;
string optionsFile = "optionsM";
string closenessesFile = "closenesses";

// cos and cerr files
string coFile, ceFile;

// getting argv inputs
if (argc==2) timenumber = argv[1];
else if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("opts")==0 || id.compare("options")==0) optionsFile = argv[2*j+2];
		else if (id.compare("close")==0 || id.compare("closenesses")==0) closenessesFile = argv[2*j+2];
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
closenesses.load(closenessesFile);

// getting argv inputs
if (argc==2) timenumber = argv[1];
else if (argc % 2 && argc>1) {
	for (uint j=0; j<(uint)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("tn")==0 || id.compare("timenumber")==0) timenumber = argv[2*j+2];
		else if (id.compare("co")==0) coFile = (string)argv[2*j+2];
		else if (id.compare("ce")==0) ceFile = (string)argv[2*j+2];
		else if (id.compare("opts")==0 || id.compare("options")==0);
		else if (id.compare("close")==0 || id.compare("closenesses")==0);
		else if (id.compare("amp")==0) opts.amp = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("open")==0) opts.open = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("alpha")==0) opts.alpha = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("zmx")==0) opts.zmx = argv[2*j+2];
		else if (id.compare("zmt")==0) opts.zmt = argv[2*j+2];
		else if (id.compare("bds")==0) opts.bds = argv[2*j+2];
		else if (id.compare("inF")==0) opts.inF = argv[2*j+2];
		else if (id.compare("minTimenumberLoad")==0 || id.compare("mintn")==0) opts.minTimenumberLoad = argv[2*j+2];
		else if (id.compare("maxTimenumberLoad")==0 || id.compare("maxtn")==0) opts.maxTimenumberLoad = argv[2*j+2];
		else if (id.compare("minfLoopLoad")==0 || id.compare("minfLoopLoad")==0) opts.minfLoopLoad = argv[2*j+2];
		else if (id.compare("maxfLoopLoad")==0 || id.compare("maxfLoopLoad")==0) opts.maxfLoopLoad = argv[2*j+2];
		else if (id.compare("minLoopLoad")==0 || id.compare("minll")==0) opts.minLoopLoad = argv[2*j+2];
		else if (id.compare("maxLoopLoad")==0 || id.compare("maxll")==0) opts.maxLoopLoad = argv[2*j+2];
		else if (id.compare("loopChoice")==0) opts.loopChoice = argv[2*j+2];
		else if (id.compare("loopMin")==0) opts.loopMin = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("loopMax")==0) opts.loopMax = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("epsiTb")==0) opts.epsiTb = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("epsiTheta")==0) opts.epsiTheta = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("loops")==0) opts.loops = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("printChoice")==0) opts.printChoice = argv[2*j+2];
		else if (id.compare("rank")==0) rank = argv[2*j+2];
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}

// filling empty co and ce filenames
if (coFile.empty()) coFile = "data/"+timenumber+"co.txt";
if (ceFile.empty()) ceFile = "data/"+timenumber+"ce.txt";

// beginning cos and ces streams
fstream cos;
cos.open(coFile.c_str(),fstream::app);

fstream ces;
ces.open(ceFile.c_str(),fstream::app);

/*----------------------------------------------------------------------------------------------------------------------------
	2. Folders
		- FilenameAttributes for defining FilenameComparator
		- FilenameComparator
		- pFolder
		- inputsFolder
		- removeUnshared
		- printing folders
		
----------------------------------------------------------------------------------------------------------------------------*/

// FilenameAttributes for defining FilenameComparator
FilenameAttributes fa_low, fa_high;
fa_low.Directory = "data";
fa_high.Directory = "data";
fa_low.Timenumber = opts.minTimenumberLoad;
fa_high.Timenumber = opts.minTimenumberLoad;
(fa_low.Extras).push_back(StringPair("fLoop",opts.minfLoopLoad));
(fa_high.Extras).push_back(StringPair("fLoop",opts.maxfLoopLoad));
(fa_low.Extras).push_back(StringPair("loop",opts.minLoopLoad));
(fa_high.Extras).push_back(StringPair("loop",opts.maxLoopLoad));
if ((opts.inF).size()>1 && ((opts.loopChoice).substr(0,5)).compare("const")==0) {
	if ((opts.inF)[1]=='n') {
		(fa_low.Extras).push_back(StringPair("step","1"));
		(fa_high.Extras).push_back(StringPair("step","1"));
	}
}

// FilenameComparator
FilenameComparator fc(fa_low,fa_high);
Folder allFiles(fc);

// pFolder
if ((opts.inF)[0]=='p') {
	fa_low.ID = "tp";
	fa_high.ID = "tp";
}
else if ((opts.inF)[0]=='m') {
	fa_low.ID = "mainp";
	fa_high.ID = "mainp";
}
else {
	ces << "inF error: " << opts.inF << " not recognised" << endl;
	return 1;
}
if ((opts.inF).size()>1) {
	if ((opts.inF)[1]=='n') {
		fa_low.Suffix = ".data";
		fa_high.Suffix = ".data";
	}
	else {
		fa_low.Suffix = ".dat";
		fa_high.Suffix = ".dat";
	}
}
else {
	fa_low.Suffix = ".dat";
	fa_high.Suffix = ".dat";
}
fc.set(fa_low,fa_high);
Folder pFolder(fc);

// inputsFolder
if ((opts.inF)[0]=='p') {
	fa_low.ID = "inputsP";
	fa_high.ID = "inputsP";
}
else if ((opts.inF)[0]=='m') {
	fa_low.ID = "inputsM";
	fa_high.ID = "inputsM";
}
(fa_low.Extras).resize(fa_low.Extras.size()-1);
(fa_high.Extras).resize(fa_high.Extras.size()-1);
fa_low.Suffix = "";
fa_high.Suffix = "";
fc.set(fa_low,fa_high);
Folder inputsFolder(fc);

// removeUnshared - not quite working, should be based soley on timenumber and loop
// removeUnshared(pFolder,inputsFolder);

// printing folders
if (pFolder.size()>0 && inputsFolder.size()>0) {
	cos << endl;
	cos << "inputs: " << endl << pFolder << inputsFolder << endl;
}
else {
	ces << endl;
	ces << "not files found for options:" << endl;
	ces << opts << endl;
	return 1;
}

/*----------------------------------------------------------------------------------------------------------------------------
	3. beginning file loop
		- beginning file loop
		- loading and printing inputs
----------------------------------------------------------------------------------------------------------------------------*/

// beginning file loop
for (uint fileLoop=0; fileLoop<pFolder.size(); fileLoop++) {

	// loading parameters
	Parameters psu;
	psu.load(inputsFolder[fileLoop]);
	//cos << "input parameters: " << endl;
	//psu.print();
	//cos << endl;

/*----------------------------------------------------------------------------------------------------------------------------
	4. beginning parameter loop
		- initializing stepper
		- defining a time
		- changing parameters (if required)
		- copying a verson of parameters with timenumber
		- printing timenumber
		- declaring checks
----------------------------------------------------------------------------------------------------------------------------*/

	// initializing stepper
	StepperOptions step_opts;
	Point2d point;
	if (((opts.loopChoice).substr(0,5)).compare("const")==0) {
		step_opts.epsi_x = opts.epsiTb;
		step_opts.epsi_y = opts.epsiTheta;
		step_opts.angle0 = pi/2.0;
		step_opts.closeness = closenesses.Step;
		step_opts.stepType = StepperOptions::constPlane;
		step_opts.directed = StepperOptions::local;
		point(psu.Tb,psu.theta);
	}
	else {
		step_opts.epsi_x =  (opts.loopMax - opts.loopMin)/(opts.loops-1.0);;
		step_opts.epsi_y = 0.0;
		step_opts.angle0 = 0.0;
		step_opts.stepType = StepperOptions::straight;
		step_opts.directed = StepperOptions::undirected;
		step_opts.closeness = closenesses.Step; // irrelevant here
		point(opts.loopMin,0.0);
	}
	Stepper stepper(step_opts,point);
	
	bool stepped = true;
	
	for (uint loop=0; loop<opts.loops; loop++) {
	
		//defining a time and starting the clock
		clock_t time;
		time = clock();
	
		// changing parameters
		Parameters ps = psu;
		if (opts.loops>1) {
			if ((opts.loopChoice)[0]=='N') {
				bool anythingChanged = ps.changeParameters(opts.loopChoice,(uint)stepper.x());
				if (loop==0 && anythingChanged) {
					cos << opts.loopChoice << " changed to " << (uint)stepper.x() << " on input" << endl;
				}
			}
			else if (((opts.loopChoice).substr(0,5)).compare("const")!=0) {
				bool anythingChanged = ps.changeParameters(opts.loopChoice,stepper.x());
				if (loop==0 && anythingChanged) {
					cos << opts.loopChoice << " changed to " << stepper.x() << " on input" << endl;
				}
			}
			else {
				ps.changeParameters("Tb",stepper.x());
				ps.changeParameters("theta",stepper.y());
			}
		}

		//copying a version of ps with timenumber
		/*Filename paramsRunFile = (string)("data/"+timenumber+"inputsM_fLoop_"+numberToString<uint>(fileLoop)\
				+"_loop_"+numberToString<uint>(loop));
		ps.save(paramsRunFile);*/
		
		//printing timenumber
		cos.close();
		FILE* cof;
		cof = fopen(coFile.c_str(),"a");
		fprintf(cof,"%12s%12s\n","timenumber: ",timenumber.c_str());
		cout << "timenumber: " << timenumber << endl;
		if (((opts.loopChoice).substr(0,5)).compare("const")==0 && loop>0) {
			double angleModTwoPi = mod(stepper.stepAngle(),-pi,pi);		
			fprintf(cof,"%12s%12.3g\n","step angle: ",angleModTwoPi);
		}
		
		// declaring Checks
		Check checkAction("action",closenesses.Action);
		Check checkSoln("solution",closenesses.Soln);
		Check checkSolnMax("solution max",closenesses.SolnMax);
		Check checkDelta("delta",closenesses.Delta);
		Check checkInv("matrix inversion",closenesses.Inv*ps.N*ps.NT);
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
	
		// do trivial or redundant checks?
		bool trivialChecks = false;
	
/*----------------------------------------------------------------------------------------------------------------------------
	5. assigning potential functions
		- assigning potential functions
		- assigning preliminary parameter structs
----------------------------------------------------------------------------------------------------------------------------*/

		// assigning potential functions
		Potential<comp> V, dV, ddV;
		if (ps.pot==1) {
			V((Potential<comp>::PotentialType)&V1<comp>,ps);
			dV((Potential<comp>::PotentialType)&dV1<comp>,ps);
			ddV((Potential<comp>::PotentialType)&ddV1<comp>,ps);
		}
		else if (ps.pot==2) {
			V((Potential<comp>::PotentialType)&V2<comp>,ps);
			dV((Potential<comp>::PotentialType)&dV2<comp>,ps);
			ddV((Potential<comp>::PotentialType)&ddV2<comp>,ps);
		}
		else if (ps.pot==3) {
			V((Potential<comp>::PotentialType)&V3<comp>,ps);
			dV((Potential<comp>::PotentialType)&dV3<comp>,ps);
			ddV((Potential<comp>::PotentialType)&ddV3<comp>,ps);
			}
		else {
			ces << "pot option not available, pot = " << ps.pot << endl;
			return 1;
		}
	
		// assigning preliminary parameter structs
		params_for_V paramsV  = {ps.epsilon, ps.A};
	
		//lambda functions for pot_r
		auto Vr = [&] (const comp& phi) {
			return -ii*ps.reg*VrFn(phi,ps.minima[0],ps.minima[1]);
		};
		auto dVr = [&] (const comp& phi) {
			return -ii*ps.reg*dVrFn(phi,ps.minima[0],ps.minima[1]);
		};
		auto ddVr = [&] (const comp& phi) {
			return -ii*ps.reg*ddVrFn(phi,ps.minima[0],ps.minima[1]);
		};

/*----------------------------------------------------------------------------------------------------------------------------
	6. omega and negVec
		- loading or constructing omega
		- loading negVec
----------------------------------------------------------------------------------------------------------------------------*/

		//deterimining omega matrices for fourier transforms in spatial direction
		vec freqs(ps.N), freqs_exp(ps.N);
		mat modes(ps.N,ps.N);
		mat omega_m1(ps.N,ps.N), omega_0(ps.N,ps.N), omega_1(ps.N,ps.N), omega_2(ps.N,ps.N);
		SaveOptions so_simple;
		so_simple.printType = SaveOptions::binary;
		so_simple.paramsIn = ps; so_simple.paramsOut = ps;
		so_simple.vectorType = SaveOptions::simple;
		so_simple.extras = SaveOptions::none;
		so_simple.printMessage = false;
		{
			Filename omegaM1F, omega0F, omega1F, omega2F, modesF, freqsF, freqsExpF; // Filename works as FilenameAttributes
			omegaM1F = (string)"data/00"+rank+"omegaM1_pot_"+numberToString<uint>(ps.pot)+"_N_"+numberToString<uint>(ps.N)\
							+"_L_"+numberToString<double>(ps.L)+".data";
			Folder omegaM1Folder(omegaM1F);
			omega0F = omegaM1F; 		omega0F.ID = "omega0"; 		Folder omega0Folder(omega0F);
			omega1F = omegaM1F; 		omega1F.ID = "omega1"; 		Folder omega1Folder(omega1F);
			omega2F = omegaM1F; 		omega2F.ID = "omega2";	 	Folder omega2Folder(omega2F);
			modesF = omegaM1F;			modesF.ID = "modes"; 		Folder modesFolder(modesF);
			freqsF = omegaM1F;			freqsF.ID = "freqs"; 		Folder freqsFolder(freqsF);
			freqsExpF = omegaM1F;		freqsExpF.ID = "freqsExp";	Folder freqsExpFolder(freqsExpF);
			if (omegaM1Folder.size()==1 && omega0Folder.size()==1 && omega1Folder.size()==1 && omega2Folder.size()==1 \
				&& modesFolder.size()==1 && freqsFolder.size()==1 && freqsExpFolder.size()==1) {
				load(omegaM1Folder[0],so_simple,omega_m1);
				load(omega0Folder[0],so_simple,omega_0);
				load(omega1Folder[0],so_simple,omega_1);
				load(omega2Folder[0],so_simple,omega_2);
				load(modesFolder[0],so_simple,modes);
				load(freqsFolder[0],so_simple,freqs);
				load(freqsExpFolder[0],so_simple,freqs_exp);
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
				save(omegaM1F,so_simple,omega_m1);
				save(omega0F,so_simple,omega_0);
				save(omega1F,so_simple,omega_1);
				save(omega2F,so_simple,omega_2);
				save(modesF,so_simple,modes);
				save(freqsF,so_simple,freqs);
				save(freqsExpF,so_simple,freqs_exp);
			}
		}

		// loading negVec
		vec negVec;
		if (opts.zmt[0]=='n' || opts.zmx[0]=='n') {
			if (ps.pot==3) {
				Filename eigVecFile = "data/00"+rank+"eigVec_pot_3_L_" + numberToString<double>(ps.L) + ".dat";	
				so_simple.printType = SaveOptions::ascii;
				load(eigVecFile,so_simple,negVec); // should automatically interpolate
				so_simple.printType = SaveOptions::binary;
			}
			else {
				Filename eigVecFile;
				string N_load;
				string Nb_load;
				Filename lower = "data/00"+rank+"eigVec_pot_" + numberToString<uint>(ps.pot)\
									 + "_N_100_Nb_100_L_" + numberToString<double>(ps.L) + ".dat";
				Filename upper = lower;
				vector<StringPair> upperExtras = lower.Extras;
				upper.Extras[1] = StringPair("N","1000");
				upper.Extras[2] = StringPair("Nb","1000");
				Folder eigVecFolder(lower,upper);
				if (eigVecFolder.size()==0) {
					ces << "no negative eigenvector files found between:" << endl << lower << endl << upper << endl;
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
				eigVecOpts.paramsOut = ps;
				eigVecOpts.printMessage = false;
				Parameters pIn = ps;
				pIn.N = stringToNumber<uint>(N_load);
				pIn.Nb = stringToNumber<uint>(Nb_load);
				pIn.NT = 1000; // a fudge so that interpolate realises the vector is only on BC
				load(eigVecFile,eigVecOpts,negVec);
			}
		}

/*----------------------------------------------------------------------------------------------------------------------------
	7. defining quantities
		- erg, linErg etc
		- p, minusDS, DDS
		- printing parameters
		- print input phi
----------------------------------------------------------------------------------------------------------------------------*/	
	
		//defining energy and number vectors
		cVec erg(ps.NT);
		cVec linErg(ps.NT);
		cVec linNum(ps.NT);
		cVec linErgOffShell(ps.NT);
		cVec linNumOffShell(ps.NT);
		cVec derivErg(ps.NT), potErg(ps.NT);
		comp linErgContm, linNumContm;
		comp linNumAB, linErgAB;
		
		//defining the action and bound and W and zero of energy
		double ergZero = (ps.pot==3? 0.0: ps.N*ps.a*real(V(ps.minima[0])) );
		comp action = ii*ps.action0;
		double bound = 0.0;
		double W;
		double E;
		double E_exact;
		double Num;
		
		//defining some quantities used to stop the Newton-Raphson loop when action stops varying
		comp action_last = action;
		uint runs_count = 0;
		uint min_runs = 3;
		uint max_runs = 100;

		//initializing phi (=p)
		vec p;
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
			load(pFolder[0],so_tp,p); //n.b there may be some problems with zero modes for binary printing
			so_tp.printType = SaveOptions::binary;
			fprintf(cof,"%12s%30s\n","input: ",(pFolder[0]()).c_str());
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
			so_tp.paramsIn = ps;
			load(lastPhi,so_tp,p);
			fprintf(cof,"%12s%30s\n","input: ",(lastPhi()).c_str());
		}
		
		// printing parameters
		ps.print(cof);
		
		//defining complexified vector Cp
		cVec Cp;
		Cp = vecComplex(p,ps);
	
		//defining DDS and minusDS
		spMat DDS(2*ps.N*ps.NT+2,2*ps.N*ps.NT+2);
		vec minusDS(2*ps.N*ps.NT+2);
			
		//very early vector print
		if ((opts.printChoice).compare("n")!=0) {
			so_tp.paramsIn = ps;
			Filename earlyPrintFile = (string)("data/"+timenumber+"mainpE_fLoop_"+numberToString<uint>(fileLoop)\
					 +"_loop_"+numberToString<uint>(loop)+"_run_" + "0.data");
			save(earlyPrintFile,so_tp,p);
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
				cerr << "main error: max_runs(" << max_runs << ") exceeded for:" << endl;
				cerr << "timenumber: " << timenumber << "; fileLoop: " << fileLoop << "; loop: " << loop << endl;
				return 1;				
			}
			
			// zero modes - fixed with chiX and chiT
			vec chiX(ps.NT*ps.N);	chiX = Eigen::VectorXd::Zero(ps.N*ps.NT); //to fix spatial zero mode
			vec chiT(ps.NT*ps.N);	chiT = Eigen::VectorXd::Zero(ps.N*ps.NT); //to fix real time zero mode
			for (uint j=0; j<ps.N; j++) {
				uint posX, posT, posCe;
				uint slicesX, slicesT;
				if (getLastInt(opts.zmx)<0) {
					ces << "getLastInt error with zmx = " << opts.zmx << endl;
					return 1;
				}
				if (getLastInt(opts.zmt)<0) {
					ces << "getLastInt error with zmt = " << opts.zmt << endl;
					return 1;
				}
				slicesX = getLastInt(opts.zmx);
				slicesT = getLastInt(opts.zmt);
				posCe = j*ps.Nb+ps.Nb-slicesT; //position C for Euclidean vector, i.e. for negVec
				map<char,unsigned int> posMap;
				posMap['A'] = j*ps.NT;
				posMap['B'] = j*ps.NT + ps.Na-1;
				posMap['C'] = j*ps.NT+ps.Na+ps.Nb-slicesT;
				posMap['D'] = j*ps.NT+(ps.NT-1)-slicesT;
				if ((opts.zmt).size()<3) ces << "zmt lacks info, zmt = " << opts.zmt << endl;
				for (uint l=0;l<((opts.zmt).size()-2);l++) {
					if (posMap.find(opts.zmt[1+l])!=posMap.end()) {
						posT = posMap.at(opts.zmt[1+l]);
						for (uint k=0;k<slicesT;k++) {
							if (opts.zmt[0]=='n' && ps.pot!=3) 	chiT(posT+k) = negVec(2*(posCe+k));
							else if (opts.zmt[0]=='n' && ps.pot==3) {
													double r = ps.r0 + j*ps.a;
													chiT(posT+k) = negVec(j)*r;
							}
							else if (opts.zmt[0]=='d')	chiT(posT+k) = p(2*(posT+k+1))-p(2*(posT+k));
							else {
								ces << "choice of zmt(" << opts.zmt << ") not allowed" << endl;
								return 1;
							}
						}
					}
				}
				posMap.erase('C');
				posMap.erase('D');
				posMap['C'] = j*ps.NT+ps.Na+ps.Nb-slicesX;
				posMap['D'] = j*ps.NT+ps.NT-slicesX;
				posCe = j*ps.Nb+ps.Nb-slicesX;
				if ((opts.zmx).size()<3) ces << "zmx lacks info, zmx = " << opts.zmx << endl;
				for (uint l=0;l<((opts.zmx).size()-2);l++) {
					if (posMap.find(opts.zmx[1+l])!=posMap.end()) {
						posX = posMap.at(opts.zmx[1+l]);
						for (uint k=0;k<slicesX;k++) {
							if (opts.zmx[0]=='n' && ps.pot!=3)		chiX(posX+k) = negVec(2*(posCe+k));
							else if (opts.zmx[0]=='n' && ps.pot==3) {
																double r = ps.r0 + j*ps.a;
																chiX(posX+k) = negVec(j)*r;
							}
							else if (opts.zmx[0]=='d' && ps.pot!=3)
								chiX(posX+k) = p(2*neigh(posX+k,1,1,ps))-p(2*neigh(posX+k,1,-1,ps));
							else {
								ces << "choice of zmx(" << opts.zmx << ") not allowed" << endl;
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
			}
			chiX = chiX/normX;
			chiT = chiT/normT;
			
			// allocating memory for DS, DDS
			minusDS = Eigen::VectorXd::Zero(2*ps.N*ps.NT+2); //initializing to zero
			DDS.setZero(); //just making sure
			Eigen::VectorXi DDS_to_reserve(2*ps.N*ps.NT+2);//number of non-zero elements per column
			DDS_to_reserve = Eigen::VectorXi::Constant(2*ps.N*ps.NT+2,13);
			if (abs(ps.theta)<2.0e-16) {
				DDS_to_reserve(0) = ps.N+3;
				DDS_to_reserve(1) = 11;
			}
			else {
				DDS_to_reserve(0) = ps.N+11;
				DDS_to_reserve(1) = ps.N+11;
			}
			DDS_to_reserve(2*ps.N*ps.NT-2) = 4;
			DDS_to_reserve(2*ps.N*ps.NT-1) = 4;
			DDS_to_reserve(2*ps.N*ps.NT) = ps.N*((opts.zmx).size()-2);
			DDS_to_reserve(2*ps.N*ps.NT+1) = 2*ps.N*((opts.zmt).size()-2);
			DDS.reserve(DDS_to_reserve);
			
			//initializing to zero
			comp kineticS = 0.0;
			comp kineticT = 0.0;
			comp potV = 0.0;
			comp pot_r = 0.0;
			bound = 0.0;
			erg = Eigen::VectorXcd::Constant(ps.NT,-ergZero);
			linErg = Eigen::VectorXcd::Zero(ps.NT);
			linErgOffShell = Eigen::VectorXcd::Zero(ps.NT);
			linNum = Eigen::VectorXcd::Zero(ps.NT);
			linNumOffShell = Eigen::VectorXcd::Zero(ps.NT);
			derivErg = Eigen::VectorXcd::Zero(ps.NT);
			potErg = Eigen::VectorXcd::Constant(ps.NT,-ergZero);
			linErgContm = 0.0;
			linNumContm = 0.0;
			
			//testing that the potential term is working for pot3
			if (ps.pot==3 && trivialChecks) {
				comp Vtrial = 0.0, Vcontrol = 0.0;
				for (unsigned int j=0; j<ps.N; j++) {
					double r = ps.r0 + j*ps.a;
					paramsV.epsi = r;
					V.setParams(paramsV);
					Vcontrol += pow(p(2*j*ps.Nb),2.0)/2.0 - pow(p(2*j*ps.Nb),4.0)/4.0/pow(r,2.0);
					Vtrial += V(p(2*j*ps.Nb));
				}
				double potTest = pow(pow(real(Vcontrol-Vtrial),2.0) + pow(imag(Vcontrol-Vtrial),2.0),0.5);
				fprintf(cof,"potTest = %8.4g\n",potTest);
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
			for (lint j = 0; j < ps.N*ps.NT; j++) {		
				uint t 					= intCoord(j,0,ps); //coordinates
				uint x 					= intCoord(j,1,ps);
				int neighPosX 			= neigh(j,1,1,ps);
				int neighNegX  			= neigh(j,1,-1,ps);
				
				comp Dt 		= DtFn(t,ps);
				comp dt 		= dtFn(t,ps);
				comp dtm 		= (t>0? dtFn(t-1,ps): ps.b);
				double Dx 		= DxFn(x,ps);
				double dx 		= dxFn(x,ps);
				double dxm 		= (x>0? dxFn(x-1,ps): ps.a);
				
				if (ps.pot==3) {
					paramsV.epsi = ps.r0+x*ps.a;
					V.setParams(paramsV);
					dV.setParams(paramsV);
					ddV.setParams(paramsV);
				}
			
				if (abs(chiX(j))>MIN_NUMBER && ps.pot!=3) { //spatial zero mode lagrange constraint
					DDS.insert(2*j,2*ps.N*ps.NT) 			= Dx*chiX(j); 
					DDS.insert(2*ps.N*ps.NT,2*j) 			= Dx*chiX(j);
					minusDS(2*j) 					+= -Dx*chiX(j)*p(2*ps.N*ps.NT);
					minusDS(2*ps.N*ps.NT) 				+= -Dx*chiX(j)*p(2*j);
				}
					
				if (abs(chiT(j))>MIN_NUMBER && t<(ps.NT-1)) {
					DDS.coeffRef(2*(j+1),2*ps.N*ps.NT+1) 	+= Dx*chiT(j); //chiT should be 0 at t=(ps.NT-1) or this line will go wrong
					DDS.coeffRef(2*ps.N*ps.NT+1,2*(j+1)) 	+= Dx*chiT(j);
					DDS.coeffRef(2*j,2*ps.N*ps.NT+1) 		+= -Dx*chiT(j);
					DDS.coeffRef(2*ps.N*ps.NT+1,2*j) 		+= -Dx*chiT(j);
		            minusDS(2*(j+1)) 				+= -Dx*chiT(j)*p(2*ps.N*ps.NT+1);
		            minusDS(2*j) 					+= Dx*chiT(j)*p(2*ps.N*ps.NT+1);
		            minusDS(2*ps.N*ps.NT+1) 				+= -Dx*chiT(j)*(p(2*(j+1))-p(2*j));
				}
					
				if (t<(ps.NT-1)) {
					for (uint k=0;k<ps.N;k++) {
						lint l = k*ps.NT+t;
						linErgOffShell(t) += 0.5*( omega_2(x,k)*( Cp(l)-ps.minima[0] )*( Cp(j)-ps.minima[0] ) \
										+ omega_0(x,k)*( Cp(l+1)-Cp(l) )*( Cp(j+1)-Cp(j) )/pow(dt,2.0));
						linNumOffShell(t) += 0.5*(omega_1(x,k)*( Cp(l)-ps.minima[0] )*( Cp(j)-ps.minima[0] ) \
										+ omega_m1(x,k)*( Cp(l+1)-Cp(l) )*( Cp(j+1)-Cp(j) )/pow(dt,2.0));
					}
				}
					
				if (abs(ps.theta)>MIN_NUMBER) {
					for (uint k=0;k<ps.N;k++) {
						lint l = k*ps.NT+t;
						linNum(t) += 2.0*ps.Gamma*omega_1(x,k)*( (p(2*l)-ps.minima[0])*(p(2*j)-ps.minima[0])/pow(1.0+ps.Gamma,2.0)\
								 + p(2*j+1)*p(2*l+1)/pow(1.0-ps.Gamma,2.0) );
						linErg(t) += 2.0*ps.Gamma*omega_2(x,k)*( (p(2*l)-ps.minima[0])*(p(2*j)-ps.minima[0])/pow(1.0+ps.Gamma,2.0)\
								+ p(2*j+1)*p(2*l+1)/pow(1.0-ps.Gamma,2.0) );
					}
				}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//boundaries
				if (ps.pot==3 && x==(ps.N-1)) {
					DDS.insert(2*j,2*j) 	= 1.0; // p=0 at r=R
					DDS.insert(2*j+1,2*j+1) = 1.0;
				}
				else if (ps.pot==3 && x==0) {
					kineticS 				+= 	Dt*pow(Cp(neighPosX),2.0)/dx/2.0;
					derivErg(t) 			+= 	pow(Cp(neighPosX),2.0)/dx/2.0; //n.b. the Dt/dt difference is ignored for erg(t)
					erg(t) 					+= 	pow(Cp(neighPosX),2.0)/dx/2.0;
					
					DDS.insert(2*j,2*j) 	= 1.0; // p=0 at r=0
					DDS.insert(2*j+1,2*j+1) = 1.0;
				}			
				else if (t==(ps.NT-1)) {
					if (neighPosX!=-1) {
						kineticS			+= Dt*pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
						derivErg(t) 		+= pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
						erg(t) 				+= pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
					}			
					potV 					+= Dt*Dx*V(Cp(j));
					pot_r 					+= Dt*Dx*Vr(Cp(j));
					erg(t) 					+= Dx*V(Cp(j)) + Dx*Vr(Cp(j));
					potErg(t) 				+= Dx*V(Cp(j)) + Dx*Vr(Cp(j));
				
					DDS.insert(2*j,2*(j-1)+1) = 1.0; //zero imaginary part of time derivative
					DDS.insert(2*j+1,2*j+1)   = 1.0; //zero imaginary part
				}
				else if (t==0) {
					kineticT 	+= Dx*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
					derivErg(t) += Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0;
					erg(t) 		+= Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0;
					
					if ((opts.bds).compare("jd")==0) {
						minusDS(2*j+1) 					+= Dx*imag(Cp(j+1)/dt);
					    DDS.coeffRef(2*j+1,2*(j+1)) 	+= -imag(Dx/dt);
					    DDS.coeffRef(2*j+1,2*(j+1)+1)	+= -real(Dx/dt);
					    minusDS(2*j+1) 					+= imag(-Dx*Cp(j)/dt);
					    DDS.coeffRef(2*j+1,2*j) 		+= imag(Dx/dt);
					    DDS.coeffRef(2*j+1,2*j+1) 		+= real(Dx/dt);
					    if (abs(ps.theta)<MIN_NUMBER) {
							/////////////////////////////////////equation R - theta=0//////////////////////////////////////
							for (uint k=0;k<ps.N;k++) {
								if (abs(omega_1(x,k))>MIN_NUMBER) {
									lint m=k*ps.NT;
									DDS.coeffRef(2*j,2*m+1) += -2.0*omega_1(x,k);
									minusDS(2*j) 			+= 2.0*omega_1(x,k)*p(2*m+1);
								}
							}
							////////////////////////////////////////////////////////////////////////////////////////
						}
						else {
							for (uint k=0;k<ps.N;k++) {
								if (abs(omega_1(x,k))>MIN_NUMBER) {
									/////////////////////equation I - theta!=0//////////////
									lint m=k*ps.NT;
									DDS.coeffRef(2*j+1,2*m) += (1.0-ps.Gamma)*omega_1(x,k)/(1.0+ps.Gamma);
									minusDS(2*j+1) 			+= -(1.0-ps.Gamma)*omega_1(x,k)*(p(2*m)-ps.minima[0])/(1.0+ps.Gamma);
									/////////////////////equation R - theta!=0//////////////
									minusDS(2*j) 			+= ps.theta*p(2*m+1)*omega_1(x,k)*(1+ps.Gamma)/(1-ps.Gamma);
									DDS.coeffRef(2*j,2*m+1)	+= -ps.theta*omega_1(x,k)*(1.0+ps.Gamma)/(1.0-ps.Gamma);
									bound 		+= -(1.0-ps.Gamma)*omega_1(x,k)*(p(2*j)-ps.minima[0])*(p(2*m)-ps.minima[0])/(1.0+ps.Gamma)\
																 + (1.0+ps.Gamma)*omega_1(x,k)*p(2*j+1)*p(2*m+1)/(1.0-ps.Gamma);
								}
							}
							//////////////////////////////////////equation R - theta!=0//////////////////////////////
				            minusDS(2*j) 					+= ps.theta*Dx*real(Cp(j+1)/dt);
			            	DDS.coeffRef(2*j,2*(j+1)) 		+= -ps.theta*real(Dx/dt);
			            	DDS.coeffRef(2*j,2*(j+1)+1) 	+= ps.theta*imag(Dx/dt);
						    minusDS(2*j) 					+= ps.theta*real(-Dx*Cp(j)/dt);
					    	DDS.coeffRef(2*j,2*j) 			+= ps.theta*real(Dx/dt);
					    	DDS.coeffRef(2*j,2*j+1) 		+= ps.theta*imag(-Dx/dt);
						}
					}
					else { ///////////////////////////////// including other terms in action at t=0 ///////////////////////////
						if (neighPosX!=-1) {
							kineticS 	+= Dt*pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
							derivErg(t) += pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
							erg(t) 		+= pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
						}
						potV 		+= Dt*Dx*V(Cp(j));
						pot_r 		+= Dt*Dx*Vr(Cp(j));
						potErg(t) 	+= Dx*V(Cp(j)) + Dx*Vr(Cp(j));
						erg(t) 		+= Dx*V(Cp(j)) + Dx*Vr(Cp(j));
						//////////////////////////////////////equation I - both///////////////////////////////////
						for (uint k=1; k<2*2; k++) {
					        int sign = pow(-1,k+1);
					        uint direc = (uint)(k/2.0);
					        int neighb = neigh(j,direc,sign,ps);
					        double dxd = (sign==1? dx: dxm);
					        if (direc == 0) {
					            minusDS(2*j+1) 					+= Dx*imag(Cp(j+sign)/dt);
					            DDS.coeffRef(2*j+1,2*(j+sign)) 	+= -imag(Dx/dt);
					            DDS.coeffRef(2*j+1,2*(j+sign)+1)+= -real(Dx/dt);
					           }
					        else if (neighb!=-1) {
					            minusDS(2*j+1) 					+= - imag(Dt*Cp(neighb))/dxd;
					            DDS.coeffRef(2*j+1,2*neighb) 	+= imag(Dt)/dxd;
					            DDS.coeffRef(2*j+1,2*neighb+1) 	+= real(Dt)/dxd;
					        }
					    }
					    comp temp0 = Dx/dt - Dt*(1.0/dx+1.0/dxm);
					    if (neighPosX==-1) 		temp0 += Dt/dx;
					    else if (neighNegX==-1) temp0 += Dt/dxm;
					    comp temp1 = Dt*Dx*( dV(Cp(j)) + dVr(Cp(j)) );//dV terms should be small
					    comp temp2 = Dt*Dx*(ddV(Cp(j)) + ddVr(Cp(j)));
					    
					    minusDS(2*j+1) 				+= imag(-temp0*Cp(j) + temp1 );
					    DDS.coeffRef(2*j+1,2*j) 	+= imag(temp0 - temp2 );
					    DDS.coeffRef(2*j+1,2*j+1) 	+= real(temp0 - temp2 );
						/////////////////////////////////////////////////////////////////////////////////////////
						if (abs(ps.theta)<MIN_NUMBER) {
							//simplest boundary conditions replaced by ones continuously connected to theta!=0 ones
							//DDS.insert(2*j+1,2*(j+1)+1) = 1.0; //zero imaginary part of time derivative
							//DDS.insert(2*j,2*j+1) = 1.0; //zero imaginary part
						
							/////////////////////////////////////equation R - theta=0//////////////////////////////////////
							for (uint k=0;k<ps.N;k++) {
								if (abs(omega_1(x,k))>MIN_NUMBER) {
									lint m=k*ps.NT;
									DDS.coeffRef(2*j,2*m+1) += -2.0*omega_1(x,k);
									minusDS(2*j) 			+= 2.0*omega_1(x,k)*p(2*m+1);
								}
							}
							////////////////////////////////////////////////////////////////////////////////////////
						}
						else {
							for (uint k=0;k<ps.N;k++) {
								if (abs(omega_1(x,k))>MIN_NUMBER) {
									/////////////////////equation I - theta!=0//////////////
									lint m=k*ps.NT;
									DDS.coeffRef(2*j+1,2*m) += (1.0-ps.Gamma)*omega_1(x,k)/(1.0+ps.Gamma);
									minusDS(2*j+1) 			+= -(1.0-ps.Gamma)*omega_1(x,k)*(p(2*m)-ps.minima[0])/(1.0+ps.Gamma);
									/////////////////////equation R - theta!=0//////////////
									minusDS(2*j) 			+= p(2*m+1)*omega_1(x,k)*(1+ps.Gamma)*ps.theta/(1-ps.Gamma);
									DDS.coeffRef(2*j,2*m+1)	+= -omega_1(x,k)*(1.0+ps.Gamma)*ps.theta/(1.0-ps.Gamma);
									bound 				+= -(1.0-ps.Gamma)*omega_1(x,k)*(p(2*j)-ps.minima[0])\
																*(p(2*m)-ps.minima[0])/(1.0+ps.Gamma)
																 + (1.0+ps.Gamma)*omega_1(x,k)*p(2*j+1)*p(2*m+1)/(1.0-ps.Gamma);
								}
							}
							//////////////////////////////////////equation R - theta!=0//////////////////////////////
							for (uint k=1; k<2*2; k++){
						        int sign = pow(-1,k+1);
						        uint direc = (uint)(k/2.0);
						        int neighb = neigh(j,direc,sign,ps);
						        double dxd = (sign==1? dx: dxm);
						        if (direc == 0) {
						            minusDS(2*j) 					+= Dx*real(Cp(j+sign)/dt)*ps.theta;
					            	DDS.coeffRef(2*j,2*(j+sign)) 	+= -real(Dx/dt)*ps.theta;
					            	DDS.coeffRef(2*j,2*(j+sign)+1) 	+= imag(Dx/dt)*ps.theta;
						        }
						        else if (neighb!=-1) {
						            minusDS(2*j) 					+= - real(Dt*Cp(neighb))*ps.theta/dxd;
						            DDS.coeffRef(2*j,2*neighb) 		+= real(Dt)*ps.theta/dxd;
					            	DDS.coeffRef(2*j,2*neighb+1) 	+= -imag(Dt)*ps.theta/dxd;
						        }
						    }
						    minusDS(2*j) 			+= ps.theta*real(-temp0*Cp(j) + temp1 );
					    	DDS.coeffRef(2*j,2*j) 	+= ps.theta*real(temp0 - temp2 );
					    	DDS.coeffRef(2*j,2*j+1) += ps.theta*imag(-temp0 + temp2 );
						}
					}
				}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//bulk
				else {
					if (neighPosX!=-1) {
						kineticS 	+= Dt*pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
						erg(t) 	 	+= pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
						derivErg(t) += pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
					}
					
					kineticT 	+= Dx*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
					potV 		+= Dt*Dx*V(Cp(j));
					pot_r 		+= Dt*Dx*Vr(Cp(j));
					erg(t) 		+= Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + Dx*V(Cp(j)) + Dx*Vr(Cp(j));
					derivErg(t) += Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0;
					potErg(t) 	+= Dx*V(Cp(j)) + Dx*Vr(Cp(j));
				
		            for (uint k=0; k<2*2; k++) {
		                int sign = pow(-1,k);
		                uint direc = (uint)(k/2.0);
		                int neighb = neigh(j,direc,sign,ps);
		                comp dtd = (sign==1? dt: dtm);
		                double dxd = (sign==1? dx: dxm);
		                if (direc == 0) {
		                    minusDS(2*j) 					+= real(Dx*Cp(j+sign)/dtd);
		                    minusDS(2*j+1) 					+= imag(Dx*Cp(j+sign)/dtd);
		                    DDS.insert(2*j,2*(j+sign)) 		= -real(Dx/dtd);
		                    DDS.insert(2*j,2*(j+sign)+1) 	= imag(Dx/dtd);
		                    DDS.insert(2*j+1,2*(j+sign)) 	= -imag(Dx/dtd);
		                    DDS.insert(2*j+1,2*(j+sign)+1) 	= -real(Dx/dtd);
		                }
		                else if (neighb!=-1) {                        
		                    minusDS(2*j) 					+= -real(Dt*Cp(neighb)/dxd);
		                    minusDS(2*j+1) 					+= -imag(Dt*Cp(neighb)/dxd);
		                    DDS.insert(2*j,2*neighb) 		= real(Dt/dxd);
		                    DDS.insert(2*j,2*neighb+1) 		= -imag(Dt/dxd);
		                    DDS.insert(2*j+1,2*neighb) 		= imag(Dt/dxd);
		                    DDS.insert(2*j+1,2*neighb+1) 	= real(Dt/dxd);
		                }
                	}
		            comp temp0 = Dx*(1.0/dt + 1.0/dtm) - Dt*(1.0/dx + 1.0/dxm);
		            if (neighPosX==-1) 		temp0 += Dt/dx;
		            else if (neighNegX==-1) temp0 += Dt/dxm;
		            comp temp1 = Dt*Dx*(dV(Cp(j)) + dVr(Cp(j)));
		            comp temp2 = Dt*Dx*(ddV(Cp(j)) + ddVr(Cp(j)));
		                
		            minusDS(2*j) 			+= real(temp1 - temp0*Cp(j));
		            minusDS(2*j+1) 			+= imag(temp1 - temp0*Cp(j));
		            DDS.insert(2*j,2*j) 	= real(-temp2 + temp0);
		            DDS.insert(2*j,2*j+1) 	= imag(temp2 - temp0);
		            DDS.insert(2*j+1,2*j) 	= imag(-temp2 + temp0);
		            DDS.insert(2*j+1,2*j+1) = real(-temp2 + temp0);
		        }
		    } // end of loop over j
		    
		    if (ps.pot==3) 		DDS.insert(2*ps.N*ps.NT,2*ps.N*ps.NT) = 1.0;
		    action = kineticT - kineticS - potV - pot_r;
		    linErgOffShell(ps.NT-1) = linErgOffShell(ps.NT-2);
	    	linNumOffShell(ps.NT-1) = linNumOffShell(ps.NT-2);
	    	
		    if (ps.pot==3) {
		    	action 		*= 4.0*pi;
		    	derivErg 	*= 4.0*pi;
		    	potErg 		*= 4.0*pi;
		    	erg			*= 4.0*pi;
		    	linErg		*= 4.0*pi;
		    	linNum		*= 4.0*pi;
		    	linNumOffShell *= 4.0*pi;
		    	linErgOffShell *= 4.0*pi;
		    }
		   
		    if (abs(ps.theta)<MIN_NUMBER) {
		    	linErg = linErgOffShell;
		    	linNum = linNumOffShell;
		    }
		    
		    //defining E, Num and W
			E = real(linErg(0));
			Num = real(linNum(0));
			W = - E*2.0*ps.Tb - ps.theta*Num - bound + 2.0*imag(action);
			
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
		
----------------------------------------------------------------------------------------------------------------------------*/
			
			//trivial test that erg=potErg+derivErg
			if (trivialChecks) {
				for (uint j=0; j<ps.NT; j++) {
					double diff = absDiff(erg(j),potErg(j)+derivErg(j));
					if (diff>1.0e-14) ces << "erg(" << j << ") != potErg + derivErg. absDiff = " << diff << endl;
				}
			}
			
			//calculating continuum approx to linErg and linNum on initial time slice - redundant
			if (ps.pot==3 && abs(ps.theta)<MIN_NUMBER) {
				for (uint k=1; k<ps.N; k++) {
					double momtm = k*pi/(ps.L-ps.r0);
					double freqSqrd = 1.0+pow(momtm,2.0);
					double Asqrd, integral1 = 0.0, integral2 = 0.0;
					for (unsigned int l=0; l<ps.N; l++) {
						double r = ps.r0 + l*ps.a;
						lint m = l*ps.NT;
						integral1 += ps.a*p(2*m)*pow(2.0/ps.L,0.5)*sin(momtm*r);
						integral2 += ps.a*(p(2*(m+1))-p(2*m))*pow(2.0/ps.L,0.5)*sin(momtm*r)/ps.b;
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
			cVec a_k(ps.N), b_k(ps.N); //ps.N.b. b_k means b*_k but can't put * in the name
			a_k = Eigen::VectorXcd::Zero(ps.N);
			b_k = Eigen::VectorXcd::Zero(ps.N);
			double T0 = 0.0;
			comp dt0 = dtFn(0,ps);
			for (uint n=0; n<ps.N; n++) {
				double w_n = freqs(n);
				double w_n_e = freqs_exp(n);
				for (uint j=0; j<ps.N; j++) {
					lint m=j*ps.NT;
					double sqrtDj = sqrt(4.0*pi*DxFn(j,ps));
					if (abs(w_n)>1.0e-16 && abs(w_n_e)>1.0e-16) {
						a_k(n) += exp(ii*w_n_e*T0)*sqrt(2.0*w_n)*modes(j,n)* \
									sqrtDj*((Cp(m+1)-ps.minima[0])-(Cp(m)-ps.minima[0])*exp(ii*w_n_e*dt0)) \
										/(exp(-ii*w_n_e*dt0)-exp(ii*w_n_e*dt0));
						b_k(n) += exp(-ii*w_n_e*T0)*sqrt(2.0*w_n)*modes(j,n)* \
									sqrtDj*((Cp(m+1)-ps.minima[0])-(Cp(m)-ps.minima[0])*exp(-ii*w_n_e*dt0)) \
										/(exp(ii*w_n_e*dt0)-exp(-ii*w_n_e*dt0));
					}
				}
			}			
			
			//using a_k and b*_k to check that a_k=Gamma*b*_k and p->p_lin as t->0 and that Sum_k(a_k*b*_k)=linNum
			// and that Sum_k(w_k*a_k*b*_k)=linErg
			double ABtest = 0.0, linRepTest;
			linNumAB = 0.0, linErgAB = 0.0;
			cVec linRep(ps.N), p0(ps.N);
			linRep = Eigen::VectorXcd::Constant(ps.N,ps.minima[0]);
			for (uint j=0; j<ps.N; j++) {
				ABtest += absDiff(a_k(j),conj(b_k(j))*ps.Gamma);
				p0(j) = Cp(j*ps.NT);
				linNumAB += a_k(j)*b_k(j);
				linErgAB += freqs(j)*a_k(j)*b_k(j);
				if (trivialChecks) {
					for (uint n=0; n<ps.N; n++) {
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
			ABtest /= (double)ps.N;
			checkAB.add(ABtest);
			double ABNtest = absDiff(linNumAB,linNum(0));
			double ABEtest = absDiff(linErgAB,linErg(0));
			ABNtest>ABEtest? checkABNE.add(ABNtest): checkABNE.add(ABEtest);
			
			//checking pot_r is much smaller than the other potential terms
			checkReg.add(abs(pot_r/potV));
			checkReg.checkMessage();
						
			//checking linearisation of linErg and linNum
			double linTestE;		double linTestN;
			double linEMax = 0.0;	double linEMin = 5.0e15; //surely it's going to be less than this
			double linNMax = 0.0;	double linNMin = 5.0e15;
			uint linearInt = (uint)(ps.Na/10);
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
			double conservTest = absDiff(erg(1),erg(ps.NT-2));
			checkCon.add(conservTest);
			
			//testing imaginary part of energy
			double imErgTest = 0.0;
			for (uint j=0; j<ps.NT; j++) if (abs(erg(j))>MIN_NUMBER) imErgTest += imag(erg(j))/abs(erg(j));
			imErgTest /= (double)ps.NT;
			checkIE.add(imErgTest);
			
			//checking agreement between erg and linErg
			E_exact = 0.0;
			for (uint j=1; j<(linearInt+1); j++) E_exact += real(erg(j));
			E_exact /= (double)linearInt;
			double trueTest = absDiff(E,E_exact);
			checkTrue.add(trueTest);
			if (!isfinite(trueTest))
				ces << "E = " << E << ", E_exact = " << E_exact << ", linearInt = " << linearInt << endl;
			
			//checking lattice small enough for E, should have parameter for this
			double momTest = E*ps.b/Num/pi; //perhaps should have a not b here
			checkLatt.add(momTest);
	
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
		- p' = p + delta
		- Cp'
----------------------------------------------------------------------------------------------------------------------------*/
			
			//solving for delta in DDS*delta=minusDS, where p' = p + delta		
			vec delta(2*ps.N*ps.NT+2);
			delta = Eigen::VectorXd::Zero(2*ps.N*ps.NT+2);
			DDS.prune(MIN_NUMBER);
			DDS.makeCompressed();
			Eigen::SparseLU<spMat> solver;
		
			solver.analyzePattern(DDS);
			if(solver.info()!=Eigen::Success) {
				ces << "DDS pattern analysis failed, solver.info() = "<< solver.info() << endl;
				return 1;
			}		
			solver.factorize(DDS);
			if(solver.info()!=Eigen::Success) {
				ces << "Factorization failed, solver.info() = "<< solver.info() << endl;
				return 1;
			}
			delta = solver.solve(minusDS);// use the factorization to solve for the given right hand side
			if(solver.info()!=Eigen::Success) {
				ces << "Solving failed, solver.info() = "<< solver.info() << endl;
				ces << "log(abs(det(DDS))) = " << solver.logAbsDeterminant() << endl;
				ces << "sign(det(DDS)) = " << solver.signDeterminant() << endl;
				return 1;
			}
		
			//independent check on whether calculation worked
			vec diff(2*ps.N*ps.NT+2);
			diff = DDS*delta-minusDS;
			double maxDiff = diff.maxCoeff();
			maxDiff = abs(maxDiff);
			checkInv.add(maxDiff);
			checkInv.checkMessage();
			if (!checkInv.good()) return 1;

			//assigning values to phi
			p += delta;
		
			//passing changes on to complex vector
			Cp = vecComplex(p,ps.N*ps.NT);
			
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
			if ((opts.printChoice).compare("p")==0 || (opts.printChoice).compare("e")==0) {
				Filename pEFile = basic;
				pEFile.ID = "mainpE";
				save(pEFile,so_tp,p);
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
			double maxP = p.maxCoeff();
			double minP = p.minCoeff();
			if (-minP>maxP) maxP = -minP;
			double normP = p.norm();
			double normDelta = delta.norm();
		
			// adding convergence checks
			checkAction.add(absDiff(action,action_last));
			action_last = action;
			checkSoln.add(normDS/normP);
			checkSolnMax.add(maxDS/maxP);
			checkDelta.add(normDelta/normP);
			
			//printing tests to see convergence
			if (runs_count==1) {
				fprintf(cof,"%5s%5s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s\n","loop","run","sol","solM","delta","linear"\
							,"true erg","on shell","AB","ABNE","conserv","latt");
			}
			fprintf(cof,"%5i%5i%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g\n",loop,runs_count,checkSoln.back(),\
				checkSolnMax.back(),checkDelta.back(),checkLin.back(),checkTrue.back(),checkOS.back(),checkAB.back(),checkABNE.back(),\
				checkCon.back(),checkLatt.back());
			
			if (!checkDelta.good()) {
				checkDelta.checkMessage();
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
			- p
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
		checkTrue.checkMessage();
		checkOS.checkMessage();
		checkABNE.checkMessage();
		if (ps.pot==3 && abs(ps.theta)<MIN_NUMBER) checkContm.checkMessage();
		checkCon.checkMessage();
		checkIE.checkMessage();
		checkLatt.checkMessage();
		
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
						timenumber.c_str(),fileLoop,loop,ps.N,ps.NT,ps.L,ps.dE,ps.Tb,ps.theta,angleToPrint,F,keep.c_str());
			fclose(stepOs);
			//fprintf(cof,"%12s%30s\n","steps:",stepFile.c_str());
		}
		else
			stepper.addResult(1.0); // choice irrelevant but a value is require to make step
	
		// printing results to terminal
		fprintf(cof,"\n");
		fprintf(cof,"%8s%8s%8s%8s%8s%8s%8s%8s%14s%14s%14s%14s\n","runs","time","ps.N","NT","L","Tb","dE","theta","Num","E","im(action)","W");
		fprintf(cof,"%8i%8g%8i%8i%8g%8g%8g%8g%14.4g%14.4g%14.4g%14.4g\n",\
				runs_count,realtime,ps.N,ps.NT,ps.L,ps.Tb,ps.dE,ps.theta,Num,E,imag(action),W);
		fprintf(cof,"\n");

		// printing results to file
		stepped = stepper.keep();
		if (stepped) {
			FILE * actionfile;
			string resultsFile = "results/"+timenumber+"mainResults.dat";
			actionfile = fopen(resultsFile.c_str(),"a");
			fprintf(actionfile,"%12s%5i%5i%5i%5i%6g%13.5g%13.5g%13.5g%13.5g%13.5g%13.5g%13.5g%8.2g%8.2g%8.2g\n",\
						timenumber.c_str(),fileLoop,loop,ps.N,ps.NT,ps.L,ps.Tb,ps.dE,ps.theta,E,Num,(2.0*imag(action)-bound)\
						,W,checkSoln.back(),checkLin.back(),checkTrue.back());
			fclose(actionfile);
			fprintf(cof,"%12s%30s\n","results:",resultsFile.c_str());
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
		ps.save(paramsRunFile);
	
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
		save(tpFile,so_tp,p);
		if (plotEverything)
			plot(tpFile,po_tp);
		
		if ((opts.printChoice).compare("n")!=0) {
			//printing linErg
			tpFile.ID = "mainlinErg";
			linErg.conservativeResize(ps.Na);
			save(tpFile,so_simple,linErg);
			plotFile = tpFile;
			plotFile.Suffix = ".png";
			po_simple.output = plotFile;
			if (plotEverything)
				plot(tpFile,po_simple);
	
			//printing erg
			tpFile.ID = "mainerg";
			//erg.conservativeResize(ps.Na);
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
			linErgOffShell.conservativeResize(ps.Na);
			save(tpFile,so_simple,linErgOffShell);
			
			//printing linNum
			tpFile.ID = "mainlinNum";
			tpFile.Suffix = ".dat";
			linNum.conservativeResize(ps.Na);
			save(tpFile,so_simple,linNum);
		
			//printing derivErg
			tpFile.ID = "mainderivErg";
			derivErg.conservativeResize(ps.Na);
			save(tpFile,so_simple,derivErg);
		
			//printing potErg
			tpFile.ID = "mainpotErg";
			potErg.conservativeResize(ps.Na);
			save(tpFile,so_simple,potErg);
		}
		
		if (!checkDelta.good()) {
				return 1;
			}
		fprintf(cof,"\n----------------------------------------------------------------------------------------------------------------------------\n\n");
		fclose(cof);
		} //ending parameter loop
	} //ending file loop

ces.close();

return 0;
}


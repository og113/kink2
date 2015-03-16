/*----------------------------------------------------------------------------------------------------------------------------
	pi
		program to solve boundary value problem on contour BC, and step out to A and D
----------------------------------------------------------------------------------------------------------------------------*/

//#define NDEBUG //NDEBUG is to remove error and bounds checking on vectors in SparseLU, for speed - only include once everything works
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <limits>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include "main.h"

//#define NDEBUG //NDEBUG is to remove error and bounds checking on vectors in SparseLU, for speed - only include once everything works

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - loading options
		2 - getting inputs
		3 - beginning parameter loop
		4 - assigning potential functions
		5 - omega and negVec
		6 - calculating input phi
		7 - defining quantities
		8 - beginning newton-raphson loop
		9 - assigning minusDS, DDS etc
		10 - checks
		11 - solving for delta
		12 - printing early
		13 - convergence
		14 - printing output
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
/*----------------------------------------------------------------------------------------------------------------------------
	1. loading options
----------------------------------------------------------------------------------------------------------------------------*/

Options opts;
opts.load("optionsP");

/*----------------------------------------------------------------------------------------------------------------------------
	2. getting inputs
		- loading inputs
		- defining timenumber
----------------------------------------------------------------------------------------------------------------------------*/

// loading inputs
Parameters psu;
psu.load("inputsP");
cout << "input parameters: " << endl;
psu.print();
cout << endl;

// defining timenumber
string timenumber;
(argc==2) ? timenumber = argv[1] : timenumber = currentDateTime();

/*----------------------------------------------------------------------------------------------------------------------------
	3. beginning parameter loop
		- defining a time
		- changing parameters (if required)
		- copying a verson of parameters with timenumber
		- printing parameters
		- declaring checks
----------------------------------------------------------------------------------------------------------------------------*/
	
for (uint loop=0; loop<opts.loops; loop++) {

	//defining a time and starting the clock
	clock_t time;
	time = clock();
	
	// changing parameters
	Parameters ps = psu;
	if (opts.loops>1) {
		if ((opts.loopChoice)[0]=='N') {
			uint pValue = (uint)opts.loopMin + (uint)(opts.loopMax - opts.loopMin)*loop/(opts.loops-1);
			bool anythingChanged = ps.changeParameters(opts.loopChoice,pValue);
			if (loop==0 && anythingChanged) {
				cout << opts.loopChoice << "changed to " << pValue << " on input" << endl;
			}
		}
		else {
			double pValue = opts.loopMin + (opts.loopMax - opts.loopMin)*loop/(opts.loops-1.0);
			bool anythingChanged = ps.changeParameters(opts.loopChoice,pValue);
			if (loop==0 && anythingChanged) {
				cout << opts.loopChoice << "changed to " << pValue << " on input" << endl;
			}
		}
	}
		
	// copying a version of ps with timenumber
	Filename paramsRunFile = (string)("./data/"+timenumber+"inputsP_fLoop_"+numberToString<uint>(fileLoop)\
			+"_loop_"+numberToString<uint>(loop));
	ps.save(paramsRunFile);
	
	//printing timenumber and parameters
	printf("%12s%12s\n","timenumber: ",timenumber.c_str());
	ps.print();
	
	// declaring Checks
	Check checkAction("action",1.0e-2);
	Check checkSoln("solution",1.0e-6);
	Check checkSolnMax("solution max",1.0e-5);
	Check checkDelta("delta",1.0);
	Check checkInv("matrix inversion",1.0e-16*ps.N*ps.NT);
	Check checkCon("energy conservation",1.0e-2);
	Check checkLatt("lattice small enough for energy",0.2);
	Check checkReg("regularisation term",1.0e-2);
	Check checkDT("time derivative of phi",1.0e-2);
	Check checkProf("phi input calculation",1.0e-5);
	
	Check checkLin("linear energy flat",5.0e-2);
	Check checkTrue("linear energy equal true energy",5.0e-2);
	Check checkIE("imaginary part of energy",1.0e-5);
	Check checkContm("linear energy equal continuum expression",5.0e-2);
	Check checkOS("linear energy equal on shell expression",5.0e-2);
	Check checkAB("a_k = Gamma*b_k",1.0e-2);
	Check checkABNE("N = Sum(a_k*b_k), E = Sum(w_k*a_k*b_k)",5.0e-2);
	Check checkLR("linear representation of phi",1.0e-12);
	
	// do trivial or redundant checks?
	bool trivialChecks = true;
	
/*----------------------------------------------------------------------------------------------------------------------------
	4. assigning potential functions
		- assigning potential functions
		- assigning preliminary parameter structs
----------------------------------------------------------------------------------------------------------------------------*/

	// assigning potential functions
	Potential V, dV, ddV;
	if (ps.pot==1) {
		V((PotentialType)&V1<comp>,ps);
		dV((PotentialType)&dV1<comp>,ps);
		ddV((PotentialType)&ddV1<comp>,ps);
	}
	else if (ps.pot==2) {
		V((PotentialType)&V2<comp>,ps);
		dV((PotentialType)&dV2<comp>,ps);
		ddV((PotentialType)&ddV2<comp>,ps);
	}
	else if (ps.pot==3) {
		V((PotentialType)&V3<comp>,ps);
		dV((PotentialType)&dV3<comp>,ps);
		ddV((PotentialType)&ddV3<comp>,ps);
		}
	else {
		cerr << "pot option not available, pot = " << ps.pot << endl;
		return 1;
	}
	
	// assigning preliminary parameter structs
	params_for_V paramsV  = {ps.epsilon, ps.A}, paramsV0  = {ps.epsilon0, ps.A};
	
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
	5. omega and negVec
		- loading or constructing omega
		- loading negVec
----------------------------------------------------------------------------------------------------------------------------*/

	//deterimining omega matrices for fourier transforms in spatial direction
	vec freqs(ps.N), freqs_exp(ps.N);
	mat modes(ps.N,ps.N);
	mat omega_m1(ps.N,ps.N), omega_0(ps.N,ps.N), omega_1(ps.N,ps.N), omega_2(ps.N,ps.N);
	SaveOptions so_simple;
	so_simple.paramsIn = ps; so_simple.paramsOut = ps;
	so_simple.vectorType = SaveOptions::simple;
	so_simple.extras = SaveOptions::none;
	so_simple.printMessage = false;
	{
		Filename omegaM1F, omega0F, omega1F, omega2F, modesF, freqsF, freqsExpF; // Filename works as FilenameAttributes
		omegaM1F = (string)"data/stable/omegaM1_pot_"+numberToString<uint>(ps.pot)+"_N_"+numberToString<uint>(ps.N)\
						+"_L_"+numberToString<double>(ps.L)+".dat";
		Folder omegaM1Folder(omegaM1F);
		omega0F = omegaM1F; 		omega0F.ID = "omega0"; 		Folder omega0Folder(omega0F);
		omega1F = omegaM1F; 		omega1F.ID = "omega1"; 		Folder omega1Folder(omega1F);
		omega2F = omegaM1F; 		omega2F.ID = "omega2";	 	Folder omega2Folder(omega2F);
		modesF = omegaM1F;			modesF.ID = "modes"; 		Folder modesFolder(modesF);
		freqsF = omegaM1F;			freqsF.ID = "freqs"; 		Folder freqsFolder(freqsF);
		freqsExpF = omegaM1F;		omegaM1F.ID = "freqsExp";	Folder freqsExpFolder(omegaM1F);
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
			bool approxOmega = false;
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
			Filename eigVecFile = "data/stable/eigVec_pot_3_L_" + numberToString<double>(ps.L) + ".dat";	
			load(eigVecFile,so_simple,negVec); // should automatically interpolate
		}
		else {
			Filename eigVecFile;
			string N_load;
			string Nb_load;
			Filename lower = "data/stable/eigVec_pot_" + numberToString<uint>(ps.pot)\
								 + "_N_100_Nb_100_L_" + numberToString<double>(ps.L) + ".dat";
			Filename upper = lower;
			vector<StringPair> upperExtras = lower.Extras;
			upper.Extras[1] = StringPair("N","1000");
			upper.Extras[2] = StringPair("Nb","1000");
			Folder eigVecFolder(lower,upper);
			if (eigVecFolder.size()==0) {
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
						NbcheckProf_load = ((eigVecFile.Extras)[2]).second;
					}
				}
			}
			SaveOptions eigVecOpts;
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
	6. defining quantities
		- erg, linErg etc
		- p, minusDS, DDS
----------------------------------------------------------------------------------------------------------------------------*/

	// erg, linErg etc
	comp action = ps.action0;
	double W;
	double E;
	cVec erg(NT);	

	//defining some quantities used to stop the Newton-Raphson loop when action stops varying
	comp action_last = action;
	uint runs_count = 0;
	uint min_runs = 3;

	//initializing phi (=p), DDS and minusDS
	vec p(2*ps.N*ps.Nb+1);
	p = Eigen::VectorXd::Zero(2*ps.N*ps.Nb+1);
	spMat DDS(2*N*Nb+1,2*N*Nb+1);
	vec minusDS(2*N*Nb+1);

/*----------------------------------------------------------------------------------------------------------------------------
	7. calculating input phi
		- loading phi
		- finding phi profile between minima
		- assigning phi
		- fixing boundary conditions
		- printing input phi
		- Cp
----------------------------------------------------------------------------------------------------------------------------*/
	
	// loading phi (if required)
	SaveOptions so_p;
	so_p.paramsIn = ps;
	so_p.paramsOut = ps;
	so_p.vectorType = SaveOptions::complexB;
	so_p.extras = SaveOptions::coords;
	so_p.zeroModes = 1;
	so_p.printMessage = false;
	if ((opts.inF).compare("f")==0) {
		Filename inputsPhiFile = (string)("data/"+opts.minTimenumberLoad+"inputsP_loop_"+\
							minLoopLoad+".dat");
		Parameters paramsPhiIn;
		paramsPhiIn.load(inputsPhiFile);
		so_p.paramsIn = paramsPhiIn;
		Filename phiFile = inputsPhiFile;
		phiFile.ID = "p";
		load(phiFile,so_p,p);
		so_p.paramsIn = ps;
	}
	else if (loop>0) {
		Filename phiFile = (string)("data/"+opts.minTimenumberLoad+"p_loop_"+\
							numberToString<uint>(loop-1)+".dat");
		load(phiFile,so_p,p);
	}
	else {
		//finding phi profile between minima
		uint profileSize = ps.Nb; //more than the minimum
		vector<double> phiProfile(profileSize);
		vector<double> rhoProfile(profileSize);
		double alphaL = ps.alpha, alphaR = ps.alpha;

		if (ps.pot!=3) {
			if (ps.pot==2) {
				double phiL = ps.minima0[1]-1.0e-2;
				double phiR = ps.minima0[0]+1.0e-2;
				for (uint j=0;j<profileSize;j++) {
					phiProfile[j] = phiL + (phiR-phiL)*j/(profileSize-1.0);
				}
	
				double profileError;
				gsl_function rho_integrand;
				rho_integrand.function = &rhoIntegrand;
				rho_integrand.params = &paramsV0;
				gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
				w = gsl_integration_workspace_alloc(1e4);
				for (uint j=0;j<profileSize;j++) {
					gsl_integration_qags(&rho_integrand, phiProfile[j], 0, 1.0e-16, 1.0e-6, 1e4, w,\
															 &(rhoProfile[j]), &profileError);
					checkProfile.add(profileError);
					checkProfile.checkMessage();
				}
				gsl_integration_workspace_free(w);
				alphaL = rhoProfile[0];
				alphaR = rhoProfile.back();
			}
		}
	}
	
	

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//assigning input phi
	if (inP.compare("f")!=0 && loop==0)
		{
		if (ps.pot==3)
			{
			vec tempPhi = loadVectorColumn("data/instanton00.dat",3);
			vec temp2Phi;
			unsigned int length = tempPhi.size();
			length = (unsigned int)(sqrt(length));
			temp2Phi = interpolateReal(tempPhi,length,length,Nb,N);
			for (unsigned int j=0; j<N*Nb; j++)
				{
				p[2*j] = temp2Phi[j];
				p[2*j+1] = 0.0;
				}
			p[2*Nb*N] = 0.5; //for the zero mode
			}
		else
			{
			if (R<alpha)
				{
				cout << "R is too small. Not possible to give thinwall input. It should be more that " << alpha;
				}
			//#pragma omp parallel for
			for (unsigned int j=0; j<N*Nb; j++)
				{
				comp t = coordB(j,0);
				comp x = coordB(j,1);
				p(2*j+1) = 0.0; //imaginary parts set to zero
				if (inP.compare("b")==0 || (inP.compare("p")==0 && Tb>R))
					{
					double rho = real(sqrt(-pow(t,2.0) + pow(x,2.0))); //should be real even without real()
					if (ps.pot==1)
						{
						if ((rho-R)<-alpha) 		p(2*j) = minima[1];
						else if ((rho-R)>alpha) 	p(2*j) = minima[0];
						else						p(2*j) = (minima[1]+minima[0])/2.0\
														 + (minima[0]-minima[1])*tanh((rho-R)/2.0)/2.0;
						}
					else if (ps.pot==2)
						{
						if ((rho-R)<=alphaL) 		p(2*j) = minima[1];
						else if ((rho-R)>=alphaR) 	p(2*j) = minima[0];
						else
							{
							vector<double> rhoPos (profileSize,rho-R);
							for (unsigned int k=0; k<profileSize; k++)
								{
								rhoPos[k] -= rhoProfile[k];
								}
							unsigned int minLoc = smallestLoc(rhoPos);
			                p(2*j) = phiProfile[minLoc];
							}
						}
					if (inP.compare("p")==0)
						{
						p(2*j) += amp*cos(sqrt(-negVal)*imag(t))*negVec(2*j);
						p(2*j+1) += amp*cos(sqrt(-negVal)*imag(t))*negVec(2*j+1);
						}
					}
				else if (inP.compare("p")==0 && Tb<R)
					{
					double rho1 = real(sqrt(-pow(t,2.0) + pow(x+R*cos(angle),2.0)));
					double rho2 = real(sqrt(-pow(t,2.0) + pow(x-R*cos(angle),2.0)));
					if ((rho1-R)<-alpha && (rho2-R)<-alpha)		p(2*j) = minima[1];
					else if ((rho1-R)>alpha || (rho2-R)>alpha)	p(2*j) = minima[0];
					else if (real(x)>0)							p(2*j) = (minima[1]+minima[0])/2.0\
																 	+ (minima[0]-minima[1])*tanh((rho1-R)/2.0)/2.0;
					else if (real(x)<0)							p(2*j) = (minima[1]+minima[0])/2.0 \
																	+ (minima[0]-minima[1])*tanh((rho2-R)/2.0)/2.0;
					else										p(2*j) = minima[1]; //i.e. if coordB(j,1) == 0
					}
				}
			}
		p(2*N*Nb) = 0.5; //initializing Lagrange parameter for removing dp/dx zero mode
		}
	else
		{
		string loadfile, inputsFile;
		if (inP.compare("f")==0)
			{
			loadfile = "./data/" + aq.inputTimeNumber + "pip_" + aq.inputLoop + ".dat";
			cout << "input: " << loadfile << endl;
			inP = "p";
			inputsFile = "./data/" + aq.inputTimeNumber + "inputsPi_" + aq.inputLoop;
			}
		else
			{
			loadfile = "./data/" + timeNumber + "pi"+inP+"_" + numberToString<int>(loop-1)+".dat";
			inputsFile = "./data/" + timeNumber + "inputsPi_" + numberToString<int>(loop-1);
			}
		unsigned int fileLength = countLines(loadfile);
		if (fileLength==(N*Nb+1))
			{
			p = loadVector(loadfile,Nb,N,1);
			}
		else if (fileLength % 2) //if its odd
			{
			unsigned int Nin, Ntin;
			//cout << "interpolating input, filelength = " << fileLength << " , Cp.size() = " << N*Nb+1 << endl;
			ifstream fin;
			fin.open(inputsFile.c_str());
			if (fin.is_open())
				{
				string line, temp;
				while(getline(fin,line))
					{
					if(line[0] == '#') continue;
					istringstream ss(line);
					ss >> Nin >> temp >> Ntin;
					break;
					}
				}
			else{cout << "unable to open " << inputsFile << endl;}
			fin.close();
			vec temp_p = loadVector(loadfile,Ntin,Nin,1);
			p = interpolate(temp_p,Ntin,Nin,Nb,N);
			}
		else
			{
			cout << "vector in file: " << loadfile << " has length " << fileLength << " so cannot load into p" << endl;
			}
		}
	
	//fixing input periodic instanton to have zero time derivative at time boundaries
	//#pragma omp parallel for
	for (unsigned int j=0;j<N;j++)
		{
	    p(2*j*Nb) = (1.0-open)*p(2*j*Nb) + open*p(2*(j*Nb+1)); //initial time real
	    p(2*(j*Nb+1)) = p(2*j*Nb);
	    p(2*j*Nb+1) = (1.0-open)*p(2*j*Nb+1) + open*p(2*(j*Nb+1)+1); //initial time imag
	    p(2*(j*Nb+1)+1) = p(2*j*Nb+1);
	    p(2*((j+1)*Nb-1)) = open*p(2*((j+1)*Nb-2)) + (1.0-open)*p(2*((j+1)*Nb-1)); //final time real
	    p(2*((j+1)*Nb-2)) = p(2*((j+1)*Nb-1));
	    p(2*((j+1)*Nb-2)+1) = open*p(2*((j+1)*Nb-1)+1) + (1.0-open)*p(2*((j+1)*Nb-2)+1); //final time imag
	    p(2*((j+1)*Nb-1)+1) = p(2*((j+1)*Nb-2)+1);
		}
    if (ps.pot==3)
    	{
		for (unsigned int j=0;j<Nb;j++)
			{
			unsigned int l0 = c(j,0,Nb);
			unsigned int m = c(j,N-1,Nb);
			p(2*l0) = 0.0; // p=0 ar r=0
			p(2*l0+1) = 0.0;
		    p(2*m) = 0.0; //p=0 ar r=R
		    p(2*m+1) = 0.0;
			}
		}
		
	//very early vector print
	string earlyPrintFile = "data/" + timeNumber + "piE"+inP+ "_" + numberToString<unsigned int>(loop) + "_0.dat";
	printVectorB(earlyPrintFile,p);
		
	//defining complexified vector Cp
	cVec Cp(Nb*N);
	Cp = vecComplex(p,N*Nb);
	
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//beginning newton-raphson loop
	while ((sol_test.back()>closenessS || solM_test.back()>closenessSM || runs_count<min_runs))
		{
		runs_count ++;
		
		//defining the zero mode at the final time boundary and the time step before
		vec Chi0(Nb*N);
		Chi0 = Eigen::VectorXd::Zero(N*Nb);
		//#pragma omp parallel for
		for (unsigned int j=0; j<N; j++)
			{
			unsigned int pos = (j+1)*Nb-1;
			long int neighPos = neigh(pos,1,1,Nb,N), neighMin = neigh(pos,1,-1,Nb,N);
			if(neighPos!=-1 && neighMin!=-1)
            	{
            	Chi0(pos) = p(2*neighPos)-p(2*neighMin); //final time slice
            	//Chi0(pos-1) = p(2*neigh(pos-1,1,1,Nb,N))-p(2*neigh(pos-1,1,-1,Nb,N)); //penultimate time slice
            	}
            }
        if (runs_count==1) //printing Chi0
        	{
        	//printVectorB("data/" + timeNumber + "Chi0.dat",Chi0);
        	}

		// allocating memory for DS, DDS
		minusDS = Eigen::VectorXd::Zero(2*N*Nb+1); //initializing to zero
		DDS.setZero(); //just making sure
		Eigen::VectorXi DDS_to_reserve(2*N*Nb+1);//number of non-zero elements per column
		DDS_to_reserve = Eigen::VectorXi::Constant(2*N*Nb+1,11);
		DDS_to_reserve(0) = 3; //these need to be changed when boundary conditions need to be more compicated
		DDS_to_reserve(1) = 3;
		DDS_to_reserve(2*N*Nb-2) = 3;
		DDS_to_reserve(2*N*Nb-1) = 3;
		DDS_to_reserve(2*N*Nb) = N;
		DDS.reserve(DDS_to_reserve);
		
		//initializing to zero
		comp kineticS = 0.0;
		comp kineticT = 0.0;
		comp pot_0 = 0.0;
		comp pot_r = 0.0;
		erg = Eigen::VectorXcd::Constant(NT,-ergZero);
		double dtTest = 0.0;
		
		//testing that the potential term is working for pot3
		if (ps.pot==3 && false)
			{
			comp Vtrial = 0.0, Vcontrol = 0.0;
			for (unsigned int j=0; j<N; j++)
				{
				double r = r0 + j*a;
				paramsV  = {r, 0.0};
				Vcontrol += pow(p(2*j*Nb),2.0)/2.0 - pow(p(2*j*Nb),4.0)/4.0/pow(r,2.0);
				Vtrial += V(p(2*j*Nb));
				}
			double potTest = pow(pow(real(Vcontrol-Vtrial),2.0) + pow(imag(Vcontrol-Vtrial),2.0),0.5);
			cout << "potTest = " << potTest << endl;
			}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//assigning values to minusDS and DDS and evaluating action
		//#pragma omp parallel for
		for (unsigned long int j = 0; j < N*Nb; j++)
			{		
			unsigned int t = intCoord(j,0,Nb); //coordinates
			unsigned int x = intCoord(j,1,Nb);
			long int neighPosX = neigh(j,1,1,Nb,N);
			if (ps.pot==3)
				{
				paramsV  = {r0+x*a, 0.0};
				}
			if (t<(Nb-1)) dtTest += abs((p(2*(j+1))-p(2*j))/b);
			
			if ((absolute(Chi0(j))>1.0e-16) && ps.pot!=3) //zero mode lagrange constraint
				{
				DDS.insert(2*j,2*N*Nb) = a*Chi0(j); 
				DDS.insert(2*N*Nb,2*j) = a*Chi0(j);
				minusDS(2*j) += -a*Chi0(j)*p(2*N*Nb);
				minusDS(2*N*Nb) += -a*Chi0(j)*p(2*j);
		    	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//boundaries
			if (ps.pot==3 && x==(N-1))
				{
				DDS.insert(2*j,2*j) = 1.0; // p=0 at r=R
				DDS.insert(2*j+1,2*j+1) = 1.0;
				}
			else if (ps.pot==3 && x==0)
				{
				DDS.insert(2*j,2*j) = 1.0; // p=0 at r=0
				DDS.insert(2*j+1,2*j+1) = 1.0;
				double csi = ((t==0 || t==(NT-1))? 0.5: 1.0);
				kineticS +=	csi*b*pow(Cp(neighPosX),2.0)/a/2.0;
				}
			else if (t==(Nb-1))
				{
				comp Dt = -b*ii/2.0;
				erg(t+Na) += pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				kineticS += Dt*pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0; //n.b. no contribution from time derivative term at the final time boundary	
				pot_0 += Dt*a*V(Cp(j));
				pot_r += Dt*a*Vr(Cp(j));
				
				
				DDS.insert(2*j,2*j) = 1.0/b; //zero time derivative
				DDS.insert(2*j,2*(j-1)) = -1.0/b;
				DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				}
			else if (t==0)
				{
				comp dt = -b*ii;
				comp Dt = -b*ii/2.0;
				kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
				kineticS += Dt*pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0;
				pot_0 += Dt*a*V(Cp(j));
				pot_r += Dt*a*Vr(Cp(j));
				erg(t+Na) += a*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				
				if (inP.compare("b")==0)
					{
					DDS.insert(2*j,2*j) = 1.0; //zero change (initial input satisfies b.c.s)
					DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
					}
				else if (inP.compare("p")==0)
					{
					DDS.insert(2*j,2*j) = -1.0/b; //zero time derivative
					DDS.insert(2*j,2*(j+1)) = 1.0/b;
					DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
					}
				}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//bulk
			else
				{
				comp dt = -b*ii;
				comp Dt = -b*ii;
				kineticT += a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
				kineticS += Dt*pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0;
				pot_0 += Dt*a*V(Cp(j));
				pot_r += Dt*a*Vr(Cp(j));
				erg(t+Na) += a*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neighPosX)-Cp(j),2.0)/a/2.0 + a*V(Cp(j)) + a*Vr(Cp(j));
				
                for (unsigned int k=0; k<2*2; k++)
                	{
                    int sign = pow(-1,k);
                    int direc = (int)(k/2.0);
                    long int neighb = neigh(j,direc,sign,Nb,N);
                    if (direc == 0)
                    	{
                        minusDS(2*j) += real(a*Cp(j+sign)/dt);
                        minusDS(2*j+1) += imag(a*Cp(j+sign)/dt);
                        DDS.insert(2*j,2*(j+sign)) = -real(a/dt);
                        DDS.insert(2*j,2*(j+sign)+1) = imag(a/dt);
                        DDS.insert(2*j+1,2*(j+sign)) = -imag(a/dt);
                        DDS.insert(2*j+1,2*(j+sign)+1) = -real(a/dt);
                        }
                    else if (neighb !=-1)
                    	{
                        minusDS(2*j) += - real(Dt*Cp(neighb)/a);
                        minusDS(2*j+1) += - imag(Dt*Cp(neighb)/a);
                        DDS.insert(2*j,2*neighb) = real(Dt/a);
                        DDS.insert(2*j,2*neighb+1) = -imag(Dt/a);
                        DDS.insert(2*j+1,2*neighb) = imag(Dt/a);
                        DDS.insert(2*j+1,2*neighb+1) = real(Dt/a);
                        }
                    }
                comp temp0 = 2.0*a/dt;
	            comp temp1 = a*Dt*(2.0*Cp(j)/pow(a,2.0) + dV(Cp(j)) + dVr(Cp(j)));
	            comp temp2 = a*Dt*(2.0/pow(a,2.0) + ddV(Cp(j)) + ddVr(Cp(j)));
	                
	            minusDS(2*j) += real(temp1 - temp0*Cp(j));
	            minusDS(2*j+1) += imag(temp1 - temp0*Cp(j));
	            DDS.insert(2*j,2*j) = real(-temp2 + temp0);
	            DDS.insert(2*j,2*j+1) = imag(temp2 - temp0);
	            DDS.insert(2*j+1,2*j) = imag(-temp2 + temp0);
	            DDS.insert(2*j+1,2*j+1) = real(-temp2 + temp0);
	            }
            }
        action = kineticT - kineticS - pot_0 - pot_r;
        if (ps.pot==3) {
        	DDS.insert(2*N*Nb,2*N*Nb) = 1.0;
        	dtTest /= (double)(Nb-1.0);
        	action *= 4.0*pi;
        }
        
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			string prefix = "./data/" + timeNumber;
			string suffix = inP+"_" + numberToString<unsigned int>(loop)+"_" + numberToString<unsigned int>(runs_count)+".dat";
			if ((print_choice.compare("v")==0 || print_choice.compare("e")==0))
				{
				string minusDSfile = prefix + "minusDSE"+suffix;
				printVectorB(minusDSfile,minusDS);
				}
			if ((print_choice.compare("p")==0 || print_choice.compare("e")==0))
				{
				string piEarlyFile = prefix + "piE"+suffix;
				printVectorB(piEarlyFile,p);
				}
			if ((print_choice.compare("m")==0 || print_choice.compare("e")==0))
				{
				string DDSfile = prefix + "DDSE"+suffix;
				printSpmat(DDSfile,DDS);
				}
			}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

	//solving for delta in DDS*delta=minusDS, where p' = p + delta		
		vec delta(2*N*Nb+1);
		delta = Eigen::VectorXd::Zero(2*N*Nb+1);
		DDS.makeCompressed();
		Eigen::SparseLU<spMat> solver;
		
		solver.analyzePattern(DDS);
		if(solver.info()!=Eigen::Success)
			{
			cerr << "DDS pattern analysis failed, solver.info() = "<< solver.info() << endl;
			return 1;
			}		
		solver.factorize(DDS);
		if(solver.info()!=Eigen::Success) 
			{
			cerr << "Factorization failed, solver.info() = "<< solver.info() << endl;
			return 1;
			}
		delta = solver.solve(minusDS);// use the factorization to solve for the given right hand side
		if(solver.info()!=Eigen::Success)
			{
			cerr << "Solving failed, solver.info() = "<< solver.info() << endl;
			cerr << "log(abs(det(DDS))) = " << solver.logAbsDeterminant() << endl;
			cerr << "sign(det(DDS)) = " << solver.signDeterminant() << endl;
			return 1;
			}
		
		//independent check on whether calculation worked
		vec diff(2*N*Nb+1);
		diff = DDS*delta-minusDS;
		double maxDiff = diff.maxCoeff();
		maxDiff = absolute(maxDiff);
		calc_test.push_back(maxDiff);
		if (calc_test.back()>closenessC)
			{
			cout << "Calculation failed" << endl;
			cout << "calc_test = " << calc_test.back() << endl;
			return 1;
			}

		//assigning values to phi
		p += delta;
		
		//passing changes on to complex vector
		Cp = vecComplex(p,N*Nb);
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		//convergence issues
		
		//checking pot_r is much smaller than the other potential terms
		reg_test.push_back(abs(pot_r/pot_0));
		if (reg_test.back()>closenessR) cerr << "regularisation term is too large, regTest = " << reg_test.back() << endl;
		
		//evaluating norms
		double normDS = minusDS.dot(minusDS);
		normDS = pow(normDS,0.5);
		double maxDS = minusDS.maxCoeff();
		double minDS = minusDS.minCoeff();
		if (-minDS>maxDS)
			{
			maxDS = -minDS;
			}
		double normP = p.dot(p);
		normP = pow(normP,0.5);
		double normDelta = delta.dot(delta);
		normDelta = pow(normDelta,0.5);
		
		//assigning test values
		//quantities used to stop newton-raphson loop
		action_test.push_back(abs(action - action_last)/abs(action_last));
		action_last = action;
		sol_test.push_back(normDS/normP);
		solM_test.push_back(maxDS);
		delta_test.push_back(normDelta/normP);
		dt_test.push_back(dtTest);
			
		//printing tests to see convergence
		if (runs_count==1)
			{
			printf("%16s%16s%16s%16s%16s%16s\n","loop","runsCount","actionTest","solTest","solMTest","deltaTest");
			}
		printf("%16i%16i%16g%16g%16g%16g\n",loop,runs_count,action_test.back(),sol_test.back(),solM_test.back(),delta_test.back());
		
		} //closing "runs" while loop
		
		if (abs(dt_test.back())<closenessDT) {
		cout << endl << "average |dp/dt| integrated over each timeslice:" << endl;
		for (unsigned int j=0; j<dt_test.size();j++) cout << dt_test[j] << endl;
		cout << endl;
		}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	//propagating solution along minkowskian time paths
	
	if (ps.pot!=3) {
		//propagating solution back in minkowskian time
		//A1. initialize mp==mphi using last point of ephi and zeros- use complex phi
		cVec ap(N*(Na+1)); //phi on section "a"
		ap = Eigen::VectorXcd::Zero(N*(Na+1));
		//#pragma omp parallel for
		for (unsigned int j=0; j<N; j++) ap(j*(Na+1)) = Cp(j*Nb);

		//A2. initialize vel - defined at half steps, first step being at t=-1/2,
		//vel(t+1/2) := (p(t+1)-p(t))/dt
		cVec velA (N*(Na+1));
		velA = Eigen::VectorXcd::Zero(N*(Na+1));
		double dtau = -b;
		double Dt0 = dtau; //b/2*(-1+1i); - this is surely wrong!!
		//#pragma omp parallel for
		//for (unsigned int j=0; j<N; j++)
			//{
		    //velA(j*(Na+1)) = 0; //due to boundary condition
			//}

		
		//A3. initialize acc using phi and expression from equation of motion and zeros-complex
		cVec accA(N*(Na+1));
		accA = Eigen::VectorXcd::Zero(N*(Na+1));
		//#pragma omp parallel for
		for (unsigned int j=0; j<N; j++)
			{
			if (ps.pot==3) 			paramsV  = {r0+j*a, A};
			unsigned int l = j*(Na+1);
			if (ps.pot==3 && j==(N-1)) 	accA(l) = 0.0;
			else if (ps.pot==3 && j==0) 	accA(l) = ((Dt0/pow(a,2.0))*(2.0*ap(neigh(l,1,1,Na+1,N))-2.0*ap(l)) \
		        								-Dt0*(dV(ap(l))+dVr(ap(l))))/dtau;
			else 							accA(l) = ((Dt0/pow(a,2.0))*(ap(neigh(l,1,1,Na+1,N))\
												+ap(neigh(l,1,-1,Na+1,N))-2.0*ap(l))-Dt0*(dV(ap(l))+dVr(ap(l))))/dtau;
			}
			
		//A4.5 starting the energy and that off
		vec linErgA(Na); linErgA = Eigen::VectorXd::Zero(Na);
		vec linNumA(Na); linNumA = Eigen::VectorXd::Zero(Na);

		//A7. run loop
		for (unsigned int t=1; t<(Na+1); t++)
			{
			//#pragma omp parallel for
		    for (unsigned int x=0; x<N; x++)
		    	{
		        unsigned int m = t+x*(Na+1);
		        velA(m) = velA(m-1) + dtau*accA(m-1);
		        ap(m) = ap(m-1) + dtau*velA(m);
		    	}
		    //#pragma omp parallel for	
		    for (unsigned int x=0; x<N; x++)
		    	{
		    	if (ps.pot==3) paramsV  = {r0+x*a, A};
		        unsigned int m = t+x*(Na+1);
		        if (ps.pot==3 && x==(N-1)) accA(m) = 0.0;
				else if (ps.pot==3 && x==0)
					{
					accA(m) = (1.0/pow(a,2.0))*(2.0*ap(neigh(m,1,1,Na+1,N))-2.0*ap(m)) \
		        		-dV(ap(m)) - dVr(ap(m));
		        	erg(Na-t) += a*pow(ap(m-1),2.0)/pow((-dtau),2.0)/2.0;
					}
				else
					{
					accA(m) = (1.0/pow(a,2.0))*(ap(neigh(m,1,1,Na+1,N))+ap(neigh(m,1,-1,Na+1,N))-2.0*ap(m)) \
		        		-dV(ap(m)) - dVr(ap(m));
		        	erg(Na-t) += a*pow(ap(m-1)-ap(m),2.0)/pow(-dtau,2.0)/2.0 + pow(ap(neigh(m,1,1,Na+1,N))-ap(m),2.0)/a/2.0 \
		        		+ a*V(ap(m)) + a*Vr(ap(m));
				    }
		        for (unsigned int y=0; y<N; y++)
		        	{
		        	unsigned int n = t + y*(Na+1);
				    if (absolute(theta)<1.0e-16)
						{
				    	linErgA(Na-t) += Eomega(x,y)*(real(ap(m))-minima[0])*(real(ap(n))-minima[0]) + Eomega(x,y)*imag(ap(m))*imag(ap(n));
						linNumA (Na-t) += omega(x,y)*(real(ap(m))-minima[0])*(real(ap(n))-minima[0]) + omega(x,y)*imag(ap(m))*imag(ap(n));
				    	}
					else
						{
						linErgA(Na-t) += 2.0*Gamma*omega(x,y)*(real(ap(m))-minima[0])*(real(ap(n))-minima[0])/pow(1.0+Gamma,2.0) + 2.0*Gamma*omega(x,y)*imag(ap(m))*imag(ap(n))/pow(1.0-Gamma,2.0);
						linNumA(Na-t) += 2.0*Gamma*Eomega(x,y)*(real(ap(m))-minima[0])*(real(ap(n))-minima[0])/pow(1.0+Gamma,2.0) + 2.0*Gamma*Eomega(x,y)*imag(ap(m))*imag(ap(n))/pow(1.0-Gamma,2.0);
						}
					}
		    	}
			}
		
		//E = 0;
		//unsigned int linearInt = (unsigned int)(Na/6);
		//#pragma omp parallel for
		//for (unsigned int j=0; j<linearInt; j++)
		//	{
		//	E += real(erg(j));
		//	}
		//E /= linearInt;
		if (ps.pot==3) {
			linErgA *= 4.0*pi;
			linNumA *= 4.0*pi;
			erg *= 4.0*pi;
		}
		E = linErgA(0);
		W = - E*2.0*Tb + 2.0*imag(action);

		//now propagating forwards along c
		//C2. initialize mp==mphi using last point of ephi and zeros- use complex phi
		cVec ccp(N*(Nc+1)); //phi on section "c"
		ccp = Eigen::VectorXcd::Zero(N*(Nc+1));
		//#pragma omp parallel for
		for (unsigned int j=0; j<N; j++) ccp(j*(Nc+1)) = Cp(j*Nb+Nb-1);

		//C3. initialize vel - defined at half steps, first step being at t=-1/2,
		//vel(t+1/2) := (p(t+1)-p(t))/dt
		cVec velC (N*(Nc+1));
		velC = Eigen::VectorXcd::Zero(N*(Nc+1));
		dtau = b;
		Dt0 = dtau; //b/2*(1-1i); - this is surely wrong!!
		//#pragma omp parallel for
		//for (unsigned int j=0; j<N; j++)
			//{
		    //velC(j*(Nc+1)) = 0; //due to boundary condition
			//}

		//C4. initialize acc using phi and expression from equation of motion and zeros-complex
		cVec accC(N*(Nc+1));
		accC = Eigen::VectorXcd::Zero(N*(Nc+1));
		//#pragma omp parallel for
		for (unsigned int j=0; j<N; j++)
			{
			unsigned int l = j*(Nc+1);
			if (ps.pot==3) paramsV  = {r0+j*a, A};
			if (ps.pot==3 && j==(N-1)) 	accC(l) = 0.0;
			else if (ps.pot==3 && j==0) 	accC(l) = ((Dt0/pow(a,2.0))*(2.0*ccp(neigh(l,1,1,Nc+1,N))-2.0*ccp(l)) \
		        							-Dt0*(dV(ccp(l))+dVr(ccp(l))))/dtau;
			else 							accC(l) = ((Dt0/pow(a,2.0))*(ccp(neigh(l,1,1,Nc+1,N))\
												+ccp(neigh(l,1,-1,Nc+1,N))-2.0*ccp(l))-Dt0*(dV(ccp(l))+dVr(ccp(l))))/dtau;
			}

		//C7. run loop
		for (unsigned int t=1; t<(Nc+1); t++)
			{
			//#pragma omp parallel for
			for (unsigned int x=0; x<N; x++)
				{
				unsigned int l = t+x*(Nc+1);
				velC(l) = velC(l-1) + dtau*accC(l-1);
				ccp(l) = ccp(l-1) + dtau*velC(l);
				}
			//#pragma omp parallel for
			for (unsigned int x=0; x<N; x++)
				{
				if (ps.pot==3) paramsV  = {r0+x*a, A};
				unsigned int l = t+x*(Nc+1);
				if (ps.pot==3 && x==(N-1)) 	accC(l) = 0.0;
				else if (ps.pot==3 && x==0) 	accC(l) = (1.0/pow(a,2.0))*(2.0*ccp(neigh(l,1,1,Nc+1,N))-2.0*ccp(l)) \
													-dV(ccp(l));
				else 							accC(l) = (1.0/pow(a,2.0))*(ccp(neigh(l,1,1,Nc+1,N)) \
													+ccp(neigh(l,1,-1,Nc+1,N))-2.0*ccp(l))-dV(ccp(l));		    	
				if (ps.pot==3 && t>1 && x!=(N-1) && x!=0)
					{
					erg (Na+Nb-2+t) += a*pow(ccp(l)-ccp(l-1),2.0)/pow(dtau,2.0)/2.0\
					 	+ pow(ccp(neigh(l-1,1,1,Nc+1,N))-ccp(l-1),2.0)/a/2.0\
						+ a*V(ccp(l-1)) + a*Vr(ccp(l-1));
					}
				else if (ps.pot==3 && t>1 && x==0) erg (Na+Nb-2+t) += pow(ccp(neigh(l-1,1,1,Nc+1,N)),2.0)/a/2.0;
				}
			}
		
		//checking energy conserved
		double ergChange = 0.0;
		double relErgChange = 0.0;
		if (absolute(real(erg(0)))>1.0e-16)
			{
			ergChange = absolute(real(erg(0))-real(erg(NT-2)));
			relErgChange = absolute((real(erg(0))-real(erg(NT-2)))/real(erg(0)));
			}
		erg_test.push_back(ergChange);
		if (erg_test.back()>closenessCon)
			{
			cout << "energy change = " << ergChange << endl;
			cout << "relative energy change = " << relErgChange << endl;
			}
		

		//12. combine phi with ap and cp and save combination to file
		cVec tCp(NT*N);
		//#pragma omp parallel for
		for (unsigned int j=0; j<NT*N; j++)
			{
		    unsigned int t = intCoord(j,0,NT);
		    unsigned int x = intCoord(j,1,NT);
		    if (t<Na)
		    	{
		        t = Na-t;
		        tCp(j) = ap(t+x*(Na+1));
		        }
		    else if (t<(Na+Nb))
		    	{
		        t = t - Na;
		        tCp(j) = Cp(t+x*Nb);
		        }
		    else
		    	{
		        t = t - Na - Nb + 1;
		        tCp(j) = ccp(t+x*(Nc+1));
		    	}
			}
			
		//making real vec from complex one
		vec tp(2*N*NT);
		tp = vecReal(tCp,NT*N);
		tp.conservativeResize(2*N*NT+1);
		tp(2*N*NT) = p(2*N*Nb);
		
		string prefix = "./data/" + timeNumber;
		string suffix = "_" + numberToString<unsigned int>(loop) + ".dat";
		
		//printing output phi on whole time contour
		string tpifile = prefix + "tpi"+inP+suffix;
		printVector(tpifile,tp);
		printf("%12s%30s\n"," ",tpifile.c_str());
		//gp(tpifile,"repi.gp");
		
		//printing linErgVec
		string linErgFile = "./data/" + timeNumber + "linErg"+inP+suffix;
		simplePrintVector(linErgFile,linErgA);
		printf("%12s%30s\n"," ",linErgFile.c_str());
	//	gpSimple(linErgFile);

    	}
    else {
    	E = 0.0;
    	W = 0.0;
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	    	//misc end of program tasks - mostly printing
    
    //stopping clock
	time = clock() - time;
	double realtime = time/1000000.0;
	
	//printing to terminal
	printf("\n");
	printf("%8s%8s%8s%8s%8s%8s%8s%12s%12s%12s\n","runs","time","N","NT","L","Tb","dE","E","im(S)","W");
	printf("%8i%8g%8i%8i%8g%8g%8g%12.4g%12.4g%12.4g\n",runs_count,realtime,N,NT,L,Tb,dE,E,imag(action),W);
	printf("\n");
	 printf("%60s\n","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

	//printing action value
	FILE * actionfile;
	actionfile = fopen("./data/action.dat","a");
	fprintf(actionfile,"%16s%8i%8i%8g%8g%8g%10.4g%10.4g%10.4g%10.4g%10.4g\n",timeNumber.c_str(),N,NT,L,Tb,dE,E,imag(action)\
	,W,sol_test.back(),erg_test.back());
	fclose(actionfile);
	
	string prefix = "./data/" + timeNumber;
	string suffix = "_" + numberToString<unsigned int>(loop) + ".dat";
	
	//copying a version of inputs with timeNumber
	string runInputs = prefix + "inputsPi" + "_" + numberToString<unsigned int>(loop);
	if (loop_choice[0] == 'N') 				changeInputs(runInputs,loop_choice, numberToString<unsigned int>(intLoopParameter));
	else if (loop_choice.compare("n")!=0) 	changeInputs(runInputs,loop_choice, numberToString<double>(doubleLoopParameter));
	else									copyFile("inputs",runInputs);
	printf("%12s%30s\n","output: ",runInputs.c_str());

	//printing output phi on Euclidean time part
	string pifile = prefix + "pi"+inP+suffix;
	printVectorB(pifile,p);
	printf("%12s%30s\n"," ",pifile.c_str());
//	gp(pifile,"repi.gp");
	
	//printing output minusDS				
	string minusDSfile = "./data/" + timeNumber + "minusDS"+suffix;
	printVectorB(minusDSfile,minusDS);
	printf("%12s%30s\n"," ",minusDSfile.c_str());
				
	//printing output DDS
	string DDSfile = prefix + "DDS"+inP+suffix;
	printSpmat(DDSfile,DDS);
	printf("%12s%30s\n"," ",DDSfile.c_str());
	
	//printing erg
	string ergFile = prefix + "erg"+inP+suffix;
	simplePrintCVector(ergFile,erg);
	printf("%12s%30s\n"," ",ergFile.c_str());
//	gpSimple(ergFile);
	
	//printing error, and eigenvalue to file
	if (ps.pot!=3)
		{
		ifstream eigenvalueIn;
		string eigenvaluefile = "data/eigValue.dat";
		eigenvalueIn.open(eigenvaluefile.c_str(),ios::in);
		string lastEigLine = getLastLine(eigenvalueIn);
		eigenvalueIn.close();
		ofstream eigenvalueOut;
		string eigValueFile = prefix + "eigValue.dat";
		eigenvalueOut.open(eigValueFile.c_str(),ios::out);
		eigenvalueOut << lastEigLine;
		eigenvalueOut.close();
		printf("%12s%30s\n"," ",eigenvaluefile.c_str());
		//printing eigenvector to file
		string eigenvectorFile = prefix + "eigVec.dat";
		copyFile("data/eigVec.dat",eigenvectorFile);
		printf("%12s%30s\n"," ",eigenvectorFile.c_str());
		}

} //closing parameter loop

return 0;
}

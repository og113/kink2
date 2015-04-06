/*----------------------------------------------------------------------------------------------------------------------------
	pi
		program to solve boundary value problem on contour BC, and step out to A and D
----------------------------------------------------------------------------------------------------------------------------*/

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
		10 - solving for delta
		11 - printing early
		12 - convergence
		13 - propagating along minkowskian time
		14 - printing output
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
/*----------------------------------------------------------------------------------------------------------------------------
	1. loading options
		- loading options
		- loading closenesses
----------------------------------------------------------------------------------------------------------------------------*/

// loading options
Options opts;
opts.load("optionsP");
//opts.print();

// loading closenesses
Closenesses closenesses;
closenesses.load("closenesses");

/*----------------------------------------------------------------------------------------------------------------------------
	2. getting inputs
		- loading inputs
		- defining timenumber
----------------------------------------------------------------------------------------------------------------------------*/

// loading inputs
Parameters psu;
psu.load("inputsP");
//cout << "input parameters: " << endl;
//psu.print();
//cout << endl;

// defining timenumber
string timenumber;
(argc==2)? timenumber = argv[1] : timenumber = currentDateTime();

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
	Parameters ps(psu);
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
	Filename paramsRunFile = (string)("data/"+timenumber+"inputsP_loop_"+numberToString<uint>(loop));
	ps.save(paramsRunFile);
	
	//printing timenumber and parameters
	printf("%12s%12s\n","timenumber: ",timenumber.c_str());
	ps.print();
	
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
	Check checkDT("1/(time derivative of phi)",closenesses.DT);
	Check checkProfile("phi input calculation",closenesses.Profile);
	
	// do trivial or redundant checks?
	bool trivialChecks = false;
	
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
		omegaM1F = (string)("data/stable/omegaM1_pot_"+numberToString<uint>(ps.pot)+"_N_"+numberToString<uint>(ps.N)\
						+"_L_"+numberToString<double>(ps.L)+".dat");
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
	double negVal;
	if (opts.zmt[0]=='n' || opts.zmx[0]=='n') {
		if (ps.pot==3) {
			Filename eigVecFile = (string)("data/stable/eigVec_pot_3_L_" + numberToString<double>(ps.L) + ".dat");	
			load(eigVecFile,so_simple,negVec); // should automatically interpolate
			negVal = -15.3;
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
						Nb_load = ((eigVecFile.Extras)[2]).second;
					}
				}
			}
			SaveOptions eigVecOpts;
			eigVecOpts.vectorType = SaveOptions::realB;
			eigVecOpts.extras = SaveOptions::coords;
			eigVecOpts.paramsOut = ps;
			eigVecOpts.printMessage = false;
			Parameters pIn(ps);
			pIn.N = stringToNumber<uint>(N_load);
			pIn.Nb = stringToNumber<uint>(Nb_load);
			pIn.NT = 1000; // a fudge so that interpolate realises the vector is only on BC
			load(eigVecFile,eigVecOpts,negVec);
			
			string eigValFile = "data/stable/eigVal_pot_"+numberToString<uint>(ps.pot)+".dat";
			string line, dross;
			double negError;
			ifstream isEigVal;
			isEigVal.open(eigValFile.c_str());
			while (!isEigVal.eof()) getline(isEigVal, line); // getting last line
			istringstream iss(line);
			iss >> dross >> dross >> dross >> negError >> negVal;
			if (negError>1.0)
				cerr << "error in eigenvalue(" << negError << ") is large, consider retrying calculation" << endl;
			isEigVal.close();
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
	double ergZero = 0.0;
	cVec erg(ps.NT);	

	//defining some quantities used to stop the Newton-Raphson loop when action stops varying
	comp action_last = action;
	uint runs_count = 0;
	uint min_runs = 3;

	//initializing phi (=p), DDS and minusDS
	vec p(2*ps.N*ps.Nb+1);
	p = Eigen::VectorXd::Zero(2*ps.N*ps.Nb+1);
	spMat DDS(2*ps.N*ps.Nb+1,2*ps.N*ps.Nb+1);
	vec minusDS(2*ps.N*ps.Nb+1);

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
							opts.minLoopLoad+".dat");
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
		if (ps.pot==3) {
			vec tempPhi, temp2Phi;
			Filename pi3GuessFile = (string)("data/sphaleronPi_L_"+numberToString<double>(ps.L)\
										+"_Tb_"+numberToString<double>(ps.Tb)+".dat");
			so_simple.column = 4;
			Parameters ps_blank;
			so_simple.paramsIn = ps_blank;
			so_simple.paramsOut = ps_blank;
			load(pi3GuessFile,so_simple,tempPhi);
			so_simple.paramsIn = ps;
			so_simple.paramsOut = ps;
			so_simple.column = 0;
			uint length = tempPhi.size();
			length = (uint)(sqrt(length));
			ps_blank.N = length;
			ps_blank.Nb = length;
			ps_blank.NT = length*4; // same fudge as before, so that it realises that the vector is just on BC
			temp2Phi = interpolateReal(tempPhi,ps_blank,ps);
			for (uint j=0; j<ps.N*ps.Nb; j++) {
				p(2*j) = temp2Phi(j);
				p(2*j+1) = 0.0;
			}
			p(2*ps.Nb*ps.N) = 0.5; // lagrange multiplier for zero mode
		}
		else {
			//finding phi profile between minima
			uint profileSize = ps.Nb; //more than the minimum
			vector<double> phiProfile(profileSize);
			vector<double> rhoProfile(profileSize);
			double alphaL = opts.alpha, alphaR = opts.alpha;
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
			if (ps.R<opts.alpha) {
				cerr << "R is too small. Not possible to give thinwall input. It should be more that " << opts.alpha;
				return 1;
			}
			for (uint j=0; j<ps.N*ps.Nb; j++) {
				comp t = coordB(j,0,ps);
				double x = real(coordB(j,1,ps));
				p(2*j+1) = 0.0; //imaginary parts set to zero
				if ((opts.inF).compare("b")==0 || ((opts.inF).compare("p")==0 && ps.Tb>ps.R)) {
					double rho = real(sqrt(-pow(t,2.0) + pow(x,2.0))); //should be real even without real()
					if (ps.pot==1) {
						if ((rho-ps.R)<-opts.alpha) 		p(2*j) = ps.minima[1];
						else if ((rho-ps.R)>opts.alpha) 	p(2*j) = ps.minima[0];
						else							p(2*j) = (ps.minima[1]+ps.minima[0])/2.0\
														 + (ps.minima[0]-ps.minima[1])*tanh((rho-ps.R)/2.0)/2.0;
					}
					else if (ps.pot==2) {
						if ((rho-ps.R)<=alphaL) 		p(2*j) = ps.minima[1];
						else if ((rho-ps.R)>=alphaR) 	p(2*j) = ps.minima[0];
						else {
							vector<double> rhoPos (profileSize,rho-ps.R);
							for (uint k=0; k<profileSize; k++) {
								rhoPos[k] -= rhoProfile[k];
							}
							uint minLoc = smallestLoc(rhoPos);
			                p(2*j) = phiProfile[minLoc];
						}
					}
					if ((opts.inF).compare("p")==0) {
						p(2*j) 	 += opts.amp*cos(sqrt(-negVal)*imag(t))*negVec(2*j);
						p(2*j+1) += opts.amp*cos(sqrt(-negVal)*imag(t))*negVec(2*j+1);
					}
				}
				else if ((opts.inF).compare("p")==0 && ps.Tb<ps.R) {
					double angle = asin(ps.Tb/ps.R);
					double rho1 = real(sqrt(-pow(t,2.0) + pow(x+ps.R*cos(angle),2.0)));
					double rho2 = real(sqrt(-pow(t,2.0) + pow(x-ps.R*cos(angle),2.0)));
					if ((rho1-ps.R)<-opts.alpha && (rho2-ps.R)<-opts.alpha)
									p(2*j) = ps.minima[1];
					else if ((rho1-ps.R)>opts.alpha || (rho2-ps.R)>opts.alpha)	
									p(2*j) = ps.minima[0];
					else if (x>0)									
									p(2*j) = (ps.minima[1]+ps.minima[0])/2.0\
												+ (ps.minima[0]-ps.minima[1])*tanh((rho1-ps.R)/2.0)/2.0;
					else if (x<0)									
									p(2*j) = (ps.minima[1]+ps.minima[0])/2.0 \
												+ (ps.minima[0]-ps.minima[1])*tanh((rho2-ps.R)/2.0)/2.0;
					else			p(2*j) = ps.minima[1]; //i.e. if coordB(j,1) == 0
				}
			}
		}
		p(2*ps.N*ps.Nb) = 0.5; // lagrange multiplier for zero mode
	}
	
	// fixing boundary conditions of input phi
	for (uint j=0;j<ps.N;j++) {
	    p(2*j*ps.Nb) = (1.0-opts.open)*p(2*j*ps.Nb) + opts.open*p(2*(j*ps.Nb+1)); //initial time real
	    p(2*(j*ps.Nb+1)) = p(2*j*ps.Nb);
	    p(2*j*ps.Nb+1) = (1.0-opts.open)*p(2*j*ps.Nb+1) + opts.open*p(2*(j*ps.Nb+1)+1); //initial time imag
	    p(2*(j*ps.Nb+1)+1) = p(2*j*ps.Nb+1);
	    p(2*((j+1)*ps.Nb-1)) = opts.open*p(2*((j+1)*ps.Nb-2)) + (1.0-opts.open)*p(2*((j+1)*ps.Nb-1)); //final time real
	    p(2*((j+1)*ps.Nb-2)) = p(2*((j+1)*ps.Nb-1));
	    p(2*((j+1)*ps.Nb-2)+1) = opts.open*p(2*((j+1)*ps.Nb-1)+1) + (1.0-opts.open)*p(2*((j+1)*ps.Nb-2)+1); //final time imag
	    p(2*((j+1)*ps.Nb-1)+1) = p(2*((j+1)*ps.Nb-2)+1);
	}
    if (ps.pot==3) {
		for (uint j=0;j<ps.Nb;j++) {
			uint l0 = c(j,0,ps.Nb);
			uint m = c(j,ps.N-1,ps.Nb);
			p(2*l0) = 0.0; // p=0 ar r=0
			p(2*l0+1) = 0.0;
		    p(2*m) = 0.0; //p=0 ar r=R
		    p(2*m+1) = 0.0;
		}
	}
	
	// printing input phi
	Filename earlyPhiFile = (string)("data/"+timenumber+"pE_loop_"+numberToString<uint>(loop)+"_run_0.dat");
	save(earlyPhiFile,so_p,p);

	// defining complex phi
	cVec Cp(ps.N*ps.Nb);
	Cp = vecComplex(p,ps.N*ps.Nb);
	
/*----------------------------------------------------------------------------------------------------------------------------
	8. beginning newton-raphson loop
		- chiX, for fixing zero mode
		- reserving memory for DDS
		- initializing erg etc to zero
		- trivial potential test
----------------------------------------------------------------------------------------------------------------------------*/

	//beginning newton-raphson loop
	while (!checkSoln.good() || !checkSolnMax.good() || runs_count<min_runs) {
		runs_count++;
		
		//defining the zero mode at the final time boundary and the time step before
		vec chiX(ps.Nb*ps.N);
		chiX = Eigen::VectorXd::Zero(ps.N*ps.Nb);
		if (ps.pot!=3) {
			for (uint j=0; j<ps.N; j++){
				lint pos = c(ps.Nb-1,j,ps.Nb);
				long int neighPos = neigh(pos,1,1,ps.Nb,ps.N,ps.pot), neighMin = neigh(pos,1,-1,ps.Nb,ps.N,ps.pot);
				if(neighPos!=-1 && neighMin!=-1) {
		        	chiX(pos) = p(2*neighPos)-p(2*neighMin); 								//final time slice
		        	//chiX(pos-1) = p(2*neigh(pos-1,1,1,ps.Nb,ps.N))-p(2*neigh(pos-1,1,-1,ps.Nb,ps.N,ps.pot)); //penultimate time slice
		        }
		    }
		}

		// allocating memory for DS, DDS
		minusDS = Eigen::VectorXd::Zero(2*ps.N*ps.Nb+1); //initializing to zero
		DDS.setZero(); //just making sure
		Eigen::VectorXi DDS_to_reserve(2*ps.N*ps.Nb+1);//number of non-zero elements per column
		DDS_to_reserve = Eigen::VectorXi::Constant(2*ps.N*ps.Nb+1,11);
		DDS_to_reserve(0) = 3; //these need to be changed when boundary conditions need to be more compicated
		DDS_to_reserve(1) = 3;
		DDS_to_reserve(2*ps.N*ps.Nb-2) = 3;
		DDS_to_reserve(2*ps.N*ps.Nb-1) = 3;
		DDS_to_reserve(2*ps.N*ps.Nb) = ps.N;
		DDS.reserve(DDS_to_reserve);
		
		//initializing to zero
		comp kineticS = 0.0;
		comp kineticT = 0.0;
		comp pot_0 = 0.0;
		comp pot_r = 0.0;
		erg = Eigen::VectorXcd::Constant(ps.NT,-ergZero);
		double dtTest = 0.0;
		
		//testing that the potential term is working for pot3
		if (ps.pot==3 && trivialChecks) {
			comp Vtrial = 0.0, Vcontrol = 0.0;
			for (uint j=0; j<ps.N; j++) {
				double r = ps.r0 + j*ps.a;
				paramsV.epsi = r;
				V.setParams(paramsV);
				Vcontrol += pow(p(2*j*ps.Nb),2.0)/2.0 - pow(p(2*j*ps.Nb),4.0)/4.0/pow(r,2.0);
				Vtrial += V(p(2*j*ps.Nb));
			}
			double potTest = pow(pow(real(Vcontrol-Vtrial),2.0) + pow(imag(Vcontrol-Vtrial),2.0),0.5);
			cout << "potTest = " << potTest << endl;
		}
/*----------------------------------------------------------------------------------------------------------------------------
	9. assigning minusDS, DDS etc
		- beginning loop over lattice points
		- fixing zero mode
		- t=(NT-1)
		- t=0
		- bulk
		- extras (*4.0*pi, etc)
----------------------------------------------------------------------------------------------------------------------------*/

		// beginning loop over lattice points
		for (lint j = 0; j < ps.N*ps.Nb; j++) {		
			uint t = intCoord(j,0,ps.Nb); //coordinates
			uint x = intCoord(j,1,ps.Nb);
			long int neighPosX = neigh(j,1,1,ps.Nb,ps.N,ps.pot);
			
			if (ps.pot==3) {
					paramsV.epsi = ps.r0+x*ps.a;
					V.setParams(paramsV);
					dV.setParams(paramsV);
					ddV.setParams(paramsV);
			}
			
			if (t<(ps.Nb-1)) dtTest += abs((p(2*(j+1))-p(2*j))/ps.b);
			
			// fixing zero mode
			if ((abs(chiX(j))>1.0e-16)) {
				DDS.insert(2*j,2*ps.N*ps.Nb) = ps.a*chiX(j); 
				DDS.insert(2*ps.N*ps.Nb,2*j) = ps.a*chiX(j);
				minusDS(2*j) += -ps.a*chiX(j)*p(2*ps.N*ps.Nb);
				minusDS(2*ps.N*ps.Nb) += -ps.a*chiX(j)*p(2*j);
		    }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//boundaries
			if (ps.pot==3 && x==(ps.N-1)) {
				DDS.insert(2*j,2*j) = 1.0; // p=0 at r=R
				DDS.insert(2*j+1,2*j+1) = 1.0;
			}
			else if (ps.pot==3 && x==0) {
				DDS.insert(2*j,2*j) = 1.0; // p=0 at r=0
				DDS.insert(2*j+1,2*j+1) = 1.0;
				double csi = ((t==0 || t==(ps.NT-1))? 0.5: 1.0);
				kineticS +=	csi*ps.b*pow(Cp(neighPosX),2.0)/ps.a/2.0;
			}
			else if (t==(ps.Nb-1)) {
				comp Dt(0.0,-ps.b/2.0);
				if (neighPosX!=-1) {
					erg(t+ps.Na) 		+= pow(Cp(neighPosX)-Cp(j),2.0)/ps.a/2.0;
					kineticS			+= Dt*pow(Cp(neighPosX)-Cp(j),2.0)/ps.a/2.0;
				}
				erg(t+ps.Na) 			+= ps.a*V(Cp(j)) + ps.a*Vr(Cp(j));
				pot_0 					+= Dt*ps.a*V(Cp(j));
				pot_r 					+= Dt*ps.a*Vr(Cp(j));				
				
				DDS.insert(2*j,2*j) 	= 1.0/ps.b; //zero time derivative
				DDS.insert(2*j,2*(j-1)) = -1.0/ps.b;
				DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
			}
			else if (t==0) {
				comp dt(0.0,-ps.b);
				comp Dt(0.0,-ps.b/2.0);
				kineticT += ps.a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
				if (neighPosX!=-1) {
					kineticS 		+= Dt*pow(Cp(neighPosX)-Cp(j),2.0)/ps.a/2.0;
					erg(t+ps.Na) 	+= pow(Cp(neighPosX)-Cp(j),2.0)/ps.a/2.0;
				}
				
				pot_0 += Dt*ps.a*V(Cp(j));
				pot_r += Dt*ps.a*Vr(Cp(j));
				erg(t+ps.Na) += ps.a*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + ps.a*V(Cp(j)) + ps.a*Vr(Cp(j));
				
				if ((opts.inF).compare("b")==0) {
					DDS.insert(2*j,2*j) = 1.0; //zero change (initial input satisfies b.c.s)
					DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				}
				else if ((opts.inF).compare("p")==0) {
					DDS.insert(2*j,2*j) = -1.0/ps.b; //zero time derivative
					DDS.insert(2*j,2*(j+1)) = 1.0/ps.b;
					DDS.insert(2*j+1,2*j+1) = 1.0; //zero imaginary part
				}
			}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//bulk
			else {
				comp dt(0.0,-ps.b);
				comp Dt(0.0,-ps.b);
				kineticT += ps.a*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
				if (neighPosX!=-1) {
					kineticS += Dt*pow(Cp(neighPosX)-Cp(j),2.0)/ps.a/2.0;
					erg(t+ps.Na) 	+= pow(Cp(neighPosX)-Cp(j),2.0)/ps.a/2.0;
				}
				pot_0 += Dt*ps.a*V(Cp(j));
				pot_r += Dt*ps.a*Vr(Cp(j));
				erg(t+ps.Na) += ps.a*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + ps.a*V(Cp(j)) + ps.a*Vr(Cp(j));
				
                for (uint k=0; k<2*2; k++) {
                    int sign = pow(-1,k);
                    uint direc = (uint)(k/2.0);
                    long int neighb = neigh(j,direc,sign,ps.Nb,ps.N,ps.pot);
                    if (direc == 0) {
                        minusDS(2*j) += real(ps.a*Cp(j+sign)/dt);
                        minusDS(2*j+1) += imag(ps.a*Cp(j+sign)/dt);
                        DDS.insert(2*j,2*(j+sign)) = -real(ps.a/dt);
                        DDS.insert(2*j,2*(j+sign)+1) = imag(ps.a/dt);
                        DDS.insert(2*j+1,2*(j+sign)) = -imag(ps.a/dt);
                        DDS.insert(2*j+1,2*(j+sign)+1) = -real(ps.a/dt);
                       }
                    else if (neighb !=-1) {
                        minusDS(2*j) += - real(Dt*Cp(neighb)/ps.a);
                        minusDS(2*j+1) += - imag(Dt*Cp(neighb)/ps.a);
                        DDS.insert(2*j,2*neighb) = real(Dt/ps.a);
                        DDS.insert(2*j,2*neighb+1) = -imag(Dt/ps.a);
                        DDS.insert(2*j+1,2*neighb) = imag(Dt/ps.a);
                        DDS.insert(2*j+1,2*neighb+1) = real(Dt/ps.a);
                    }
                }
                comp temp0 = 2.0*ps.a/dt;
	            comp temp1 = ps.a*Dt*(2.0*Cp(j)/pow(ps.a,2.0) + dV(Cp(j)) + dVr(Cp(j)));
	            comp temp2 = ps.a*Dt*(2.0/pow(ps.a,2.0) + ddV(Cp(j)) + ddVr(Cp(j)));
	                
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
        	DDS.insert(2*ps.N*ps.Nb,2*ps.N*ps.Nb) = 1.0;	// redundant lagrange multiplier
        	dtTest /= (double)(ps.Nb-1.0);
        	dtTest = 1.0/dtTest;
        	action 	*= 4.0*pi;
        	erg 	*= 4.0*pi;
        }
        
/*----------------------------------------------------------------------------------------------------------------------------
	9. solving for delta
		- defining delta
		- analyzing pattern
		- factorizing
		- solving
		- check on inversion
		- p' = p + delta
		- Cp'
----------------------------------------------------------------------------------------------------------------------------*/

	//solving for delta in DDS*delta=minusDS, where p' = p + delta		
		vec delta(2*ps.N*ps.Nb+1);
		delta = Eigen::VectorXd::Zero(2*ps.N*ps.Nb+1);
		DDS.makeCompressed();
		Eigen::SparseLU<spMat> solver;
		
		solver.analyzePattern(DDS);
		if(solver.info()!=Eigen::Success) {
			cerr << "DDS pattern analysis failed, solver.info() = "<< solver.info() << endl;
			return 1;
		}		
		solver.factorize(DDS);
		if(solver.info()!=Eigen::Success) {
			cerr << "Factorization failed, solver.info() = "<< solver.info() << endl;
			return 1;
		}
		delta = solver.solve(minusDS);// use the factorization to solve for the given right hand side
		if(solver.info()!=Eigen::Success) {
			cerr << "Solving failed, solver.info() = "<< solver.info() << endl;
			cerr << "log(abs(det(DDS))) = " << solver.logAbsDeterminant() << endl;
			cerr << "sign(det(DDS)) = " << solver.signDeterminant() << endl;
			return 1;
		}
		
		//independent check on whether calculation worked
		vec diff(2*ps.N*ps.Nb+1);
		diff = DDS*delta-minusDS;
		double maxDiff = diff.maxCoeff();
		maxDiff = abs(maxDiff);
		checkInv.add(maxDiff);
		checkInv.checkMessage();
		if (!checkInv.good()) {
			return 1;
		}

		//assigning values to phi
		p += delta;
		
		//passing changes on to complex vector
		Cp = vecComplex(p,ps.N*ps.Nb);
		
/*----------------------------------------------------------------------------------------------------------------------------
	10. printing early
		- p
		- minusDS
		- DDS
		- delta
		- chiX
		- erg
----------------------------------------------------------------------------------------------------------------------------*/
		
		//printing early if desired	
		if ((opts.printChoice).compare("n")!=0) {
			Filename basic = (string)("data/"+timenumber+"basic_loop_"+numberToString<uint>(loop)\
								+"_run_"+numberToString<uint>(runs_count)+".dat");
			if ((opts.printChoice).compare("p")==0 || (opts.printChoice).compare("e")==0) {
				Filename pEFile = basic;
				pEFile.ID = "pE";
				save(pEFile,so_p,p);
			}
			if ((opts.printChoice).compare("v")==0 || (opts.printChoice).compare("e")==0) {
				Filename vEFile = basic;
				vEFile.ID = "minusDSE";
				save(vEFile,so_p,minusDS);
			}
			if ((opts.printChoice).compare("m")==0 || (opts.printChoice).compare("e")==0) {
				Filename mEFile = basic;
				mEFile.ID = "DDSE";
				save(mEFile,so_simple,DDS);
			}
			if ((opts.printChoice).compare("d")==0 || (opts.printChoice).compare("e")==0) {
				Filename dEFile = basic;
				dEFile.ID = "deltaE";
				save(dEFile,so_p,delta);
			}
			if ((opts.printChoice).compare("z")==0 || (opts.printChoice).compare("e")==0) {
				Filename zEFile = basic;
				zEFile.ID = "chiXE";
				save(zEFile,so_simple,chiX);
			}
			if ((opts.printChoice).compare("erg")==0 || (opts.printChoice).compare("e")==0) {
				Filename lEFile = basic;
				lEFile.ID = "ergE";
				save(lEFile,so_simple,erg);
			}
		}
		
/*----------------------------------------------------------------------------------------------------------------------------
	10. convergence issues
		- checkReg
		- evaluating norms
		- evaluating checks to stop n-r loop
		- printing convergence tests
		- end of n-r loop
		- checkDT
----------------------------------------------------------------------------------------------------------------------------*/

		// checking pot_r is much smaller than the other potential terms
		checkReg.add(abs(pot_r/pot_0));
		checkReg.checkMessage();
		
		// evaluating norms
		double normDS = minusDS.dot(minusDS);
		normDS = pow(normDS,0.5);
		double maxDS = minusDS.maxCoeff();
		double minDS = minusDS.minCoeff();
		if (-minDS>maxDS) maxDS = -minDS;
		double normP = p.dot(p);
		normP = pow(normP,0.5);
		double normDelta = delta.dot(delta);
		normDelta = pow(normDelta,0.5);
		
		// evaluating checks to stop n-r loop
		checkAction.add(abs(action - action_last)/abs(action_last));
		action_last = action;
		checkSoln.add(normDS/normP);
		checkSolnMax.add(maxDS);
		checkDelta.add(normDelta/normP);
		checkDT.add(dtTest);
			
		// printing tests to see convergence
		if (runs_count==1) {
			printf("%16s%16s%16s%16s%16s%16s\n","loop","runsCount","actionTest","solTest","solMTest","deltaTest");
		}
		printf("%16i%16i%16g%16g%16g%16g\n",loop,runs_count,checkAction.back(),checkSoln.back()\
									,checkSolnMax.back(),checkDelta.back());
		
		} //closing newton-raphson loop
		
		checkDT.checkMessage();
		
/*----------------------------------------------------------------------------------------------------------------------------
	11. propagating solution along minkowskian time
		- B->A
			- initialise ap, velA, accA
			- initialise linErgA, linNumA
			- run B->A loop
			- evaluate E, N and W
		- C->D
			- initialise cp, velC, accC
			- run C->D loop 
		- check energy conserved
		- form tp
		- print tp, linErg, linNum
		
----------------------------------------------------------------------------------------------------------------------------*/
	
	//propagating solution back in minkowskian time
	if (ps.pot!=3) {
		//A1. initialize mp==mphi using last point of ephi and zeros- use complex phi
		cVec ap(ps.N*(ps.Na+1)); //phi on section "a"
		ap = Eigen::VectorXcd::Zero(ps.N*(ps.Na+1));
		for (uint j=0; j<ps.N; j++) ap(j*(ps.Na+1)) = Cp(j*ps.Nb);

		//A2. initialize vel - defined at half steps, first step being at t=-1/2,
		//vel(t+1/2) := (p(t+1)-p(t))/dt
		cVec velA (ps.N*(ps.Na+1));
		velA = Eigen::VectorXcd::Zero(ps.N*(ps.Na+1));
		double dtau = -ps.b;
		double Dt0 = dtau; //b/2*(-1+1i); - this is surely wrong!!
		
		//A3. initialize acc using phi and expression from equation of motion and zeros-complex
		cVec accA(ps.N*(ps.Na+1));
		accA = Eigen::VectorXcd::Zero(ps.N*(ps.Na+1));
		//#pragma omp parallel for
		for (uint j=0; j<ps.N; j++) {
			lint l = j*(ps.Na+1);
			accA(l) = ((Dt0/pow(ps.a,2.0))*(ap(neigh(l,1,1,ps.Na+1,ps.N,ps.pot))\
												+ap(neigh(l,1,-1,ps.Na+1,ps.N,ps.pot))-2.0*ap(l))-Dt0*(dV(ap(l))+dVr(ap(l))))/dtau;
		}
			
		//A4.5 starting the energy and that off
		vec linErgA(ps.Na); linErgA = Eigen::VectorXd::Zero(ps.Na);
		vec linNumA(ps.Na); linNumA = Eigen::VectorXd::Zero(ps.Na);

		//A7. run loop
		for (uint t=1; t<(ps.Na+1); t++) {
		    for (uint x=0; x<ps.N; x++) {
		        lint m = t+x*(ps.Na+1);
		        velA(m) = velA(m-1) + dtau*accA(m-1);
		        ap(m) = ap(m-1) + dtau*velA(m);
		    }
		    for (uint x=0; x<ps.N; x++) {
		        uint m = t+x*(ps.Na+1);
				accA(m) = (1.0/pow(ps.a,2.0))*(ap(neigh(m,1,1,ps.Na+1,ps.N,ps.pot))+ap(neigh(m,1,-1,ps.Na+1,ps.N,ps.pot))-2.0*ap(m)) \
	        		-dV(ap(m)) - dVr(ap(m));
	        	erg(ps.Na-t) += ps.a*pow(ap(m-1)-ap(m),2.0)/pow(-dtau,2.0)/2.0\
	        		 + pow(ap(neigh(m,1,1,ps.Na+1,ps.N,ps.pot))-ap(m),2.0)/ps.a/2.0 + ps.a*V(ap(m)) + ps.a*Vr(ap(m));
		        for (uint y=0; y<ps.N; y++) {
		        	lint n = t + y*(ps.Na+1);
				    	linErgA(ps.Na-t) += omega_2(x,y)*(real(ap(m))-ps.minima[0])*(real(ap(n))-ps.minima[0])\
				    						 + omega_2(x,y)*imag(ap(m))*imag(ap(n));
						linNumA (ps.Na-t) += omega_1(x,y)*(real(ap(m))-ps.minima[0])*(real(ap(n))-ps.minima[0])\
											 + omega_1(x,y)*imag(ap(m))*imag(ap(n));
				}
			}
		}
		
		if (ps.pot==3) {
			linErgA *= 4.0*pi;
			linNumA *= 4.0*pi;
		}
		E = linErgA(0);
		W = - E*2.0*ps.Tb + 2.0*imag(action);

		//now propagating forwards along c
		//C2. initialize mp==mphi using last point of ephi and zeros- use complex phi
		cVec ccp(ps.N*(ps.Nc+1)); //phi on section "c"
		ccp = Eigen::VectorXcd::Zero(ps.N*(ps.Nc+1));
		for (uint j=0; j<ps.N; j++) ccp(j*(ps.Nc+1)) = Cp(j*ps.Nb+ps.Nb-1);

		//C3. initialize vel - defined at half steps, first step being at t=-1/2,
		//vel(t+1/2) := (p(t+1)-p(t))/dt
		cVec velC (ps.N*(ps.Nc+1));
		velC = Eigen::VectorXcd::Zero(ps.N*(ps.Nc+1));
		dtau = ps.b;
		Dt0 = dtau; //b/2*(1-1i); - this is surely wrong!!

		//C4. initialize acc using phi and expression from equation of motion and zeros-complex
		cVec accC(ps.N*(ps.Nc+1));
		accC = Eigen::VectorXcd::Zero(ps.N*(ps.Nc+1));
		for (uint j=0; j<ps.N; j++){
			lint l = j*(ps.Nc+1);
			accC(l) = ((Dt0/pow(ps.a,2.0))*(ccp(neigh(l,1,1,ps.Nc+1,ps.N,ps.pot))+ccp(neigh(l,1,-1,ps.Nc+1,ps.N,ps.pot))-2.0*ccp(l))\
						-Dt0*(dV(ccp(l))+dVr(ccp(l))))/dtau;
		}

		//C7. run loop
		for (uint t=1; t<(ps.Nc+1); t++) {
			for (uint x=0; x<ps.N; x++) {
				lint l = t+x*(ps.Nc+1);
				velC(l) = velC(l-1) + dtau*accC(l-1);
				ccp(l) = ccp(l-1) + dtau*velC(l);
			}
			for (uint x=0; x<ps.N; x++) {
				lint l = t+x*(ps.Nc+1);
				accC(l) = (1.0/pow(ps.a,2.0))*(ccp(neigh(l,1,1,ps.Nc+1,ps.N,ps.pot))+ccp(neigh(l,1,-1,ps.Nc+1,ps.N,ps.pot))-2.0*ccp(l))\
								-dV(ccp(l));		    	
				erg (ps.Na+ps.Nb-2+t) += ps.a*pow(ccp(l)-ccp(l-1),2.0)/pow(dtau,2.0)/2.0\
				 	+ pow(ccp(neigh(l-1,1,1,ps.Nc+1,ps.N,ps.pot))-ccp(l-1),2.0)/ps.a/2.0\
					+ ps.a*V(ccp(l-1)) + ps.a*Vr(ccp(l-1));
			}
		}
		
		//checking energy conserved
		double relErgChange = 0.0;
		if (abs(real(erg(0)))>MIN_NUMBER) {
			relErgChange = absDiff(erg(0),erg(ps.NT-2));
		}
		checkCon.add(relErgChange);
		checkCon.checkMessage();

		//12. combine phi with ap and cp and save combination to file
		cVec tCp(ps.NT*ps.N);
		//#pragma omp parallel for
		for (uint j=0; j<ps.NT*ps.N; j++) {
		    uint t = intCoord(j,0,ps.NT);
		    uint x = intCoord(j,1,ps.NT);
		    if (t<ps.Na) {
		        t = ps.Na-t;
		        tCp(j) = ap(t+x*(ps.Na+1));
		       }
		    else if (t<(ps.Na+ps.Nb)) {
		        t = t - ps.Na;
		        tCp(j) = Cp(t+x*ps.Nb);
		       }
		    else {
		        t = t - ps.Na - ps.Nb + 1;
		        tCp(j) = ccp(t+x*(ps.Nc+1));
		    }
		}
			
		//making real vec from complex one
		vec tp(2*ps.N*ps.NT);
		tp = vecReal(tCp,ps.NT*ps.N);
		tp.conservativeResize(2*ps.N*ps.NT+1);
		tp(2*ps.N*ps.NT) = p(2*ps.N*ps.Nb);
		
		string prefix = "data/"+timenumber;
		string suffix = "_loop_"+numberToString<uint>(loop)+".dat";
		
		// printing tp
		Filename tpFile = (string)(prefix+"tp"+suffix);
		SaveOptions so_tp = so_p;
		so_tp.vectorType = SaveOptions::complex;
		so_tp.printMessage = true;
		save(tpFile,so_tp,tp);
		
		//printing linErgA
		Filename linErgFile = (string)(prefix+"linErg"+suffix);
		so_simple.printMessage = true;
		save(linErgFile,so_simple,linErgA);
    }
    else {
    	E = 0.0;
    	W = 0.0;
    }
    
/*----------------------------------------------------------------------------------------------------------------------------
	14. printing output
		- stopping clock
		- printing results to terminal
		- printing results to file
		- printing (and plotting if vectors):
			- p
			- minusDS
			- DDS
			- erg
		- return
		
----------------------------------------------------------------------------------------------------------------------------*/
    
    //stopping clock
	time = clock() - time;
	double realtime = time/1000000.0;
	
	//printing results to terminal
	printf("\n");
	printf("%8s%8s%8s%8s%8s%8s%8s%12s%12s%12s\n","runs","time","N","NT","L","Tb","dE","E","im(S)","W");
	printf("%8i%8g%8i%8i%8g%8g%8g%12.4g%12.4g%12.4g\n",runs_count,realtime,ps.N,ps.NT,ps.L,ps.Tb,ps.dE,E,imag(action),W);
	printf("\n");
	printf("%60s\n","----------------------------------------------------------------------------------------------------");

	//printing results to file
	FILE * actionfile;
	actionfile = fopen("./data/action.dat","a");
	fprintf(actionfile,"%16s%8i%8i%8g%8g%8g%10.4g%10.4g%10.4g%10.4g%10.4g\n",timenumber.c_str()\
				,ps.N,ps.NT,ps.L,ps.Tb,ps.dE,E,imag(action)\
				,W,checkSoln.back(),checkCon.back());
	fclose(actionfile);
	
	string prefix = "./data/"+timenumber;
	string suffix = "_loop_"+numberToString<uint>(loop)+".dat";
	
	// plot options
	PlotOptions po_p;
	po_p.gp = "gp/repi.gp";
	po_p.style = "points";
	Filename plotFile = (string)("data/"+timenumber+"p_loop_"+numberToString<uint>(loop)+".png");
	po_p.output = plotFile;
	po_p.printMessage = true;
	
	PlotOptions po_simple;
	po_simple.column = 1;
	po_simple.style = "linespoints";
	po_simple.printMessage = true;

	//printing output phi on Euclidean time part
	so_p.printMessage = true;
	Filename pFile = (string)(prefix+"p"+suffix);
	save(pFile,so_p,p);
	plot(pFile,po_p);
	
	//printing output minusDS
	pFile.ID = "minudDS";
	save(pFile,so_p,minusDS);
	//plotFile = pFile;
	//plotFile.Suffix = ".png";
	//po_p.output = plotFile;
	//plot(pFile,po_p);
			
	//printing output DDS
	so_simple.printMessage = true;
	pFile.ID = "DDS";
	save(pFile,so_simple,DDS);
	
	//printing erg
	pFile.ID = "erg";
	//erg.conservativeResize(ps.Na);
	save(pFile,so_simple,erg);
	plotFile = pFile;
	plotFile.Suffix = ".png";
	po_simple.output = plotFile;
	plot(pFile,po_simple);

	if (!checkDelta.good()) {
		return 1;
	}

} //closing parameter loop

return 0;
}

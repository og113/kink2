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
#include "parameters.h"
#include "lattice.h"
#include "folder.h"
#include "check.h"
#include "print.h"
#include "omega.h"
#include "main.h"

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - loading options
		2 - Folders
		3 - beginning file loop
		4 - beginning parameter loop
		5 - assigning potential functions
		6 - omega and negVec
		7 - defining quantitites
		8 - beginning newton-raphson loop
		9 - assigning minusDS, DDS etc
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main()
{
/*----------------------------------------------------------------------------------------------------------------------------
	1. loading options
----------------------------------------------------------------------------------------------------------------------------*/

Options opts;
opts.load("options");

/*----------------------------------------------------------------------------------------------------------------------------
	2. Folders
		- FilenameAttributes for defining FilenameComparator
		- FilenameComparator
		- pFolder
		- inputsFolder
		- removeUnshared
		- printing folders
		- defining timenumber
		
----------------------------------------------------------------------------------------------------------------------------*/

// FilenameAttributes for defining FilenameComparator
FilenameAttributes fa_low, fa_high;
fa_low.Timenumber = opts.minTimenumberLoad;
fa_high.Timenumber = opts.minTimenumberLoad;
(fa_low.Extras).push_back(stringPair("Loop",opts.minLoopLoad);
(fa_high.Extras).push_back(stringPair("Loop",opts.maxLoopLoad);

// FilenameComparator
FilenameComparator fc(fa_low,fa_high);
Folder allFiles(fc);

// pFolder
if ((opts.inF).compare("p")==0) {
	fa_low.ID = "tpip";
	fa_high.ID = "tpip";
}
else if ((opts.inF).compare("m")==0) {
	fa_low.ID = "mainpi";
	fa_high.ID = "mainpi";
}
else {
	cerr << "inF error: " << opts.inF << " not recognised" << endl;
	return 1;
}
fc(fa_low,fa_high);
Folder pFolder(fc);

// inputsFolder
if ((opts.inF).compare("p")==0) {
	fa_low.ID = "inputsPi";
	fa_high.ID = "inputsPi";
}
else if ((opts.inF).compare("m")==0) {
	fa_low.ID = "inputsM";
	fa_high.ID = "inputsM";
}
fc(fa_low,fa_high);
Folder inputsFolder(fc);

// removeUnshared
removeUnshared(pFolder,inputsFolder);

// printing folders
cout << "inputs: " << pFolder << inputsFolder << endl;

//defining timenumber
string timenumber = currentDateTime();

/*----------------------------------------------------------------------------------------------------------------------------
	3. beginning file loop
		- beginning file loop
		- loading and printing inputs
----------------------------------------------------------------------------------------------------------------------------*/

// beginning file loop
for (unsigned int fileLoop=0; fileLoop<pFolder.size(); fileLoop++) {

	// loading parameters
	Parameters psu;
	psu.load(inputsFolder[fileLoop]);
	cout << "input parameters: " << endl;
	psu.print();
	cout << endl;

/*----------------------------------------------------------------------------------------------------------------------------
	4. beginning parameter loop
		- changing parameters (if required)
		- copying a verson of parameters with timenumber
		- declaring checks
----------------------------------------------------------------------------------------------------------------------------*/
	
	for (unsigned int loop=0; loop<opts.loops; loop++) {
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

	//copying a version of ps with timenumber
	Filename paramsRunFile = (string)("./data/" + timenumber + "inputsM_" + numberToString<uint>(loop));
	ps.save(paramsRunFile);
	
	// declaring Checks
	Check checkAction("action",1.0e-2);
	Check checkSoln("solution",1.0e-6);
	Check checkSolnMax("solution max",1.0e-5);
	Check checkDelta("delta",1.0);
	Check checkInv("matrix inversion",1.0e-16*ps.N*ps.NT);
	Check checkCon("energy conservation",1.0e-2);
	Check checkLin("linear energy flat",5.0e-2);
	Check checkTrue("linear energy equal true energy",5.0e-2);
	Check checkLatt("lattice small enough for energy",0.2);
	Check checkReg("regularisation term",1.0e-2);
	Check checkIE("imaginary part of energy",1.0e-5);
	Check checkContm("linear energy equal continuum expression",5.0e-2);
	Check checkOS("linear energy equal on shell expression",5.0e-2);
	Check checkAB("a_k = Gamma*b_k",1.0e-2);
	Check checkABNE("N = Sum(a_k*b_k), E = Sum(w_k*a_k*b_k)",5.0e-2);
	Check checkLR("linear representation of phi",1.0e-12);
	
/*----------------------------------------------------------------------------------------------------------------------------
	5. assigning potential functions
		- typedef PotentialType
		- assigning potential functions
		- assigning preliminary parameter structs
----------------------------------------------------------------------------------------------------------------------------*/
	
	// typedef PotentialType
	typedef comp(*PotentialType)(const comp&, const struct ps_for_V&);

	// assigning potential functions
	Potential V, dV, ddV;
// assigning potential functions
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
	ps_for_V paramsV  = {ps.epsilon, ps.A};
	
	//lambda functions for pot_r
	auto Vr = [&] (const comp& phi) {
		return -ii*reg*VrFn(phi,ps.minima[0],ps.minima[1]);
	};
	auto dVr = [&] (const comp& phi) {
		return -ii*reg*dVrFn(phi,ps.minima[0],ps.minima[1]);
	};
	auto ddVr = [&] (const comp& phi) {
		return -ii*reg*ddVrFn(phi,ps.minima[0],ps.minima[1]);
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
	{
		Filename omegaM1F, omega0F, omega1F, omega2F, modesF, freqsF, freqsExpF; // Filename works as FilenameAttributes
		omegaM1F = (string)"data/stable/omegaM1_pot_"+numberToString<uint>(ps.pot)+"_N_"+numberToString<uint>(ps.N)\
						+"_L_"+numberToString<double>(ps.L)+".dat";
		Folder omegaM1Folder(omegaM1F);
		omega0F = omegaM1F; 		omegaF.ID = "omega0"; 		Folder omega0Folder(omega0F);
		omega1F = omegaM1F; 		omega1F.ID = "omega1"; 		Folder omega1Folder(omega1F);
		omega2F = omegaM1F; 		omega2F.ID = "omega2";	 	Folder omega2Folder(omega2F);
		modesF = omegaM1F;			modes.ID = "modes"; 		Folder modesFolder(modesF);
		freqsF = omegaM1F;			freqs.ID = "freqs"; 		Folder freqsFolder(freqsF);
		freqsExpF = omegaM1F;		omegaM1F.ID = "freqsExp";	Folder freqsExpFolder(omegaM1F);
		SaveOptions so;
		so.ParamsIn = ps;
		so.ParamsOut = ps;
		so.vectorType = SaveOptions::simple; // not sure that this is necessary
		so.extras = SaveOptions::none;
		if (omegaM1Folder.size()==1 && omega0Folder.size()==1 && omega1Folder.size()==1 && omega2Folder.size()>=1 \
			modesFolder.size()==1 && freqsFolder.size()==1 && freqsExpFolder.size()==1) {
			load(omegaM1Folder[0],so,omega_m1);
			load(omega0Folder[0],so,omega_0);
			load(omega1Folder[0],so,omega_1);
			load(omega2Folder[0],so,omega_2);
			load(modesFolder[0],so,modes);
			load(freqsFolder[0],so,freqs);
			load(freqsExpFolder[0],so,freqs_exp);
		}
		else {
			bool approxOmega = Flse;
			if (!approxOmega) {
				numericalModes(modes,freqs,freqs_exp,ps);
			}
			else {
				analyticModes(modes,freqs,freqs_exp,ps);
			}
			omegasFn(approxOmega,modes,freqs,omega_m1,omega_0,omega_1,omega_2,ps);
			save(omegaM1F,so,omega_m1);
			save(omega0F,so,omega_0);
			save(omega1F,so,omega_1);
			save(omega2F,so,omega_2);
			save(modesF,so,modes);
			save(freqsF,so,freqs);
			save(freqsExpF,so,freqs_exp);
		}
	}

	// getting negVec
	vec negVec;
	if (opts.zmt[0]=='n' || opts.zmx[0]=='n') {
		if (ps.pot==3) {
			Filename eigVecFile = "data/stable/eigVec_pot_1_L_" + numberToString<double>(ps.L) + ".dat";	
			SaveOptions eigVecOpts;
			eigVecOpts.vectorType = SaveOptions::simple;
			eigVecOpts.extras = SaveOptions::none;
			eigVecOpts.paramsOut = ps;
			load(eigVecFile,eigVecOpts,&negVec); // should automatically interpolate
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
			Parameters pIn = ps;
			pIn.N = stringToNumber<uint>(N_load);
			pIn.Nb = stringToNumber<uint>(Nb_load);
			pIn.NT = 1000; // a fudge so that interpolate realises the vector is only on BC
			load(eigVecFile,eigVecOpts,&negVec);
		}

/*----------------------------------------------------------------------------------------------------------------------------
	7. defining quantities
		- time
		- erg, linErg etc
		- p, minusDS, DDS
		- print parameters
		- print input phi
----------------------------------------------------------------------------------------------------------------------------*/
		
		//defining a time and starting the clock
		clock_t time;
		time = clock();
	
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
		double ergZero = (ps.pot==3? 0.0: p.N*ps.a*real(V(ps.minima[0])) );
		comp action = ii*ps.action0;
		double bound;
		double W;
		double E;
		double E_exact;
		double Num;
		
		//defining some quantities used to stop the Newton-Raphson loop when action stops varying
		comp action_last = action;
		uint runs_count = 0;
		uint min_runs = 3;

		//initializing phi (=p)
		vec p;
		SaveOptions so_tp;
		so_tp.ParamsIn = psu;
		so_tp.ParamsOut = ps;
		so_tp.vectorType = SaveOptions::complex;
		so_tp.extras = SaveOptions::coords;
		so_tp.zeroModes = 2;
		if (loop==0) {
			load(pFolder[0],so_tp,p);
			printf("%12s%30s\n","input: ",(pFolder[0]()).c_str());
		}
		else {
			Filename lastPhi = (string)("./data/" + timeNumber + "mainpi_fLoop_" + numberToString<uint>(fileLoop) + "_loop_"\
								 + numberToString<uint>(loop-1)+".dat");
			so_tp.ParamsIn = ps;
			load(lastPhi,so_tp,p);
			printf("%12s%30s\n","input: ",(lastPhi()).c_str());
			}
		
		//defining complexified vector Cp
		cVec Cp;
		Cp = vecComplex(p,ps);
	
		//defining DDS and minusDS
		spMat DDS(2*ps.N*ps.NT+2,2*ps.N*ps.NT+2);
		vec minusDS(2*ps.N*ps.NT+2);
			
		//printing loop name and parameters
		printf("%12s%12s\n","timenumber: ",timenumber.c_str());
		printParameters();
			
		//very early vector print
		so_tp.ParamsIn = ps;
		Filename earlyPrintFile = (string)("data/" + timeNumber + "mainpiE_fLoop_" + numberToString<uint>(fileLoop)\
				 + "_loop_" + numberToString<unsigned int>(loop) + "_run_" + "0.dat");
		save(earlyPrintFile,so_tp,p);
	
/*----------------------------------------------------------------------------------------------------------------------------
	8. beginning newton-raphson loop
		- beginning newton-raphson loop	
		- chiT, chiX
		- reserving space in DDS
		- zeroing erg etc
----------------------------------------------------------------------------------------------------------------------------*/
	
		//beginning newton-raphson loop	
		while (!checkSoln.good() || !checkSolnMax.good() || runs_count<min_runs)) {
			runs_count++;
			
			//defining the zero mode at the final time boundary and the time step before
			vec chiX(ps.NT*ps.N);	chiX = Eigen::VectorXd::Zero(ps.N*ps.NT); //to fix spatial zero mode
			vec chiT(ps.NT*ps.N);	chiT = Eigen::VectorXd::Zero(ps.N*ps.NT); //to fix real time zero mode
			for (uint j=0; j<ps.N; j++) {
				uint posX, posT, posCe;
				uint slicesX, slicesT;
				if (getLastInt(opts.zmx)<0) {
					cerr << "getLastInt error with zmx = " << opts.zmx << endl;
					return 1;
				}
				if (getLastInt(opts.zmt)<0) {
					cerr << "getLastInt error with zmt = " << opts.zmt << endl;
					return 1;
				}
				slicesX = getLastInt(zmx);
				slicesT = getLastInt(zmt);
				posCe = j*Nb+Nb-slicesT; //position C for Euclidean vector, i.e. for negVec
				map<char,unsigned int> posMap;
				posMap['A'] = j*NT;
				posMap['B'] = j*NT + (Na-1);
				posMap['C'] = j*NT+Na+Nb-slicesT;
				posMap['D'] = j*NT+(NT-1)-slicesT;
				if (zmt.size()<3) { cout << "zmt lacks info, zmt = " << zmt << endl; }
				for (unsigned int l=0;l<(zmt.size()-2);l++) {
					if (posMap.find(zmt[1+l])!=posMap.end()) {
						posT = posMap.at(zmt[1+l]);
						for (unsigned int k=0;k<slicesT;k++) {
							if (zmt[0]=='n' && ps.pot!=3) 	chiT(posT+k) = negVec(2*(posCe+k));
							else if (zmt[0]=='n' && ps.pot==3) {
								double r = r0 + j*a;
								chiT(posT+k) = negVec(j)*r;
							}
							else if (zmt[0]=='d')				chiT(posT+k) = p(2*(posT+k+1))-p(2*(posT+k));
							else {
								cerr << "choice of zmt not allowed" << endl;
								return 1;
							}
						}
					}
				}
				posMap.erase('C');
				posMap.erase('D');
				posMap['C'] = j*NT+Na+Nb-slicesX;
				posMap['D'] = j*NT+NT-slicesX;
				posCe = j*Nb+Nb-slicesX;
				if (zmx.size()<3) { cout << "zmx lacks info, zmx = " << zmx << endl; }
				for (unsigned int l=0;l<(zmx.size()-2);l++) {
					if (posMap.find(zmx[1+l])!=posMap.end()) {
						posX = posMap.at(zmx[1+l]);
						for (unsigned int k=0;k<slicesX;k++) {
							if (zmx[0]=='n' && ps.pot!=3)		chiX(posX+k) = negVec(2*(posCe+k));
							else if (zmx[0]=='n' && ps.pot==3) {
								double r = r0 + j*a;
								chiX(posX+k) = negVec(j)*r;
							}
							else if (zmx[0]=='d' && ps.pot!=3) chiX(posX+k) = p(2*neigh(posX+k,1,1,NT,N))-p(2*neigh(posX+k,1,-1,NT,N));
							else {
								cerr << "choice of zmx not allowed" << endl;
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
				cerr << "norm of chiX = " << normX << ", norm of chiT = " << normT << endl;
			}
			chiX = chiX/normX;
			chiT = chiT/normT;
			
			// allocating memory for DS, DDS
			minusDS = Eigen::VectorXd::Zero(2*N*NT+2); //initializing to zero
			DDS.setZero(); //just making sure
			Eigen::VectorXi DDS_to_reserve(2*N*NT+2);//number of non-zero elements per column
			DDS_to_reserve = Eigen::VectorXi::Constant(2*N*NT+2,13);
			if (abs(theta)<2.0e-16) {
				DDS_to_reserve(0) = N+3;
				DDS_to_reserve(1) = 11;
			}
			else {
				DDS_to_reserve(0) = N+11;
				DDS_to_reserve(1) = N+11;
			}
			DDS_to_reserve(2*N*NT-2) = 4;
			DDS_to_reserve(2*N*NT-1) = 4;
			DDS_to_reserve(2*N*NT) = N*(zmx.size()-2);
			DDS_to_reserve(2*N*NT+1) = 2*N*(zmt.size()-2);
			DDS.reserve(DDS_to_reserve);
			
			//initializing to zero
			comp kineticS = 0.0;
			comp kineticT = 0.0;
			comp potV = 0.0;
			comp pot_r = 0.0;
			bound = 0.0;
			erg = Eigen::VectorXcd::Constant(NT,-ergZero);
			linErg = Eigen::VectorXcd::Zero(NT);
			linErgOffShell = Eigen::VectorXcd::Zero(NT);
			linNum = Eigen::VectorXcd::Zero(NT);
			linNumOffShell = Eigen::VectorXcd::Zero(NT);
			derivErg = Eigen::VectorXcd::Zero(NT);
			potErg = Eigen::VectorXcd::Constant(NT,-ergZero);
			linErgContm = 0.0;
			linNumContm = 0.0;
			
			//testing that the potential term is working for pot3
			if (ps.pot==3 && false) {
				comp Vtrial = 0.0, Vcontrol = 0.0;
				for (unsigned int j=0; j<N; j++) {
					double r = r0 + j*a;
					paramsV  = {r, 0.0};
					Vcontrol += pow(p(2*j*Nb),2.0)/2.0 - pow(p(2*j*Nb),4.0)/4.0/pow(r,2.0);
					Vtrial += V(p(2*j*Nb));
				}
				double potTest = pow(pow(real(Vcontrol-Vtrial),2.0) + pow(imag(Vcontrol-Vtrial),2.0),0.5);
				cout << "potTest = " << potTest << endl;
			}

/*----------------------------------------------------------------------------------------------------------------------------
	9. assigning minusDS, DDS etc
		- beginning loop over lattice points
----------------------------------------------------------------------------------------------------------------------------*/
			
			// beginning loop over lattice points
			for (unsigned long int j = 0; j < N*NT; j++)
				{		
				unsigned int t 			= intCoord(j,0,NT); //coordinates
				unsigned int x 			= intCoord(j,1,NT);
				int neighPosX 			= neigh(j,1,1,NT,N);
				int neighNegX  			= neigh(j,1,-1,NT,N);
				
				comp Dt 		= DtFn(t);
				comp dt 		= dtFn(t);
				comp dtm 		= (t>0? dtFn(t-1): b);
				double Dx 		= DxFn(x);
				double dx 		= dxFn(x);
				double dxm 		= (x>0? dxFn(x-1): a);
				
				if (ps.pot==3) paramsV  = {r0+x*a, 0.0};
			
				if (abs(chiX(j))>MIN_NUMBER && ps.pot!=3) //spatial zero mode lagrange constraint
					{
					DDS.insert(2*j,2*N*NT) 			= Dx*chiX(j); 
					DDS.insert(2*N*NT,2*j) 			= Dx*chiX(j);
					minusDS(2*j) 					+= -Dx*chiX(j)*p(2*N*NT);
					minusDS(2*N*NT) 				+= -Dx*chiX(j)*p(2*j);
					}
					
				if (abs(chiT(j))>MIN_NUMBER && t<(NT-1))
					{
					DDS.coeffRef(2*(j+1),2*N*NT+1) 	+= Dx*chiT(j); //chiT should be 0 at t=(NT-1) or this line will go wrong
					DDS.coeffRef(2*N*NT+1,2*(j+1)) 	+= Dx*chiT(j);
					DDS.coeffRef(2*j,2*N*NT+1) 		+= -Dx*chiT(j);
					DDS.coeffRef(2*N*NT+1,2*j) 		+= -Dx*chiT(j);
		            minusDS(2*(j+1)) 				+= -Dx*chiT(j)*p(2*N*NT+1);
		            minusDS(2*j) 					+= Dx*chiT(j)*p(2*N*NT+1);
		            minusDS(2*N*NT+1) 				+= -Dx*chiT(j)*(p(2*(j+1))-p(2*j));
					}
					
				if (t<(NT-1))
					{
					for (unsigned int k=0;k<N;k++)
						{
						unsigned int l = k*NT+t;
						linErgOffShell(t) += 0.5*( omega_2(x,k)*( Cp(l)-minima[0] )*( Cp(j)-minima[0] ) \
										+ omega_0(x,k)*( Cp(l+1)-Cp(l) )*( Cp(j+1)-Cp(j) )/pow(dt,2.0));
						linNumOffShell(t) += 0.5*(omega_1(x,k)*( Cp(l)-minima[0] )*( Cp(j)-minima[0] ) \
										+ omega_m1(x,k)*( Cp(l+1)-Cp(l) )*( Cp(j+1)-Cp(j) )/pow(dt,2.0));
						}
					}
					
				if (abs(theta)>MIN_NUMBER)
					{
					for (unsigned int k=0;k<N;k++)
						{
						unsigned int l = k*NT+t;
						linNum(t) += 2.0*Gamma*omega_1(x,k)*( (p(2*l)-minima[0])*(p(2*j)-minima[0])/pow(1.0+Gamma,2.0)\
								 + p(2*j+1)*p(2*l+1)/pow(1.0-Gamma,2.0) );
						linErg(t) += 2.0*Gamma*omega_2(x,k)*( (p(2*l)-minima[0])*(p(2*j)-minima[0])/pow(1.0+Gamma,2.0)\
								+ p(2*j+1)*p(2*l+1)/pow(1.0-Gamma,2.0) );
						}
					}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//boundaries
				if (ps.pot==3 && x==(N-1))
					{
					DDS.insert(2*j,2*j) 	= 1.0; // p=0 at r=R
					DDS.insert(2*j+1,2*j+1) = 1.0;
					}
				else if (ps.pot==3 && x==0)
					{
					kineticS 				+= 	Dt*pow(Cp(neighPosX),2.0)/dx/2.0;
					derivErg(t) 			+= 	pow(Cp(neighPosX),2.0)/dx/2.0; //n.b. the Dt/dt difference is ignored for erg(t)
					erg(t) 					+= 	pow(Cp(neighPosX),2.0)/dx/2.0;
					
					DDS.insert(2*j,2*j) 	= 1.0; // p=0 at r=0
					DDS.insert(2*j+1,2*j+1) = 1.0;
					}			
				else if (t==(NT-1))
					{
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
				else if (t==0)
					{
					kineticT 	+= Dx*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
					derivErg(t) += Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0;
					erg(t) 		+= Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0;
					
					if (bds.compare("jd")==0)
						{
						minusDS(2*j+1) 					+= Dx*imag(Cp(j+1)/dt);
					    DDS.coeffRef(2*j+1,2*(j+1)) 	+= -imag(Dx/dt);
					    DDS.coeffRef(2*j+1,2*(j+1)+1)	+= -real(Dx/dt);
					    minusDS(2*j+1) 					+= imag(-Dx*Cp(j)/dt);
					    DDS.coeffRef(2*j+1,2*j) 		+= imag(Dx/dt);
					    DDS.coeffRef(2*j+1,2*j+1) 		+= real(Dx/dt);
					    if (abs(theta)<MIN_NUMBER)
							{
							/////////////////////////////////////equation R - theta=0//////////////////////////////////////
							for (unsigned int k=0;k<N;k++)
								{
								if (abs(omega_1(x,k))>MIN_NUMBER)
									{
									unsigned int m=k*NT;
									DDS.coeffRef(2*j,2*m+1) += -2.0*omega_1(x,k);
									minusDS(2*j) 			+= 2.0*omega_1(x,k)*p(2*m+1);
									}
								}
							////////////////////////////////////////////////////////////////////////////////////////
							}
						else
							{
							for (unsigned int k=0;k<N;k++)
								{
								if (abs(omega_1(x,k))>MIN_NUMBER)
									{
									/////////////////////equation I - theta!=0//////////////
									unsigned int m=k*NT;
									DDS.coeffRef(2*j+1,2*m) += (1.0-Gamma)*omega_1(x,k)/(1.0+Gamma);
									minusDS(2*j+1) 			+= -(1.0-Gamma)*omega_1(x,k)*(p(2*m)-minima[0])/(1.0+Gamma);
									/////////////////////equation R - theta!=0//////////////
									minusDS(2*j) 			+= theta*p(2*m+1)*omega_1(x,k)*(1+Gamma)/(1-Gamma);
									DDS.coeffRef(2*j,2*m+1)	+= -theta*omega_1(x,k)*(1.0+Gamma)*theta/(1.0-Gamma);
									bound 		+= -(1.0-Gamma)*omega_1(x,k)*(p(2*j)-minima[0])*(p(2*m)-minima[0])/(1.0+Gamma)\
																 + (1.0+Gamma)*omega_1(x,k)*p(2*j+1)*p(2*m+1)/(1.0-Gamma);
									}
								}
							//////////////////////////////////////equation R - theta!=0//////////////////////////////
				            minusDS(2*j) 					+= theta*Dx*real(Cp(j+1)/dt);
			            	DDS.coeffRef(2*j,2*(j+1)) 		+= -theta*real(Dx/dt);
			            	DDS.coeffRef(2*j,2*(j+1)+1) 	+= theta*imag(Dx/dt);
						    minusDS(2*j) 					+= theta*real(-Dx*Cp(j)/dt);
					    	DDS.coeffRef(2*j,2*j) 			+= theta*real(Dx/dt);
					    	DDS.coeffRef(2*j,2*j+1) 		+= theta*imag(-Dx/dt);
							}
						}
					else ///////////////////////////////// including other terms in action at t=0 ///////////////////////////
						{
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
						for (unsigned int k=1; k<2*2; k++)
					    	{
					        int sign = pow(-1,k+1);
					        int direc = (int)(k/2.0);
					        int neighb = neigh(j,direc,sign,NT,N);
					        double dxd = (sign==1? dx: dxm);
					        if (direc == 0)
					        	{
					            minusDS(2*j+1) 					+= Dx*imag(Cp(j+sign)/dt);
					            DDS.coeffRef(2*j+1,2*(j+sign)) 	+= -imag(Dx/dt);
					            DDS.coeffRef(2*j+1,2*(j+sign)+1)+= -real(Dx/dt);
					            }
					        else if (neighb!=-1)
					        	{
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
						if (abs(theta)<MIN_NUMBER)
							{
							//simplest boundary conditions replaced by ones continuously connected to theta!=0 ones
							//DDS.insert(2*j+1,2*(j+1)+1) = 1.0; //zero imaginary part of time derivative
							//DDS.insert(2*j,2*j+1) = 1.0; //zero imaginary part
						
							/////////////////////////////////////equation R - theta=0//////////////////////////////////////
							for (unsigned int k=0;k<N;k++)
								{
								if (abs(omega_1(x,k))>MIN_NUMBER)
									{
									unsigned int m=k*NT;
									DDS.coeffRef(2*j,2*m+1) += -2.0*omega_1(x,k);
									minusDS(2*j) 			+= 2.0*omega_1(x,k)*p(2*m+1);
									}
								}
							////////////////////////////////////////////////////////////////////////////////////////
							}
						else
							{
							for (unsigned int k=0;k<N;k++)
								{
								if (abs(omega_1(x,k))>MIN_NUMBER)
									{
									/////////////////////equation I - theta!=0//////////////
									unsigned int m=k*NT;
									DDS.coeffRef(2*j+1,2*m) += (1.0-Gamma)*omega_1(x,k)/(1.0+Gamma);
									minusDS(2*j+1) 			+= -(1.0-Gamma)*omega_1(x,k)*(p(2*m)-minima[0])/(1.0+Gamma);
									/////////////////////equation R - theta!=0//////////////
									minusDS(2*j) 			+= p(2*m+1)*omega_1(x,k)*(1+Gamma)*theta/(1-Gamma);
									DDS.coeffRef(2*j,2*m+1)	+= -omega_1(x,k)*(1.0+Gamma)*theta/(1.0-Gamma);
									bound 				+= -(1.0-Gamma)*omega_1(x,k)*(p(2*j)-minima[0])*(p(2*m)-minima[0])/(1.0+Gamma)\
																 + (1.0+Gamma)*omega_1(x,k)*p(2*j+1)*p(2*m+1)/(1.0-Gamma);
									}
								}
							//////////////////////////////////////equation R - theta!=0//////////////////////////////
							for (unsigned int k=1; k<2*2; k++)
						    	{
						        int sign = pow(-1,k+1);
						        int direc = (int)(k/2.0);
						        int neighb = neigh(j,direc,sign,NT,N);
						        double dxd = (sign==1? dx: dxm);
						        if (direc == 0)
						        	{
						            minusDS(2*j) 					+= Dx*real(Cp(j+sign)/dt)*theta;
					            	DDS.coeffRef(2*j,2*(j+sign)) 	+= -real(Dx/dt)*theta;
					            	DDS.coeffRef(2*j,2*(j+sign)+1) 	+= imag(Dx/dt)*theta;
						            }
						        else if (neighb!=-1)
						        	{
						            minusDS(2*j) 					+= - real(Dt*Cp(neighb))*theta/dxd;
						            DDS.coeffRef(2*j,2*neighb) 		+= real(Dt)*theta/dxd;
					            	DDS.coeffRef(2*j,2*neighb+1) 	+= -imag(Dt)*theta/dxd;
						            }
						        }
						    minusDS(2*j) 			+= theta*real(-temp0*Cp(j) + temp1 );
					    	DDS.coeffRef(2*j,2*j) 	+= theta*real(temp0 - temp2 );
					    	DDS.coeffRef(2*j,2*j+1) += theta*imag(-temp0 + temp2 );
							}
						}
					}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//bulk
				else
					{
					if (neighPosX!=-1) {
						kineticS 	+= Dt*pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
						erg(t) 	 	+= pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
						derivErg(t) += pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
					}
					
					kineticT 	+= Dx*pow(Cp(j+1)-Cp(j),2.0)/dt/2.0;
					potV 		+= Dt*Dx*V(Cp(j));
					pot_r 		+= Dt*Dx*Vr(Cp(j));
<<<<<<< HEAD
					erg(t) 		+= Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + Dx*V(Cp(j)) + Dx*Vr(Cp(j));
					derivErg(t) += Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0;
					potErg(t) 	+= Dx*V(Cp(j)) + Dx*Vr(Cp(j));
=======
					erg(t) 		+= Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0 + Dx*V(Cp(j)) + Dx*Vr(Cp(j));
					derivErg(t) += Dx*pow(Cp(j+1)-Cp(j),2.0)/pow(dt,2.0)/2.0 + pow(Cp(neighPosX)-Cp(j),2.0)/dx/2.0;
					potErg(t) 	+= Dx*(V(Cp(j)) + Vr(Cp(j)));
>>>>>>> 3502b10cf58fbb10751d0562ef85a0f133e7ed72
				
		            for (unsigned int k=0; k<2*2; k++)
                	{
                    int sign = pow(-1,k);
                    int direc = (int)(k/2.0);
                    int neighb = neigh(j,direc,sign,NT,N);
                    comp dtd = (sign==1? dt: dtm);
                    double dxd = (sign==1? dx: dxm);
                    if (direc == 0)
                    	{
                        minusDS(2*j) 					+= real(Dx*Cp(j+sign)/dtd);
                        minusDS(2*j+1) 					+= imag(Dx*Cp(j+sign)/dtd);
                        DDS.insert(2*j,2*(j+sign)) 		= -real(Dx/dtd);
                        DDS.insert(2*j,2*(j+sign)+1) 	= imag(Dx/dtd);
                        DDS.insert(2*j+1,2*(j+sign)) 	= -imag(Dx/dtd);
                        DDS.insert(2*j+1,2*(j+sign)+1) 	= -real(Dx/dtd);
                        }
                    else if (neighb!=-1)
                    	{                        
                        minusDS(2*j) 					+= -real(Dt*Cp(neighb)/dxd);
                        minusDS(2*j+1) 					+= -imag(Dt*Cp(neighb)/dxd);
                        DDS.insert(2*j,2*neighb) 		= real(Dt/dxd);
                        DDS.insert(2*j,2*neighb+1) 		= -imag(Dt/dxd);
                        DDS.insert(2*j+1,2*neighb) 		= imag(Dt/dxd);
                        DDS.insert(2*j+1,2*neighb+1) 	= real(Dt/dxd);
                        }
                    }
		            comp temp0 = Dx*(1.0/dt + 1.0/dtm) - Dt*(1.0/dx + 1.0/dxm);
		            if (neighPosX==-1) temp0 += Dt/dx;
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
		        }
		    if (ps.pot==3) DDS.insert(2*N*NT,2*N*NT) = 1.0;
		    action = kineticT - kineticS - potV - pot_r;
		    linErgOffShell(NT-1) = linErgOffShell(NT-2);
	    	linNumOffShell(NT-1) = linNumOffShell(NT-2);
	    	
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
		   
	    	
		    if (abs(theta)<MIN_NUMBER) {
		    	linErg = linErgOffShell;
		    	linNum = linNumOffShell;
		    }
		    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
			//checking linErg, linNum, bound, computing W, E
			
			//trivial test that erg=potErg+derivErg
			if (false) {
				for (unsigned int j=0; j<NT; j++) {
					double diff = absDiff(erg(j),potErg(j)+derivErg(j));
					if (diff>1.0e-14) cerr << "erg(" << j << ") != potErg + derivErg. absDiff = " << diff << endl;
				}
			}
			
			//calculating continuum approx to linErg and linNum on initial time slice
			if (ps.pot==3) {
				for (unsigned int k=1; k<N; k++) {
					double momtm = k*pi/(L-r0);
					double freqSqrd = 1.0+pow(momtm,2.0);
					double Asqrd, integral1 = 0.0, integral2 = 0.0;
					for (unsigned int l=0; l<N; l++) {
						double r = r0 + l*a;
						unsigned int m = l*NT;
						integral1 += a*p(2*m)*pow(2.0/L,0.5)*sin(momtm*r);
						integral2 += a*(p(2*(m+1))-p(2*m))*pow(2.0/L,0.5)*sin(momtm*r)/b;
					}
					Asqrd = pow(integral1,2.0) + pow(integral2,2.0)/freqSqrd;
					linErgContm += 2.0*pi*Asqrd*freqSqrd;
					linNumContm += 2.0*pi*Asqrd*pow(freqSqrd,0.5);
				}
			}
			
			//calculating a_k and b*_k
			cVec a_k(N), b_k(N); //N.b. b_k means b*_k but can't put * in the name
			a_k = Eigen::VectorXcd::Zero(N);
			b_k = Eigen::VectorXcd::Zero(N);
			double T0 = real(coord(0,0))+Ta;
			comp dt0 = dtFn(0);
			for (unsigned int n=0; n<N; n++) {
				double w_n = sqrt(real(eigenValues(n)));
				double w_n_e = w_n_exp(n);
				for (unsigned int j=0; j<N; j++) {
					unsigned int m=j*NT;
					double sqrtDj = sqrt(4.0*pi*DxFn(j));
					if (abs(w_n)>1.0e-16 && abs(w_n_e)>1.0e-16) {
						a_k(n) += exp(ii*w_n_e*T0)*sqrt(2.0*w_n)*eigenVectors(j,n)* \
									sqrtDj*((Cp(m+1)-minima[0])-(Cp(m)-minima[0])*exp(ii*w_n_e*dt0)) \
										/(exp(-ii*w_n_e*dt0)-exp(ii*w_n_e*dt0));
						b_k(n) += exp(-ii*w_n_e*T0)*sqrt(2.0*w_n)*eigenVectors(j,n)* \
									sqrtDj*((Cp(m+1)-minima[0])-(Cp(m)-minima[0])*exp(-ii*w_n_e*dt0)) \
										/(exp(ii*w_n_e*dt0)-exp(-ii*w_n_e*dt0));
					}
				}
			}			
			
			//using a_k and b*_k to check that a_k=Gamma*b*_k and p->p_lin as t->0 and that Sum_k(a_k*b*_k)=linNum
			// and that Sum_k(w_k*a_k*b*_k)=linErg
			double ABtest = 0.0, linRepTest, ABNEtest;
			linNumAB = 0.0, linErgAB = 0.0;
			cVec linRep(N), p0(N);
			linRep = Eigen::VectorXcd::Constant(N,minima[0]);
			for (unsigned int j=0; j<N; j++) {
				ABtest += absDiff(a_k(j),conj(b_k(j))*Gamma);
				p0(j) = Cp(j*NT);
				if (j<N) {
					linNumAB += a_k(j)*b_k(j);
					linErgAB += sqrt(eigenValues(j))*a_k(j)*b_k(j);
				}
				for (unsigned int n=0; n<N; n++) {
					double w_n = sqrt(real(eigenValues(n)));
					double w_n_e = w_n_exp(n);
					double sqrtDj = sqrt(4.0*pi*DxFn(j));
					if (abs(w_n)>1.0e-16) {
						linRep(j) += eigenVectors(j,n)*(a_k(n)*exp(-ii*w_n_e*T0)+b_k(n)*exp(ii*w_n_e*T0)) \
										/sqrt(2.0*w_n)/sqrtDj;
					}
				}
			}
			ABtest /= (double)N;
			linRepTest = absDiff(p0,linRep);
			double ABNtest = absDiff(linNumAB,linNum(0));
			double ABEtest = absDiff(linErgAB,linErg(0));
			ABNEtest = (ABNtest>ABEtest? ABNtest: ABEtest);
			AB_test.push_back(ABtest);
			linRep_test.push_back(linRepTest);
			ABNE_test.push_back(ABNEtest);
			/*if (AB_test.back()>closenessAB || !isfinite(AB_test.back())) cerr << "AB_test = " << AB_test.back() << endl;
			if (linRep_test.back()>closenessLR || !isfinite(linRep_test.back()))cerr << "linRep_test = " << linRep_test.back() << endl;
			if (ABNE_test.back()>closenessABNE || !isfinite(ABNE_test.back())) {
				cerr << "linNumAB = " << linNumAB << ", linNum(0) = " << \
					linNum(0) << ", linNumOffShell(0) = " << linNumOffShell(0) << endl;
				cerr << "linErg = " << linErgAB << ", linErg(0) = " << linErg(0) \
					<< ", linErgOffShell(0) = " << linErgOffShell(0) << endl;
			}*/
			
			//checking pot_r is much smaller than the other potential terms
			reg_test.push_back(abs(pot_r/potV));
			if (reg_test.back()>closenessR)
				{
				cerr << "regularisation term is too large, regTest = " << reg_test.back() << endl;
				}
			
			//checking linearisation of linErg and linNum
			double linTestE;	double linTestN;
			double linEMax = 0.0;	double linEMin = 5.0e15; //surely it's going to be less than this
			double linNMax = 0.0;	double linNMin = 5.0e15;
			unsigned int linearInt = (unsigned int)(Na/10);
			for (unsigned int j=1;j<(linearInt+1);j++)
				{
				if (abs(linErgOffShell(j))>linEMax) linEMax = abs(linErgOffShell(j));
				if (abs(linErgOffShell(j))<linEMin) linEMin = abs(linErgOffShell(j));
				if (abs(linNumOffShell(j))>linNMax) linNMax = abs(linNumOffShell(j));
				if (abs(linNumOffShell(j))<linNMin) linNMin = abs(linNumOffShell(j));
				}
			linTestE = absDiff(linEMax,linEMin);
			linTestN = absDiff(linNMax,linNMin);
			lin_test.push_back(linTestE);
			if (linTestN>linTestE) lin_test.back() = linTestN;
			
			//checking agreement of on and off shell linear energy at initial time
			double onShellTest = absDiff(linErg(0),linErgOffShell(0));
			onShell_test.push_back(onShellTest);
			
			//checking conservation of E
			double conservTest = absDiff(erg(1),erg(NT-2));
			conserv_test.push_back(conservTest);
						
			//defining E, Num and cW
			E = real(linErg(0));
			Num = real(linNum(0));
			E_exact = 0.0;
			for (unsigned int j=1; j<(linearInt+1); j++) E_exact += real(erg(j));
			E_exact /= (double)linearInt;
			W = - E*2.0*Tb - theta*Num - bound + 2.0*imag(action);
			
			//testing imaginary part of energy
			double imErgTest = 0.0;
			for (unsigned int j=0; j<NT; j++) if (abs(erg(j))>MIN_NUMBER) imErgTest += imag(erg(j))/abs(erg(j));
			imErgTest /= (double)NT;
			imErg_test.push_back(imErgTest);
			
			//checking agreement between erg and linErg
			double trueTest = absDiff(E,E_exact);
			true_test.push_back(trueTest);
			if (!isfinite(trueTest)) cout << "E = " << E << ", E_exact = " << E_exact << ", linearInt = " << linearInt << endl;
			
			//checking lattice small enough for E, should have parameter for this
			double momTest = E*b/Num/pi; //perhaps should have a not b here
			mom_test.push_back(momTest);
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			//solving for delta in DDS*delta=minusDS, where p' = p + delta		
			vec delta(2*N*NT+2);
			delta = Eigen::VectorXd::Zero(2*N*NT+2);
			DDS.prune(MIN_NUMBER);
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
			vec diff(2*N*NT+2);
			diff = DDS*delta-minusDS;
			double maxDiff = diff.maxCoeff();
			maxDiff = abs(maxDiff);
			calc_test.push_back(maxDiff);
			if (calc_test.back()>closenessC)
				{
				cerr << "Calculation failed" << endl;
				cerr << "calc_test = " << calc_test.back() << endl;
				return 1;
				}

			//assigning values to phi
			p += delta;
		
			//passing changes on to complex vector
			Cp = vecComplex(p,N*NT);
			
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		//printing early if desired
		if (runs_count == aq.printRun || aq.printRun == 0)
			{
			string prefix = "./data/" + timeNumber;
			string suffix = "_" + numberToString<unsigned int>(fileLoop) + "_" + numberToString<unsigned int>(loop)+"_"+numberToString<unsigned int>(runs_count)+".dat";
			if ((print_choice.compare("v")==0 || print_choice.compare("e")==0))
				{
				string minusDSfile = prefix + "mainminusDSE"+suffix;
				printVector(minusDSfile,minusDS);
				}
			if ((print_choice.compare("p")==0 || print_choice.compare("e")==0) || delta_test.back()>0.2)
				{
				string piEarlyFile = prefix + "mainpiE"+suffix;
				printVector(piEarlyFile,p);
				//gp(piEarlyFile,"repi.gp");
				}
			if ((print_choice.compare("m")==0 || print_choice.compare("e")==0))
				{
				string DDSfile = prefix + "mainDDSE"+suffix;
				printSpmat(DDSfile,DDS);
				}
			if ((print_choice.compare("z")==0 || print_choice.compare("e")==0))
				{
				string earlychiXFile = prefix + "mainchiXE" + suffix;
				printVector(earlychiXFile,chiX);
				//gp(earlychiXFile,"repi.gp");
				string earlychiTFile = prefix + "mainchiTE" + suffix;
				printVector(earlychiTFile,chiT);
				//gp(earlychiTFile,"repi.gp");
				}
			if ((print_choice.compare("d")==0 || print_choice.compare("e")==0))
				{
				string earlydeltaFile = prefix + "maindeltaE" + suffix;
				printVector(earlydeltaFile,delta);
				}
			if ((print_choice.compare("l")==0 || print_choice.compare("e")==0))
				{
				string earlyLinErgFile = prefix + "mainlinErgE"+suffix;
				simplePrintCVector(earlyLinErgFile,linErg);
				//gpSimple(earlyLinErgFile);
				string earlyErgFile = prefix + "mainergE" + suffix;
				simplePrintCVector(earlyErgFile,erg);
				//gpSimple(earlyErgFile);
				string earlyDerivErgFile = prefix + "mainderivErgE"+suffix;
				//derivErg.conservativeResize(Na);
				simplePrintCVector(earlyDerivErgFile,derivErg);
				string earlyPotErgFile = prefix + "mainpotErgE"+suffix;
				//potErg.conservativeResize(Na);
				simplePrintCVector(earlyPotErgFile,potErg);
				}
			if ((print_choice.compare("ab")==0 || print_choice.compare("e")==0))
				{
				string earlya_kFile = prefix + "maina_kE"+suffix;
				simplePrintCVector(earlya_kFile,a_k);
				string earlyb_kFile = prefix + "mainb_kE" + suffix;
				simplePrintCVector(earlyb_kFile,b_k);
				}
			}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
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
		
			//assigning test values
			//quantities used to stop newton-raphson loop
			action_test.push_back(absDiff(action,action_last));
			action_last = action;
			sol_test.push_back(normDS/normP);
			solM_test.push_back(maxDS/maxP);
			delta_test.push_back(normDelta/normP);
			
			//printing tests to see convergence
			if (runs_count==1)
				{
				printf("%5s%5s%11s%11s%11s%11s%11s%11s%11s%11s%11s%11s\n","loop","run","solTest","solMTest","deltaTest","linTest","trueTest","shellTest","ABTest","ABNETest","cnsrvTest","momTest");
				}
			printf("%5i%5i%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g%11.4g\n",loop,runs_count,sol_test.back(),solM_test.back(),delta_test.back(),lin_test.back(),true_test.back(),onShell_test.back(),AB_test.back(),ABNE_test.back(),conserv_test.back(),mom_test.back());
			
			if (delta_test.back()>closenessD) {
				cerr << "deltaTest > " << closenessD << ", stopping NR loop." << endl;
				break;
			}
			
			} //ending while loop
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	    		//misc end of program tasks - mostly printing

		//checking different measures of energy agree
		if (true_test.back()>closenessT || onShell_test.back()>closenessON || ABNE_test.back()>closenessABNE || true) {
			cout << "erg(1)             = " << erg(1) << endl;
			cout << "linErgOffShell(0)  = " << linErgOffShell(0) << endl;
			cout << "linErgAB           = " << linErgAB << endl;
			cout << "linErgContm        = " << linErgContm << endl;
			cout << "linErg(0)          = " << linErg(0) << endl;
		}

		//checking energy conserved
		if (conserv_test.back()>closenessCon) cout << endl << "ergTest = " << conserv_test.back() << endl;
		
		//checking energy real
		if (imErg_test.back()>closenessIE) cout << endl << "imErg = "<< imErg_test.back()  << endl;
		
		//checking lattice small enough
		if (mom_test.back()>closenessP) cout << endl << "momTest = "<< mom_test.back()  << endl;
		
		if (ps.pot==3 && abs(theta)<MIN_NUMBER) {
			comp temp = linErg(0);
			if (absDiff(temp,linErgContm)>closenessCL) 
				cout << "linErg(0) = " << temp << " and linErgContm = " << linErgContm << " don't agree. E_exact = " << E_exact << endl;
			temp = linNum(0);
			if (absDiff(temp,linNumContm)>closenessCL) 
				cout << "linNum(0) = " << temp << " and linNumContm = " << linNumContm << " don't agree. " << endl;
		}
		
		//stopping clock
		time = clock() - time;
		double realtime = time/1000000.0;
	
		//printing to terminal
		printf("\n");
		printf("%8s%8s%8s%8s%8s%8s%8s%8s%14s%14s%14s%14s\n","runs","time","N","NT","L","Tb","dE","theta","Num","E","im(action)","W");
		printf("%8i%8g%8i%8i%8g%8g%8g%8g%14.4g%14.4g%14.4g%14.4g\n",runs_count,realtime,N,NT,L,Tb,dE,theta,Num,E,imag(action),real(W));
		printf("\n");
		 printf("%60s\n","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

		//printing action value
		FILE * actionfile;
		actionfile = fopen("./data/mainAction.dat","a");
		fprintf(actionfile,"%12s%6i%6i%8g%8g%8g%8g%10.4g%10.4g%10.4g%10.4g%10.4g%10.4g%10.4g\n",timeNumber.c_str(),N,NT,L,Tb,dE,theta,E,Num,imag(action)\
		,real(W),sol_test.back(),lin_test.back(),true_test.back());
		fclose(actionfile);
		
		string prefix = "./data/" + timeNumber;
		string suffix = "_" + numberToString<unsigned int>(fileLoop)+"_" + numberToString<unsigned int>(loop)+ ".dat";
	
		//copying a version of inputs with timeNumber and theta changed
		string runInputs = prefix + "inputsM_"+ numberToString<unsigned int>(fileLoop) + "_" + numberToString<unsigned int>(loop); //different suffix
		if (abs(maxTheta-minTheta)>MIN_NUMBER) 	changeInputs(runInputs,"theta", numberToString<double>(theta),inputsFiles[fileLoop]);
		else if (abs(maxTb-minTb)>MIN_NUMBER)	changeInputs(runInputs,"Tb", numberToString<double>(Tb),inputsFiles[fileLoop]);
		else 									copyFile(inputsFiles[fileLoop],runInputs);
		printf("%12s%30s\n","output: ",runInputs.c_str());
	
		//printing output phi
		string tpifile =  prefix + "mainpi"+suffix;
		printVector(tpifile,p);
		//gp(tpifile,"repi.gp");
		printf("%12s%30s\n"," ",tpifile.c_str());
	
		//printing output minusDS				
		string minusDSfile = prefix + "mainminusDS"+suffix;
		printVector(minusDSfile,minusDS);
		printf("%12s%30s\n"," ",minusDSfile.c_str());
				
		//printing output DDS
		string DDSfile = prefix + "mainDDS"+suffix;
		printSpmat(DDSfile,DDS);
	    printf("%12s%30s\n"," ",DDSfile.c_str());
	    
		//printing linNum
		string linNumFile = prefix + "mainlinNum"+suffix;
		linNum.conservativeResize(Na);
		simplePrintCVector(linNumFile,linNum);
		printf("%12s%30s\n"," ",linNumFile.c_str());
		//gpSimple(linNumFile);
	
		//printing linErg
		string linErgFile = prefix + "mainlinErg"+suffix;
		linErg.conservativeResize(Na);
		simplePrintCVector(linErgFile,linErg);
		//gpSimple(linErgFile);
		printf("%12s%30s\n"," ",linErgFile.c_str());
		
		//printing linErgOffShell
		string linErgOffShellFile = prefix + "mainlinErgOffShell"+suffix;
		linErgOffShell.conservativeResize(Na);
		simplePrintCVector(linErgOffShellFile,linErgOffShell);
		//gpSimple(linErgFileOffShell);
		printf("%12s%30s\n"," ",linErgOffShellFile.c_str());
		
		//printing derivErg
		string derivErgFile = prefix + "mainderivErg"+suffix;
		//derivErg.conservativeResize(Na);
		simplePrintCVector(derivErgFile,derivErg);
		//gpSimple(linErgFile);
		printf("%12s%30s\n"," ",derivErgFile.c_str());
		
		//printing potErg
		string potErgFile = prefix + "mainpotErg"+suffix;
		//potErg.conservativeResize(Na);
		simplePrintCVector(potErgFile,potErg);
		//gpSimple(linErgFile);
		printf("%12s%30s\n"," ",potErgFile.c_str());
	
		//printing erg
		string ergFile = prefix + "mainerg" + suffix;
		erg.conservativeResize(Na);
		simplePrintCVector(ergFile,erg);
		//gpSimple(ergFile);
		printf("%12s%30s\n"," ",ergFile.c_str());
		cout << endl;
		
		if (delta_test.back()>closenessD) {
				return 1;
			}
		
		} //ending parameter loop
	} //ending file loop

return 0;
}


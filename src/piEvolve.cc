/*---------------------------------------------------------------------------------------------
	piEvolve
		- something to compile results to get periodic instanton for pot=3
---------------------------------------------------------------------------------------------*/
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <string>
#include "simple.h"
#include "folder.h"
#include "parameters.h"
#include "print.h"


//#define NDEBUG //NDEBUG is to remove error and bounds checking on vectors in SparseLU, for speed - only include once everything works

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - defining key parameters
		2 - argv inputs
		3 - getting parameters from specific inputs
		4 - loading phi on BC
		5 - omega
		6 - propagating euclidean solution along B->A and C->D
		7 - compiling and printing results
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/


int main(int argc, char** argv) {

/* ---------------------------------------------------------------------------------------------
	1. defining key parameters
		- parameters
		- specific options
---------------------------------------------------------------------------------------------*/

// parameters
Parameters ps_run;
Parameters ps_in;

// specific options
int direction = 1; // direction of time evolution
double sigma = 1.0; // set sigma=-1 for euclidean evolution; sigma=1 for minkowskian
bool testTunnel = false, testLinear = false, changeParams = false, approxOmega = false;
double closenessLin = 1.0e-2, Tlin = 0.0;
double closenessEdge = 1.0e-3, Redge = 0.0;
double closenessMom = 1.9e-2, momTest;

/* ---------------------------------------------------------------------------------------------
	2. argv inputs
		- getting argv inputs
		- printing type of test if chosen 
---------------------------------------------------------------------------------------------*/

// getting argv inputs
string timenumber;
string timenumberIn = "150112114306";
string loopIn = "0";
if (argc == 2) timenumberIn = argv[1];
else if (argc % 2 && argc>1) {
	for (unsigned int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("tn")==0) timenumberIn = argv[2*j+2];
		else if (id.compare("r1")==0) timenumberIn = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("loop")==0 || id.compare("l")==0) loopIn = argv[2*j+2];
		else if (id.compare("test")==0 || id.compare("tunnel")==0 || id.compare("tt")==0) testTunnel = (bool)atoi(argv[2*j+2]);
		else if (id.compare("linearization")==0 || id.compare("lin")==0) testLinear = (bool)atoi(argv[2*j+2]);
		else if (id.compare("closeness")==0 || id.compare("close")==0) closenessLin = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("changeParams")==0) changeParams = (bool)atoi(argv[2*j+2]);
		else if (id.compare("approxOmega")==0) approxOmega = (bool)atoi(argv[2*j+2]);
		else if (id.compare("N")==0) N = atoi(argv[2*j+2]);
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}
else if (argc == 1) {
	cerr << "must provide input to pi3.cc" << endl;
	return 1;
}
else {
	cerr << "must provide an even number of inputs in format '-name value':" << endl;
	for (int j=1; j<argc; j++) cerr << argv[j] << " ";
	cerr << endl;
	return 1;
}


// printing type of test if chosen 
timenumber = timenumberIn;
if (testTunnel) {
	cout << "pi3 testing if tunnelled" << endl;
}
if (testLinear) {
	testTunnel = false;
	cout << "pi3 testing linearization" << endl;
}

/* ---------------------------------------------------------------------------------------------
	3. getting parameters from specific inputs
		- getting inputs corresponding to timenumberIn
		- chosing p_run
		- checking parameters
---------------------------------------------------------------------------------------------*/

// getting inputs corresponding to timenumberIn
string prefix = "./data/"+timenumberIn;
string suffix = "_loop_"+numberToString<uint>(loopIn)+".dat";
Filename inputsFile = (string)(prefix+"inputsP"+suffix);
inputsFile.Suffix = "";
ps_in.load(inputsFile);

// chosing p_run
ps_run = ps_in;
ps_run.changeParameters("N",300);
	{
	double dt = ps_run.a/5.0;
	ps_run.changeParameters("Na",(unsigned int)(ps_in.Ta/dt));
	ps_run.changeParameters("Nb",(unsigned int)(ps_in.Tb/dt));
	ps_run.changeParameters("Nc",(unsigned int)(ps_in.Tc/dt));
	}

// checking parameters
if (abs(ps_in.Ta)>1.1*ps_in.L && !testTunnel && !testLinear) {
	cerr << "R is too small compared to Ta. R = " << ps_in.L << ", Ta = " << ps_in.Ta << endl;
	return 1;
}
if (abs(ps_in.Tc)>1.1*ps_in.L && !testTunnel && !testLinear) {
	cerr << "R is too small compared to Tc. R = " << ps_in.L << ", Tc = " << ps_in.Tc << endl;
	return 1;
}
if (abs(ps_run.b)>0.5*ps_run.a) {
	cerr << "dt too large. dt = " << ps_run.b << ", dr = " << ps_run.a << endl;
	return 1;
}
if (abs(ps_in.Ta)<2.5 && !testTunnel && !testLinear) {
	cerr << "Ta too small. Ta = " << ps_in.Ta << endl;
	return 1;
}

/* ---------------------------------------------------------------------------------------------
	4. loading phi on BC
		- load options
		- defining phiBC
		- loading phiBC
---------------------------------------------------------------------------------------------*/

// load options
Filename pFile = (string)(prefix+"p"+suffix);
SaveOptions so_in;
so_in.paramsIn = ps_in;
so_in.paramsOut = ps_in;
so_in.vectorType = SaveOptions::complexB;
so_in.extras = SaveOptions::coords;
so_in.zeroModes = 1;
so_in.printMessage = true;

// defining phiBC
vec phiBC;

// loading phiBC
load(pFile,so_p,phiBC);

/* ---------------------------------------------------------------------------------------------
	5. getting omega matrices
		- loading if possible
		- otherwise constructing omegas and saving results
---------------------------------------------------------------------------------------------*/

//deterimining omega matrices for fourier transforms in spatial direction
	vec freqs(ps_run.N), freqs_exp(ps_run.N);
	mat modes(ps_run.N,ps_run.N);
	mat omega_m1(ps_run.N,ps_run.N), omega_0(ps_run.N,ps_run.N), omega_1(ps_run.N,ps_run.N), omega_2(ps_run.N,ps_run.N);
	SaveOptions so_simple;
	so_simple.paramsIn = ps_run; so_simple.paramsOut = ps_run;
	so_simple.vectorType = SaveOptions::simple;
	so_simple.extras = SaveOptions::none;
	so_simple.printMessage = false;
	{
		Filename omegaM1F, omega0F, omega1F, omega2F, modesF, freqsF, freqsExpF; // Filename works as FilenameAttributes
		omegaM1F = (string)("data/stable/omegaM1_pot_"+numberToString<uint>(ps_run.pot)+"_N_"+numberToString<uint>(ps_run.N)\
						+"_L_"+numberToString<double>(ps_run.L)+".dat");
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
				numericalModes(modes,freqs,freqs_exp,ps_run);
			}
			else {
				analyticModes(modes,freqs,freqs_exp,ps_run);
			}
			omegasFn(approxOmega,modes,freqs,omega_m1,omega_0,omega_1,omega_2,ps_run);
			save(omegaM1F,so_simple,omega_m1);
			save(omega0F,so_simple,omega_0);
			save(omega1F,so_simple,omega_1);
			save(omega2F,so_simple,omega_2);
			save(modesF,so_simple,modes);
			save(freqsF,so_simple,freqs);
			save(freqsExpF,so_simple,freqs_exp);
		}
	}

/* ---------------------------------------------------------------------------------------------
	6. propagating euclidean solution along B->A and C->D
---------------------------------------------------------------------------------------------*/

vec phiA, phiC, linearizationA;
double linErgContm, linNumContm, nonLinErgA, linErgFieldA, ergA, linErgA, linNumA;
double kineticT, kineticS, massTerm;
double Tinf = 0.0;

uint j=0;
while(j<2) {

	if (testLinear) {
		Nt = 10*ps_run.N;
		j=1;
		direction = -1;
	}
	else if (testTunnel) {
		Nt = 10*ps_run.N;
		j=1;
		direction = 1;
	}
	else if (j==0) {
		Nt = ps_run.Nc;
		direction = 1;
	}
	else if (j==1) {
		Nt = ps_run.Na;
		direction = -1;
	}
	
	vec phi((Nt+1)*ps_run.N), vel((Nt+1)*ps_run.N), acc((Nt+1)*ps_run.N);
	vec nonLinErg(Nt+1), linErgField(Nt+1), erg(Nt+1);
	linErgField = Eigen::VectorXd::Zero(Nt+1);
	nonLinErg = Eigen::VectorXd::Zero(Nt+1);
	erg = Eigen::VectorXd::Zero(Nt+1);
	linErgContm = 0.0, linNumContm = 0.0;
	linErgA = 0.0, linNumA = 0.0;
	kineticT = 0.0, kineticS = 0.0, massTerm = 0.0;
	vec initial;
	initial = interpolate(phiBC,ps_in.Nb,ps_in.N,Nt+1,ps_run.N);

	//initial condition 1)
	//defining phi on initial time slice
	for (unsigned int j=0;j<ps_run.N;j++) {
		unsigned int l = j*(Nt+1);
		unsigned int m;
		(direction == 1) ? m = Nt + l : m = l;
		double r = r0 + j*dr;
		phi(l) = initial(m)/r;
	}
	
	//initialising linErg and linNum	
	for (unsigned int x=0;x<ps_run.N;x++) {
		unsigned int j = x*(Nt+1);
		double r = r0 + x*dr, eta;
		eta = ((x==0 || x==(ps_run.N-1)) ? 0.5 : 1.0);
		nonLinErg(0) += 4.0*pi*pow(r,2.0)*eta*0.25*pow(phi(j),4.0)*dr;
		linErgField(0) += 4.0*pi*pow(r,2.0)*eta*0.5*pow(phi(j),2.0)*dr;
		if (x<(ps_run.N-1)) linErgField(0) += 4.0*pi*r*(r+dr)*0.5*pow(phi(j+(Nt+1))-phi(j),2.0)/dr;
	}

	//initial condition 2)
	//intitialize velocity, dphi/dt(x,0) = 0;
	for (unsigned int j=0; j<ps_run.N; j++) {
		unsigned int l = j*(Nt+1);
		vel(l) = 0.0;
	}

	//boundary condition 1)
	//phi=0.0 at r=L
	for (unsigned int j=0;j<(Nt+1);j++) {
		unsigned int l = (ps_run.N-1)*(Nt+1) + j;
		phi(l) = 0.0;
	}

	//boundary condition 2)
	//initialize acc using phi and expression from equation of motion
	/*the unusual d^2phi/dr^2 term and the absence of the first derivative term
	are due to the boundary condition 2) which, to second order, is phi(t,-dr) = phi(t,dr)*/
	acc(0) = 2.0*(phi(Nt+1) - phi(0))/pow(dr,2.0) - phi(0) + pow(phi(0),3.0);
	acc(0) *= sigma*0.5; //as initial time slice, generated from taylor expansion and equation of motion
	for (unsigned int j=1; j<(ps_run.N-1); j++) {
		unsigned int l = j*(Nt+1);
		double r = r0 + dr*j;
		acc(l) = (phi(l+(Nt+1)) + phi(l-(Nt+1)) - 2.0*phi(l))/pow(dr,2.0) + (phi(l+(Nt+1))-phi(l-(Nt+1)))/r/dr - phi(l) + pow(phi(l),3.0);
		acc(l) *= sigma*0.5;
	}

	//A7. run loop
	for (unsigned int u=1; u<(Nt+1); u++) {
		for (unsigned int x=0; x<(ps_run.N-1); x++) { //don't loop over last x position as fixed by boundary condition 1)
		    unsigned int m = u+x*(Nt+1);
		    vel(m) = vel(m-1) + dt*acc(m-1);
		    phi(m) = phi(m-1) + dt*vel(m);
		    if (testTunnel) {
		    	double testInf = phi(m);
		    	if (!isfinite(testInf) && abs(Tinf)<1.0e-16) Tinf = u*dt;
		    }
		}
		acc(u) = 2.0*(phi(u+Nt+1) - phi(u))/pow(dr,2.0) - phi(u) + pow(phi(u),3.0);
		acc(u) *= sigma;
		linErgField(u-1) += 4.0*pi*dr*pow(r0,2.0)*0.5*pow(phi(u)-phi(u-1),2.0)/pow(dt,2.0);
		linErgField(u-1) += 4.0*pi*dr*pow(r1,2.0)*0.5*pow(phi(u+(ps_run.N-1)*(Nt+1))-phi(u+(ps_run.N-1)*(Nt+1)-1),2.0)/pow(dt,2.0);
		linErgField(u) += 4.0*pi*( r0*(r0+dr)*0.5*pow(phi(u+(Nt+1))-phi(u),2.0)/dr + dr*pow(r0,2.0)*0.5*pow(phi(u),2.0) );
		linErgField(u) += 4.0*pi*dr*pow(r1,2.0)*0.5*pow(phi(u+(ps_run.N-1)*(Nt+1)),2.0);
		nonLinErg(u) += 4.0*pi*dr*pow(r0,2.0)*0.25*pow(phi(u),4.0);
		nonLinErg(u) += 4.0*pi*dr*pow(r1,2.0)*0.25*pow(phi(u+(ps_run.N-1)*(Nt+1)),4.0);
		for (unsigned int x=1; x<(ps_run.N-1); x++) {
		    unsigned int m = u+x*(Nt+1);
		    double r = r0 + x*dr;
		    acc(m) = (phi(m+(Nt+1)) + phi(m-(Nt+1)) - 2.0*phi(m))/pow(dr,2.0) + (phi(m+(Nt+1))-phi(m-(Nt+1)))/r/dr - phi(m) + pow(phi(m),3.0);
		    acc(m) *= sigma;
		    linErgField(u-1) +=  4.0*pi*pow(r,2.0)*dr*0.5*pow(phi(m)-phi(m-1),2.0)/pow(dt,2.0);
			linErgField(u) += 4.0*pi*(r*(r+dr)*0.5*pow(phi(m+(Nt+1))-phi(m),2.0)/dr + pow(r,2.0)*0.5*dr*pow(phi(m),2.0));
			nonLinErg(u) += 4.0*pi*pow(r,2.0)*0.25*pow(phi(m),4.0)*dr;
		}
	}
	
	for (unsigned int k=0; k<(Nt+1); k++) {
		erg(k) = linErgField(k) + nonLinErg(k);
	}
	
	for (unsigned int k=0; k<ps_run.N; k++) {
		unsigned int u = Nt + k*(Nt+1);
		double r = r0 + k*dr;
		kineticT += 4.0*pi*dr*pow(r,2.0)*0.5*pow(phi(u)-phi(u-1),2.0)/pow(dt,2.0);
		if (k<N) kineticS += 4.0*pi*r*(r+dr)*0.5*pow(phi(u+(Nt+1))-phi(u),2.0)/dr;
		massTerm += 4.0*pi*dr*pow(r,2.0)*0.5*pow(phi(u-1),2.0);
	}
	
	for (unsigned int k=1; k<ps_run.N; k++) {
		double momtm = k*pi/(r1-r0);
		double freqSqrd = 1.0+pow(momtm,2.0);
		double Asqrd, integral1 = 0.0, integral2 = 0.0;
		for (unsigned int l=0; l<ps_run.N; l++) {
			double r = r0 + l*dr;
			unsigned int m = (Nt-1) + l*(Nt+1);
			integral1 += dr*r*phi(m)*pow(2.0/(r1-r0),0.5)*sin(momtm*r);
			integral2 += dr*r*(phi(m+1)-phi(m))*pow(2.0/(r1-r0),0.5)*sin(momtm*r)/dt;
		}
		Asqrd = pow(integral1,2.0) + pow(integral2,2.0)/freqSqrd;
		linErgContm += 2.0*pi*Asqrd*freqSqrd;
		linNumContm += 2.0*pi*Asqrd*pow(freqSqrd,0.5);
	}
	
	if (j==0) {
		phiC = phi;
	}
	else if (j==1) {
		phiA = phi;
		ergA = erg(Nt-1);
		nonLinErgA = nonLinErg(Nt-1);
		linErgFieldA = linErgField(Nt-1);
		if (testLinear) {
			linearizationA = Eigen::VectorXd::Zero(Nt+1);
			bool nonLin = true;
			for (unsigned int k=0; k<(Nt+1); k++) {
				linearizationA(k) = absDiff(erg(k),linErgField(k));
				if (linearizationA(k)>closenessLin) nonLin = true;
				if (linearizationA(k)<closenessLin && nonLin) {
					Tlin = k*dt;
					nonLin = false;
				}
			}
			for (unsigned int j=0;j<ps_run.N;j++){
				unsigned int n = (ps_run.N-1-j)*(Nt+1)+(unsigned int)(Tlin/dt);
				double rj = r0+(ps_run.N-1-j)*dr;
				if (abs(phi(n))>closenessEdge && abs(Redge)<1.0e-16)
		        	Redge = rj;
				for (unsigned int k=0;k<ps_run.N;k++){
					unsigned int l = j*(Nt+1)+Nt-1;
					unsigned int m = k*(Nt+1)+Nt-1;
					double r = r0 + j*dr;
					double s = r0 + k*dr;
					// the would be imaginary parts of the following expression cancel as omega and Eomega are symmetric
					// so wlog we have dropped them
						linErgA += 0.5*r*s*\
							(omega_2(j,k)*phi(l)*phi(m)+omega_0(j,k)*(phi(l+1)-phi(l))*(phi(m+1)-phi(m))/pow(dt,2.0));
						linNumA += 0.5*r*s*\
							(omega_1(j,k)*phi(l)*phi(m)+omega_m1(j,k)*(phi(l+1)-phi(l))*(phi(m+1)-phi(m))/pow(dt,2.0));
				}
			}
		}
	}
	j++;
		
} // end of while j<2 loop

/* ---------------------------------------------------------------------------------------------
	7. compiling and printing results
---------------------------------------------------------------------------------------------*/

if (testLinear) {
	Filename linFile = (string)(prefix + "linearization" + suffix);
	SaveOptions so_simple;
	so_simple.paramsIn = ps_run; so_simple.paramsOut = ps_run;
	so_simple.vectorType = SaveOptions::simple;
	so_simple.extras = SaveOptions::none;
	so_simple.printMessage = false;

	vec tVec = Eigen::VectorXd::Zero(Nt+1);
	for (unsigned int t=0;t<(Nt+1);t++){
		tVec(t) = t*dt;
	}
	save(linFile,so_simple,tVec);
	so_simple.vectorType = SaveOptions::append;
	save(linFile,so_simple,linearizationA);
	so_simple.vectorType = SaveOptions::simple;
	
	momTest = linErgA*dr/linNumA/pi;
	printf("linearization printed:  %39s\n",linearizationFile.c_str());
	printf("Tlin(%6.4f)          = %6.4f\n",closenessLin,Tlin);
	printf("Redge(%6.4f)         = %6.4f\n",closenessEdge,Redge);
	printf("momTest               = %6.4f\n",momTest);
	if (changeParams) {
		//string inputsOut = "data/" + timenumber + "inputsPi3_" + loopIn;
		string inputsOut = "inputs";
		int Nmom = (int)(linErgA*Redge/pi/linNumA/closenessMom) + 1;
		int Nbmom = (int)(Nmom*4*Tbin/Redge);
		int Nalin = (int)(Nbmom*Tlin/Tbin);
		printf("L changed to :        %6.4g\n",Redge);
		printf("N changed to :        %6i\n",Nmom);
		printf("Na changed to:        %6i\n",Nalin); 
		printf("Nb changed to:        %6i\n\n",Nbmom); 
		changeInputs("data/temp1","LoR",numberToString<double>(Redge/10.0),inputsF);
		changeInputs("data/temp2","N",numberToString<int>(Nmom),"data/temp1");
		changeInputs("data/temp1","Na",numberToString<int>(Nalin),"data/temp2");
		changeInputs("data/temp2","Nb",numberToString<int>(Nbmom),"data/temp1");
		copyFile("data/temp2",inputsOut);
	}
}
else if (!testTunnel) {	
	vec tVec((Nain+Nbin+Ncin)*Nin), rVec((Nain+Nbin+Ncin)*Nin);
	for (unsigned int t=0;t<(Nain+Nbin+Ncin);t++) {
		for (unsigned int r=0; r<Nin; r++) {
			unsigned int j= t + r*(Nain+Nbin+Ncin);
			tVec(j) = t*dtin;
			rVec(j) = r0 + r*drin;
		}
	}
	
	// constructing input to main
	vec mainIn(ps_in.NT*ps_in.N);
	vec phiAOut, phiCOut;
	phiAOut = interpolate(phiA,ps_run.Na+1,ps_run.N,ps_in.Na+1,ps_in.N);
	phiCOut = interpolate(phiC,ps_run.Nc+1,ps_run.N,ps_in.Nc+1,ps_in.N);
	for (unsigned int j=0;j<ps_in.NT;j++) {
		for (unsigned int k=0; k<ps_in.N; k++) {
			unsigned int l = j+k*ps_in.NT, m;
			double r = r0 + k*drin;
			if (j<ps_in.Na) {
				m = (ps_in.Na-j)+k*(ps_in.Na+1);
				mainIn(l) = phiAOut(m)*r;
			}
			else if (j<(ps_in.Na+ps_in.Nb)) {
				m = (j-ps_in.Na)+k*ps_in.Nb;
				mainIn(l) = phiBC(m);
			}
			else {
		        m = j - ps_in.Na - ps_in.Nb + 1 + k*(ps_in.Nc+1);
				mainIn(l) = phiCOut(m)*r;
			}
		}
	}
	string mainInFile = "data/" + timenumber + "tpip_" + loopIn + ".dat";
	printThreeVectors(mainInFile,tVec,rVec,mainIn);
	gp(mainInFile,"pi3.gp");

	printf("%8s%8s%8s%8s%8s%8s\n","N","Na","Nb","Nc","L","Tb");
	printf("%8i%8i%8i%8i%8g%8g\n",ps_in.N,ps_in.Na,ps_in.Nb,ps_in.Nc,ps_in.L,ps_in.Tb);
	printf("\n");
	
	printf("Input:                  %39s\n",filename.c_str());
	printf("tpip printed:                %39s pics/pi3.png\n",mainInFile.c_str());
}
else {
	printf("Input:              %39s\n",filename.c_str());
	printf("Tinf              = %8.4f\n",Tinf);
	}
				printf("kineticT          = %8.4f\n",kineticT);
				printf("kineticS          = %8.4f\n",kineticS);
				printf("massTerm          = %8.4f\n",massTerm);
				printf("erg(0)            = %8.4f\n",ergA);
				printf("nonLinErgA(0)     = %8.4f\n",nonLinErgA);
				printf("linErgFieldA(0)   = %8.4f\n",linErgFieldA);
				printf("linErgContmA(0)   = %8.4f\n",linErgContm);
if (testLinear) printf("linErgA           = %8.4f\n",linErgA);
				printf("linNumContmA(0)   = %8.4f\n",linNumContm);
if (testLinear) printf("linNumA           = %8.4f\n\n",linNumA);

double finalTest = linErgFieldA;
if (testTunnel) {
	if ( !isfinite(finalTest) ) return 0;
	else return 1;
}
else return 0;
}

/*---------------------------------------------------------------------------------------------
	piEvolve
		- something to compile results to get periodic instanton for Pot=3
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
#include "omega.h"
#include "parameters.h"
#include "print.h"
#include "print3.h"


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
bool testTunnel = false, testLinear = false, changeParams = false;
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
string timenumberIn = "000";
string loopIn = "0";
if (argc == 2) timenumberIn = argv[1];
else if (argc % 2 && argc>1) {
	for (unsigned int j=0; j<(unsigned int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("tn")==0) timenumberIn = argv[2*j+2];
		else if (id.compare("loop")==0 || id.compare("l")==0) loopIn = argv[2*j+2];
		else if (id.compare("test")==0||id.compare("tunnel")==0||id.compare("tt")==0) testTunnel=(bool)stringToNumber<uint>(argv[2*j+2]);
		else if ((id.substr(0,3)).compare("lin")==0) testLinear = (bool)stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("closeness")==0 || id.compare("close")==0) closenessLin = stringToNumber<double>(argv[2*j+2]);
		else if (id.compare("changeParams")==0) changeParams = (bool)stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("N")==0) ps_run.changeParameters("N",stringToNumber<uint>(argv[2*j+2]));
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
string prefix = "data/"+timenumberIn;
string suffix = "_loop_"+loopIn+".dat";
Filename inputsFile = (string)(prefix+"inputsP"+suffix);
inputsFile.Suffix = "";
ps_in.load(inputsFile);

// chosing p_run
ps_run = ps_in;
/*ps_run.changeParameters("N",(uint)300);
	{
	double dt = ps_run.a/5.0;
	ps_run.changeParameters("Na",(uint)(ps_in.Ta/dt));
	ps_run.changeParameters("Nb",(uint)(ps_in.Tb/dt));
	ps_run.changeParameters("Nc",(uint)(ps_in.Tc/dt));
	}*/

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
	cerr << "ps_run.b too large. ps_run.b = " << ps_run.b << ", ps_run.a = " << ps_run.a << endl;
	return 1;
}
if (abs(ps_in.Ta)<2.5 && !testTunnel && !testLinear) {
	cerr << "Ta too small. Ta = " << ps_in.Ta << endl;
	return 1;
}

/* ---------------------------------------------------------------------------------------------
	4. loading phi on BC
		- load options
		- phiBC
		- printing parameters in
---------------------------------------------------------------------------------------------*/

// load options
Filename pFile = (string)(prefix+"p"+suffix);
SaveOptions so_in;
so_in.paramsIn = ps_in;
so_in.paramsOut = ps_in;
so_in.vectorType = SaveOptions::realB; // input is actually complex but imaginary parts are zero, don't need to load them
so_in.extras = SaveOptions::coords;
so_in.zeroModes = 1;
so_in.printMessage = true;
so_in.printType = SaveOptions::ascii;
so_in.column = 0;

// phiBC
vec phiBC;
load(pFile,so_in,phiBC);

// printing parameters in
cout << "parameters in:" << endl;
ps_in.print();

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
	so_simple.printType = SaveOptions::ascii;
	{
		Filename omegaM1F, omega0F, omega1F, omega2F, modesF, freqsF, freqsExpF; // Filename works as FilenameAttributes
		omegaM1F = (string)("data/stable/omegaM1_Pot_"+numberToString<uint>(ps_run.Pot)+"_N_"+numberToString<uint>(ps_run.N)\
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
		- declaring vectors and doubles holding results
		- setting direction and Nt for run
		- setting vector sizes and initializing results to zero
		- initializing phi, linErg etc from phiBC
		- applying initial conditions
		- stepping in time
		- evaluating energy expressions, including continuum approx to linErg
		- compiling outputs from run
		
---------------------------------------------------------------------------------------------*/

// declaring vectors and doubles holding results
vec phiA, phiC, linearizationA;
double linErgContm, linNumContm, nonLinErgA, linErgFieldA, ergA, linErgA, linNumA;
double kineticT, kineticS, massTerm;
double Tinf = 0.0;

uint run=0, Nt;
while(run<2) {
	// setting direction and Nt for run
	if (testLinear) {
		Nt = 10*ps_run.N;
		run=1;
		direction = -1;
	}
	else if (testTunnel) {
		Nt = 10*ps_run.N;
		run=1;
		direction = 1;
	}
	else if (run==0) {
		Nt = ps_run.Nc;
		direction = 1;
	}
	else if (run==1) {
		Nt = ps_run.Na;
		direction = -1;
	}
	
	if (Nt==0) {
		cerr << "Nt==0 on run " << run << endl;
		run++;
		continue;
	}
	
	// setting vector sizes and initializing results to zero
	vec phi((Nt+1)*ps_run.N), vel((Nt+1)*ps_run.N), acc((Nt+1)*ps_run.N);
	vec nonLinErg(Nt+1), linErgField(Nt+1), erg(Nt+1);
	linErgField = Eigen::VectorXd::Zero(Nt+1);
	nonLinErg = Eigen::VectorXd::Zero(Nt+1);
	erg = Eigen::VectorXd::Zero(Nt+1);
	linErgContm = 0.0, linNumContm = 0.0;
	linErgA = 0.0, linNumA = 0.0;
	kineticT = 0.0, kineticS = 0.0, massTerm = 0.0;
	
	// initializing phi from phiBC
	vec initial;
	// phiBC1d
	vec phiBC1d(ps_in.N);
	for (uint j=0; j<ps_in.N; j++) {
		lint l = j*ps_in.Nb;
		lint m = (direction == 1 ? ps_in.Nb-1 + l : l);
		phiBC1d(j) = phiBC(m);
	}
	initial = interpolate1d(phiBC1d,ps_in.N,ps_run.N);

	//initial condition 1)
	//defining phi on initial time slice
	for (unsigned int j=0;j<ps_run.N;j++) {
		lint l = j*(Nt+1);
		double r = ps_run.r0 + j*ps_run.a;
		phi(l) = initial(j)/r;
	}
	
	//initialising linErg and linNum	
	for (unsigned int x=0;x<ps_run.N;x++) {
		lint j = x*(Nt+1);
		double r = ps_run.r0 + x*ps_run.a, eta;
		eta = ((x==0 || x==(ps_run.N-1)) ? 0.5 : 1.0);
		nonLinErg(0) += 4.0*PI*pow(r,2.0)*eta*0.25*pow(phi(j),4.0)*ps_run.a;
		linErgField(0) += 4.0*PI*pow(r,2.0)*eta*0.5*pow(phi(j),2.0)*ps_run.a;
		if (x<(ps_run.N-1)) linErgField(0) += 4.0*PI*r*(r+ps_run.a)*0.5*pow(phi(j+(Nt+1))-phi(j),2.0)/ps_run.a;
	}

	//initial condition 2)
	//intitialize velocity, dphi/dt(x,0) = 0;
	for (unsigned int j=0; j<ps_run.N; j++) {
		lint l = j*(Nt+1);
		vel(l) = 0.0;
	}

	//boundary condition 1)
	//phi=0.0 at r=L
	for (unsigned int j=0;j<(Nt+1);j++) {
		lint l = (ps_run.N-1)*(Nt+1) + j;
		phi(l) = 0.0;
	}

	//boundary condition 2)
	//initialize acc using phi and expression from equation of motion
	/*the unusual d^2phi/dr^2 term and the absence of the first derivative term
	are due to the boundary condition 2) which, to second order, is phi(t,-ps_run.a) = phi(t,ps_run.a)*/
	acc(0) = 2.0*(phi(Nt+1) - phi(0))/pow(ps_run.a,2.0) - phi(0) + pow(phi(0),3.0);
	acc(0) *= sigma*0.5; //as initial time slice, generated from taylor expansion and equation of motion
	for (unsigned int j=1; j<(ps_run.N-1); j++) {
		lint l = j*(Nt+1);
		double r = ps_run.r0 + ps_run.a*j;
		acc(l) = (phi(l+(Nt+1)) + phi(l-(Nt+1)) - 2.0*phi(l))/pow(ps_run.a,2.0) + (phi(l+(Nt+1))-phi(l-(Nt+1)))/r/ps_run.a - phi(l) + pow(phi(l),3.0);
		acc(l) *= sigma*0.5;
	}

	//run loop
	for (unsigned int u=1; u<(Nt+1); u++) {
		for (unsigned int x=0; x<(ps_run.N-1); x++) { //don't loop over last x position as fixed by boundary condition 1)
		    lint m = u+x*(Nt+1);
		    vel(m) = vel(m-1) + ps_run.b*acc(m-1);
		    phi(m) = phi(m-1) + ps_run.b*vel(m);
		    if (testTunnel) {
		    	double testInf = phi(m);
		    	if (!isfinite(testInf) && abs(Tinf)<MIN_NUMBER) Tinf = u*ps_run.b;
		    }
		}
		acc(u) = 2.0*(phi(u+Nt+1) - phi(u))/pow(ps_run.a,2.0) - phi(u) + pow(phi(u),3.0);
		acc(u) *= sigma;
		linErgField(u-1) += 4.0*PI*ps_run.a*pow(ps_run.r0,2.0)*0.5*pow(phi(u)-phi(u-1),2.0)/pow(ps_run.b,2.0);
		linErgField(u-1) += 4.0*PI*ps_run.a*pow(ps_run.L,2.0)*0.5*pow(phi(u+(ps_run.N-1)*(Nt+1))-phi(u+(ps_run.N-1)*(Nt+1)-1),2.0)/pow(ps_run.b,2.0);
		linErgField(u) += 4.0*PI*( ps_run.r0*(ps_run.r0+ps_run.a)*0.5*pow(phi(u+(Nt+1))-phi(u),2.0)/ps_run.a + ps_run.a*pow(ps_run.r0,2.0)*0.5*pow(phi(u),2.0) );
		linErgField(u) += 4.0*PI*ps_run.a*pow(ps_run.L,2.0)*0.5*pow(phi(u+(ps_run.N-1)*(Nt+1)),2.0);
		nonLinErg(u) += 4.0*PI*ps_run.a*pow(ps_run.r0,2.0)*0.25*pow(phi(u),4.0);
		nonLinErg(u) += 4.0*PI*ps_run.a*pow(ps_run.L,2.0)*0.25*pow(phi(u+(ps_run.N-1)*(Nt+1)),4.0);
		for (unsigned int x=1; x<(ps_run.N-1); x++) {
		    unsigned int m = u+x*(Nt+1);
		    double r = ps_run.r0 + x*ps_run.a;
		    acc(m) = (phi(m+(Nt+1)) + phi(m-(Nt+1)) - 2.0*phi(m))/pow(ps_run.a,2.0) + (phi(m+(Nt+1))-phi(m-(Nt+1)))/r/ps_run.a - phi(m) + pow(phi(m),3.0);
		    acc(m) *= sigma;
		    linErgField(u-1) +=  4.0*PI*pow(r,2.0)*ps_run.a*0.5*pow(phi(m)-phi(m-1),2.0)/pow(ps_run.b,2.0);
			linErgField(u) += 4.0*PI*(r*(r+ps_run.a)*0.5*pow(phi(m+(Nt+1))-phi(m),2.0)/ps_run.a + pow(r,2.0)*0.5*ps_run.a*pow(phi(m),2.0));
			nonLinErg(u) += 4.0*PI*pow(r,2.0)*0.25*pow(phi(m),4.0)*ps_run.a;
		}
	}
	
	// forming erg
	for (unsigned int k=0; k<(Nt+1); k++) {
		erg(k) = linErgField(k) + nonLinErg(k);
	}
	
	// getting kineticT and massTerm
	for (unsigned int k=0; k<ps_run.N; k++) {
		lint u = Nt + k*(Nt+1);
		double r = ps_run.r0 + k*ps_run.a;
		kineticT += 4.0*PI*ps_run.a*pow(r,2.0)*0.5*pow(phi(u)-phi(u-1),2.0)/pow(ps_run.b,2.0);
		if (k<(ps_run.N-1)) kineticS += 4.0*PI*r*(r+ps_run.a)*0.5*pow(phi(u+(Nt+1))-phi(u),2.0)/ps_run.a;
		massTerm += 4.0*PI*ps_run.a*pow(r,2.0)*0.5*pow(phi(u-1),2.0);
	}
	
	// finding contm expression for linErg and linNum
	for (unsigned int k=1; k<ps_run.N; k++) {
		double momtm = k*PI/(ps_run.L-ps_run.r0);
		double freqSqrd = 1.0+pow(momtm,2.0);
		double Asqrd, integral1 = 0.0, integral2 = 0.0;
		for (unsigned int l=0; l<ps_run.N; l++) {
			double r = ps_run.r0 + l*ps_run.a;
			lint m = (Nt-1) + l*(Nt+1);
			integral1 += ps_run.a*r*phi(m)*pow(2.0/(ps_run.L-ps_run.r0),0.5)*sin(momtm*r);
			integral2 += ps_run.a*r*(phi(m+1)-phi(m))*pow(2.0/(ps_run.L-ps_run.r0),0.5)*sin(momtm*r)/ps_run.b;
		}
		Asqrd = pow(integral1,2.0) + pow(integral2,2.0)/freqSqrd;
		linErgContm += 2.0*PI*Asqrd*freqSqrd;
		linNumContm += 2.0*PI*Asqrd*pow(freqSqrd,0.5);
	}
	
	// assigning phi to phiC or phiA
	if (run==0) {
		phiC = phi;
	}
	else if (run==1) {
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
					Tlin = k*ps_run.b;
					nonLin = false;
				}
			}
			for (unsigned int j=0;j<ps_run.N;j++){
				lint n = (ps_run.N-1-j)*(Nt+1)+(unsigned int)(Tlin/ps_run.b);
				double rj = ps_run.r0+(ps_run.N-1-j)*ps_run.a;
				if (abs(phi(n))>closenessEdge && abs(Redge)<1.0e-16)
		        	Redge = rj;
				for (unsigned int k=0;k<ps_run.N;k++){
					unsigned int l = j*(Nt+1)+Nt-1;
					unsigned int m = k*(Nt+1)+Nt-1;
					double r = ps_run.r0 + j*ps_run.a;
					double s = ps_run.r0 + k*ps_run.a;
					// the would be imaginary parts of the following expression cancel as omega and Eomega are symmetric
					// so wlog we have dropped them
						linErgA += 0.5*r*s*\
							(omega_2(j,k)*phi(l)*phi(m)+omega_0(j,k)*(phi(l+1)-phi(l))*(phi(m+1)-phi(m))/pow(ps_run.b,2.0));
						linNumA += 0.5*r*s*\
							(omega_1(j,k)*phi(l)*phi(m)+omega_m1(j,k)*(phi(l+1)-phi(l))*(phi(m+1)-phi(m))/pow(ps_run.b,2.0));
				}
			}
		}
	}
	run++;		
} // end of while run<2 loop

/* ---------------------------------------------------------------------------------------------
	7. compiling and printing results
		- doing one of the following:
			- testLinear
			- printing for main
			- testTunnel
		- printing erg, linErg etc.
---------------------------------------------------------------------------------------------*/

if (testLinear) {
	Filename linFile = (string)(prefix + "linearization" + suffix);
	SaveOptions so_simple;
	so_simple.paramsIn = ps_run; so_simple.paramsOut = ps_run;
	so_simple.vectorType = SaveOptions::simple;
	so_simple.extras = SaveOptions::none;
	so_simple.printMessage = true;

	vec tVec = Eigen::VectorXd::Zero(Nt+1);
	for (unsigned int t=0;t<(Nt+1);t++){
		tVec(t) = t*ps_run.b;
	}
	save(linFile,so_simple,tVec);
	so_simple.printMessage = false;
	so_simple.vectorType = SaveOptions::append;
	save(linFile,so_simple,linearizationA);
	so_simple.vectorType = SaveOptions::simple;
	
	PlotOptions po_simple;
	Filename linPlotFile = linFile;
	linPlotFile.Suffix = ".png";
	po_simple.output = linPlotFile;
	po_simple.printMessage = true;
	po_simple.column = 1;
	po_simple.column2 = 2;
	plot(linFile,po_simple);
	
	
	momTest = linErgA*ps_run.a/linNumA/PI;
	printf("Tlin(%6.4f)      =     %6.4f\n",closenessLin,Tlin);
	printf("Redge(%6.4f)     =     %6.4f\n",closenessEdge,Redge);
	printf("momTest           =   %6.4f\n",momTest);
	if (changeParams) {
		Filename inputsOut = (string)"inputsP";
		uint Nmom = (uint)(linErgA*Redge/PI/linNumA/closenessMom) + 1;
		uint Nbmom = (uint)(Nmom*4*ps_in.Tb/Redge);
		uint Nalin = (uint)(Nbmom*Tlin/ps_in.Tb);
		printf("L changed to :        %6.4g\n",Redge);
		printf("N changed to :        %6i\n",Nmom);
		printf("Na changed to:        %6i\n",Nalin); 
		printf("Nb changed to:        %6i\n\n",Nbmom); 
		ps_in.changeParameters("LoR",Redge/10.0);
		ps_in.changeParameters("N",Nmom);
		ps_in.changeParameters("Na",Nalin);
		ps_in.changeParameters("Nb",Nbmom);
		ps_in.save(inputsOut);
	}
}
else if (!testTunnel) {	
	
	cout << "parameters run:" << endl;
	ps_run.print();
	
	// constructing input to main
	vec mainIn(2*ps_in.NT*ps_in.N+2);
	vec phiAOut, phiCOut;
	Parameters ps_run_a = ps_run, ps_run_c = ps_run, ps_in_a = ps_in, ps_in_c = ps_in;
	ps_run_a.NT = ps_run.Na+1;
	ps_run_c.NT = ps_run.Nc+1;
	ps_in_a.NT = ps_in.Na+1;
	ps_in_c.NT = ps_in.Nc+1;
	if (ps_in.Na>0)
		phiAOut = interpolateReal(phiA,ps_run_a,ps_in_a);
	if (ps_in.Nc>0)
		phiCOut = interpolateReal(phiC,ps_run_c,ps_in_c);
	for (unsigned int j=0;j<ps_in.NT;j++) {
		for (unsigned int k=0; k<ps_in.N; k++) {
			unsigned int l = j+k*ps_in.NT, m;
			double r = ps_in.r0 + k*ps_in.a;
			mainIn(2*l+1) = 0.0;
			if (j<ps_in.Na) {
				m = (ps_in.Na-j)+k*(ps_in.Na+1);
				mainIn(2*l) = phiAOut(m)*r;
			}
			else if (j<(ps_in.Na+ps_in.Nb)) {
				m = (j-ps_in.Na)+k*ps_in.Nb;
				mainIn(2*l) = phiBC(m);
			}
			else {
		        m = j - ps_in.Na - ps_in.Nb + 1 + k*(ps_in.Nc+1);
				mainIn(2*l) = phiCOut(m)*r;
			}
		}
	}
	mainIn(2*ps_in.NT*ps_in.N) = 0.5;
	mainIn(2*ps_in.NT*ps_in.N+1) = 0.5;
	Filename mainInFile = (string)(prefix + "tp" + suffix);
	SaveOptions so_tp;
	so_tp.paramsIn = ps_in;
	so_tp.paramsOut = ps_in;
	so_tp.vectorType = SaveOptions::complex;
	so_tp.extras = SaveOptions::coords;
	so_tp.zeroModes = 0;
	so_tp.printMessage = true;
	so_tp.printType = SaveOptions::ascii;
	save(mainInFile,so_tp,mainIn);
	mainInFile.ID = "tpnew";
	saveVectorAscii(mainInFile,mainIn);
	mainInFile.ID = "tp";
	
	PlotOptions po_tp;
	po_tp.gp = "gp/repi.gp";
	po_tp.style = "points";
	Filename plotFile = mainInFile;
	plotFile.Suffix = ".png";
	po_tp.output = plotFile;
	po_tp.printMessage = true;
	
	plot(mainInFile,po_tp);

}
else {
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

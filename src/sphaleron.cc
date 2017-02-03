/* -------------------------------------------------------------------------------------------------------------------------
	sphaleron
		- program to find sphaleron solution by solving boundary value ode by shooting method
------------------------------------------------------------------------------------------------------------------------- */
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
#include <ctype.h>
#include <cstring> //for memcpy
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "simple.h"
#include "folder.h"
#include "omega.h"
#include "parameters.h"
#include "print.h"
#include "sphaleron_fns.h"

using namespace std;

/*----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------
	CONTENTS
		1 - defining key quantities
		2 - getting argv inputs
		3 - shooting method
		4 - calculating energy
		5 - matrix of linear fluctuations
		6 - omega
		7 - propagating along minkowskian time
		8 - printing results
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/

int main(int argc, char** argv) {

/*----------------------------------------------------------------------------------------------------------------------------
	1. defining key quantities
----------------------------------------------------------------------------------------------------------------------------*/

void_params paramsVoid;
gsl_odeiv2_system sys = {func, jac, 4, &paramsVoid};

double F = 1.0, dF;
double aim = 0.0;
double closeness = 1.0e-8;
double r0 = 1.0e-15, r1 = 5.0;
unsigned int N = 1e3;
unsigned int runsCount = 0;

/* -------------------------------------------------------------------------------------------------------------------------
	2. getting argv inputs
-------------------------------------------------------------------------------------------------------------------------*/
bool allIWantIsD1AndD2 = true;

if (argc>2 && argc%2) {
	for (int j=0; j<(int)(argc/2); j++) {
		string temp1 = argv[2*j+1];
		string temp2 = argv[2*j+2];
		if (temp1[0]=='-') temp1 = temp1.substr(1);
		if (temp1.compare("r1")==0) r1 = stringToNumber<double>(temp2);
		else if (temp1.compare("r0")==0) r0 = stringToNumber<double>(temp2);
		else if (temp1.compare("N")==0) N = stringToNumber<unsigned int>(temp2);
		else if (temp1.compare("allIWantIsD1AndD2")==0) allIWantIsD1AndD2 = (stn<uint>(temp2)!=0);
		else {
			cerr << "input " << temp1 << " not understood" << endl;
			return 1;
		}
	}
}
else if (!argc%2)
	cerr << "inputs not understood" << endl;

double dr = r1-r0;
dr /= (double)N;

Parameters ps;
ps.Pot = 3;
ps.r0 = r0;
ps.L = r1;
ps.N = N+1;
ps.a = dr;

/* -------------------------------------------------------------------------------------------------------------------------
	3. shooting method
		- initializing solution vectors
		- initial guess
		- shooting loop
		- printing solution
-------------------------------------------------------------------------------------------------------------------------*/

// initializing solution vectors
vec y0Vec(N+1), y2Vec(N+1);

// initial guess
double Y0 = 4.337;
//cout << "initial y0: ";
//cin >> Y0;
cout << "r1 = " << r1 << endl << endl;
printf("%16s%16s%16s%16s%16s%16s%16s%16s\n","run","N","y(r1)","yMin","F-aim","Y0Old","Y0New","-F/dF");

// shooting loop
while (abs(F-aim)>closeness) {
	runsCount++;
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_yp_new (&sys,  gsl_odeiv2_step_rk8pd, 1.0e-6, 1.0e-6, 0.0);

	double r = r0, ri;
	double y0[4] = { Y0, 0.0, 1.0, 0.0};
	double y[4];
	memcpy(y, y0, sizeof(y0));
	double yMin = abs(y[0]-aim);
	y0Vec[0] = y[0];
	y2Vec[0] = y[2];
	int status;
	unsigned int i, iMin = 0;
	
	for (i = 1; i <= N; i++) {
		ri =(double)i;
		ri *= dr;
		ri += r0;
		status = gsl_odeiv2_driver_apply (d, &r, ri, y);
		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			printf ("i = %3i, r = %3g\n",i,r);
			break;
		}
		if ((y[0]-aim)<yMin && (y[0]-aim)>0.0) {
			yMin = y[0]-aim;
			iMin = i;
		}
		y0Vec[i] = y[0];
		y2Vec[i] = y[2];
		//printf ("%.5e %.5e %.5e\n", r, y[0], y[1]);
		if ((y[0]-aim)<0.0) {
			iMin = i;
			if ((y[0]-aim)<-0.2) break;
		}
	}
	if (status != GSL_SUCCESS) break;
		
	F = y0Vec[iMin]-aim; //as final boundary condition is y(r1)=0.0;
	dF = y2Vec[iMin];
	printf("%16i%16i%16g%16g%16g%16.12g",runsCount,i,y[0],yMin,F,Y0);
	if (abs(dF)>2.0e-16) {
		Y0 += -F/dF;
		printf("%16.12g%16g\n",Y0,-F/dF);
	}
	gsl_odeiv2_driver_free (d);
	if (i==(N+1)) F = y[0]-aim;
}

// printing solution
string filename = "data/stable/sphaleron_L_"+numberToString<double>(r1)+".dat";
string picname = "pics/sphaleron_L_"+numberToString<double>(r1)+".png";
SaveOptions so_simple;
so_simple.printType = SaveOptions::ascii;
so_simple.printMessage = true;
so_simple.vectorType = SaveOptions::simple;
so_simple.extras = SaveOptions::none;
so_simple.column = 0;
so_simple.zeroModes = 0;
PlotOptions po_simple;
po_simple.column = 1;
po_simple.style = "linespoints";
po_simple.output = picname;
po_simple.printMessage = true;
save(filename,so_simple,y0Vec);
plot(filename,po_simple);

/*----------------------------------------------------------------------------------------------------------------------------
	4. calculating energy
----------------------------------------------------------------------------------------------------------------------------*/

//finding E
double E, Eerror;
E_params paramsE;
paramsE.Y_0 = Y0;
gsl_function E_integrand_gsl;
E_integrand_gsl.function = &E_integrand;
E_integrand_gsl.params = &paramsE;
gsl_integration_workspace *w = gsl_integration_workspace_alloc(1e4);
gsl_integration_qag(&E_integrand_gsl, r0, r1, 1.0e-10, 1.0e-9, 1e4, 4, w, &E, &Eerror);
gsl_integration_workspace_free(w);
if (Eerror>1.0e-8) 	cout << "E error = " << Eerror << endl;
else 				cout << "E = " << E << endl << endl;

/*----------------------------------------------------------------------------------------------------------------------------
	5. matrix of linear fluctuations
----------------------------------------------------------------------------------------------------------------------------*/

// matrix of linear fluctuations
bool printD1AndD2 = true;
if (printD1AndD2) {
	spMat D1(N+1,N+1), D2(N+1,N+1);
	D1.setZero();
	D2.setZero();
	double r = r0;
	for (unsigned int j=0; j<(N+1); j++) {
		if (j==0) {
			D1.insert(j,j) = 1.0/pow(dr,2.0) + (1.0 - 3.0*pow(y0Vec[j],2.0))/2.0;
			D1.insert(j,j+1) = -1.0/pow(dr,2.0);
			D2.insert(j,j) = 1.0/pow(dr,2.0) + 1.0/dr/r + 1.0 - 3.0*pow(y0Vec[j],2.0);
			D2.insert(j,j+1) = -1.0/pow(dr,2.0) - 1.0/dr/r;
			/*D1.insert(j,j) = -1.0; //derivative of phi is zero at r=0
			D1.insert(j,j+1) = 1.0;
			D2.insert(j,j) = -1.0;
			D2.insert(j,j+1) = 1.0;*/
		}
		else if (j==N) {
			D1.insert(j,j) = 1.0/pow(dr,2.0) - 2.0*(1.0-dr/r)/r/dr + (1.0 - 3.0*pow(y0Vec[j],2.0))/2.0;
			D1.insert(j,j-1) = -1.0/pow(dr,2.0) + 2.0*(1.0-dr/r)/r/dr;
			D2.insert(j,j) = 1.0/pow(dr,2.0) - 1.0/dr/r + 1.0 - 3.0*pow(y0Vec[j],2.0);
			D2.insert(j,j-1) = -1.0/pow(dr,2.0) + 1.0/dr/r;
			/*D1.insert(j,j) = 1.0; //phi goes to zero as r->infty
			D2.insert(j,j) = 1.0;*/
		}
		else {
			D1.insert(j,j) = 2.0/pow(dr,2.0) - 2.0*(1.0-dr/r)/r/dr + 1.0 - 3.0*pow(y0Vec[j],2.0);
			D1.insert(j,j+1) = -1.0/pow(dr,2.0);
			D1.insert(j,j-1) = -1.0/pow(dr,2.0) + 2.0*(1.0-dr/r)/r/dr;
			D2.insert(j,j) = 2.0/pow(dr,2.0) + 1.0 - 3.0*pow(y0Vec[j],2.0);
			D2.insert(j,j+1) = -1.0/pow(dr,2.0) - 1.0/dr/r;
			D2.insert(j,j-1) = -1.0/pow(dr,2.0) + 1.0/dr/r;
		}
		r += dr;
	}
	D1.makeCompressed();
	D2.makeCompressed();
	so_simple.extras = SaveOptions::loc;
	save("data/stable/D1_L_"+numberToString<double>(r1)+".dat",so_simple,D1);
	save("data/stable/D2_L_"+numberToString<double>(r1)+".dat",so_simple,D2);
	so_simple.extras = SaveOptions::none;
}
printf("From Matlab: D1 gives omega^2_- = -15.31,\n");
printf("             D2 gives omega^2_- = -15.34\n\n");

if (allIWantIsD1AndD2 && printD1AndD2) return 0;

/*----------------------------------------------------------------------------------------------------------------------------
	6. omega
----------------------------------------------------------------------------------------------------------------------------*/

//deterimining omega matrices for fourier transforms in spatial direction
mat omega(N+1,N+1);
mat Eomega(N+1,N+1);
{
	vec freqs(ps.N), freqs_exp(ps.N);
	mat modes(ps.N,ps.N);
	mat omega_m1(ps.N,ps.N), omega_0(ps.N,ps.N), omega_1(ps.N,ps.N), omega_2(ps.N,ps.N);
	Filename omegaM1F, omega0F, omega1F, omega2F, modesF, freqsF, freqsExpF; // Filename works as FilenameAttributes
	omegaM1F = (string)"data/stable/omegaM1_Pot_3_N_"+numberToString<uint>(ps.N)\
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
	omega = omega_1;
	Eomega = omega_2;
}
	
/*----------------------------------------------------------------------------------------------------------------------------
	7. propagating along minkowskian time
		- defining time related parameters
		- defining solution vectors
		- getting eigVec from file
	
n.b. in order to calculate linErg and linNum
----------------------------------------------------------------------------------------------------------------------------*/

// defining time related parameters
unsigned int Nt = N*3;
double T = 0.8*(r1-r0), amp = -1.0e-2, sigma = 1.0; //set sigma=-1 for euclidean evolution
if ((T-4.0)>1.1*(r1-r0))
		cout << "R is too small compared to T. R = " << r1-r0 << ", T = " << T << endl;
//cout << "amp of negative mode: ";
//cin >> amp;
double dt = T/Nt; //equals Dt
if (dt>0.5*dr) {
	cout << "dt too large. dt = "<< dt << ", dr = " << dr << endl;
	return 0;
}
ps.Tb = T;
ps.Nb = Nt+1;
ps.b = dt;
	
// defining solution vectors
vec phi((Nt+1)*(N+1)), vel((Nt+1)*(N+1)), acc((Nt+1)*(N+1)), eigVec, linNum(Nt+1), linErg(Nt+1), linErgField(Nt+1);
vec nonLinErg(Nt+1), erg(Nt+1);
linErg = Eigen::VectorXd::Zero(Nt+1);
linErgField = Eigen::VectorXd::Zero(Nt+1);
linNum = Eigen::VectorXd::Zero(Nt+1);
nonLinErg = Eigen::VectorXd::Zero(Nt+1);
erg = Eigen::VectorXd::Zero(Nt+1);

//getting eigVec from file
so_simple.paramsOut = ps;
string eigVecFile = "data/stable/eigVec_Pot_3_L_"+numberToString<double>(r1)+".dat";
load(eigVecFile,so_simple,eigVec);
if (eigVec[0]<0) eigVec *= -1.0;

//initial condition 1)
//defining phi on initial time slice
for (unsigned int j=0;j<(N+1);j++) {
	unsigned int l = j*(Nt+1);
	phi(l) = y0Vec[j] + amp*eigVec(j);
}
	
//initialising linErg and linNum	
for (unsigned int x=0;x<(N+1);x++) {
	unsigned int j = x*(Nt+1);
	double r = r0 + x*dr, sigma = 1.0;
	if (x==0 || x==N) sigma = 0.5;
	nonLinErg(0) += 4.0*PI*pow(r,2.0)*sigma*0.25*pow(phi(j),4.0)*dr;
	linErgField(0) += 4.0*PI*pow(r,2.0)*sigma*0.5*pow(phi(j),2.0)*dr;
	if (x<N) {
		linErgField(0) += 4.0*PI*r*(r+dr)*0.5*pow(phi(j+(Nt+1))-phi(j),2.0)/dr;
	}
	for (unsigned int k=0;k<(N+1);k++) {
		unsigned int l = k*(Nt+1);
		linErg(0) += Eomega(x,k)*phi(l)*phi(j);
		linNum(0) += omega(x,k)*phi(l)*phi(j);	
	}
}

//initial condition 2)
//intitialize velocity, dphi/dt(x,0) = 0;
for (unsigned int j=0; j<(N+1); j++) {
	unsigned int l = j*(Nt+1);
    vel(l) = 0.0;
}

//boundary condition 1)
//phi=0.0 at r=L
for (unsigned int j=0;j<(Nt+1);j++) {
	unsigned int l = N*(Nt+1) + j;
	phi(l) = 0.0;
}

//boundary condition 2)
//initialize acc using phi and expression from equation of motion
/*the unusual d^2phi/dr^2 term and the absence of the first derivative term
are due to the boundary condition 2) which, to second order, is phi(t,-dr) = phi(t,dr)*/
acc(0) = 2.0*(phi(Nt+1) - phi(0))/pow(dr,2.0) - phi(0) + pow(phi(0),3.0);
acc(0) *= sigma*0.5; //as initial time slice, generated from taylor expansion and equation of motion
for (unsigned int j=1; j<N; j++) {
	unsigned int l = j*(Nt+1);
	double r = r0 + dr*j;
    acc(l) = (phi(l+(Nt+1)) + phi(l-(Nt+1)) - 2.0*phi(l))/pow(dr,2.0) + (phi(l+(Nt+1))-phi(l-(Nt+1)))/r/dr - phi(l) + pow(phi(l),3.0);
    acc(l) *= sigma*0.5;
}

//A7. run loop
for (unsigned int u=1; u<(Nt+1); u++) {
    for (unsigned int x=0; x<N; x++) { //don't loop over last x position as fixed by boundary condition 1)
        unsigned int m = u+x*(Nt+1);
        vel(m) = vel(m-1) + dt*acc(m-1);
        phi(m) = phi(m-1) + dt*vel(m);
    }
    acc(u) = 2.0*(phi(u+Nt+1) - phi(u))/pow(dr,2.0) - phi(u) + pow(phi(u),3.0);
    acc(u) *= sigma;
    linErgField(u-1) += 4.0*PI*dr*pow(r0,2.0)*0.5*pow(phi(u)-phi(u-1),2.0)/pow(dt,2.0);
    linErgField(u-1) += 4.0*PI*dr*pow(r1,2.0)*0.5*pow(phi(u+N*(Nt+1))-phi(u+N*(Nt+1)-1),2.0)/pow(dt,2.0);
    linErgField(u) += 4.0*PI*( r0*(r0+dr)*0.5*pow(phi(u+(Nt+1))-phi(u),2.0)/dr + dr*pow(r0,2.0)*0.5*pow(phi(u),2.0) );
    linErgField(u) += 4.0*PI*dr*pow(r1,2.0)*0.5*pow(phi(u+N*(Nt+1)),2.0);
    nonLinErg(u) += 4.0*PI*dr*pow(r0,2.0)*0.25*pow(phi(u),4.0);
    nonLinErg(u) += 4.0*PI*dr*pow(r1,2.0)*0.25*pow(phi(u+N*(Nt+1)),4.0);
    for (unsigned int x=1; x<N; x++) {
        unsigned int m = u+x*(Nt+1);
        double r = r0 + x*dr;
        acc(m) = (phi(m+(Nt+1)) + phi(m-(Nt+1)) - 2.0*phi(m))/pow(dr,2.0) + (phi(m+(Nt+1))-phi(m-(Nt+1)))/r/dr - phi(m) + pow(phi(m),3.0);
        acc(m) *= sigma;
        linErgField(u-1) +=  4.0*PI*pow(r,2.0)*dr*0.5*pow(phi(m)-phi(m-1),2.0)/pow(dt,2.0);
		linErgField(u) += 4.0*PI*(r*(r+dr)*0.5*pow(phi(m+(Nt+1))-phi(m),2.0)/dr + pow(r,2.0)*0.5*dr*pow(phi(m),2.0));
		nonLinErg(u) += 4.0*PI*pow(r,2.0)*0.25*pow(phi(m),4.0)*dr;
        for (unsigned int k=0;k<(N+1);k++) {
			unsigned int j = u+k*(Nt+1);
			linErg(u) += Eomega(x,k)*phi(m)*phi(j);
			linNum(u) += omega(x,k)*phi(m)*phi(j);
		}
    }
}
	
for (unsigned int k=0; k<(Nt+1); k++) {
	erg(k) = linErgField(k) + nonLinErg(k);
}
	
double linErgContm = 0.0, linNumContm = 0.0;
for (unsigned int k=1; k<(N+1); k++) {
	double momtm = k*PI/(r1-r0);
	double freqSqrd = 1.0+pow(momtm,2.0);
	double Asqrd, integral1 = 0.0, integral2 = 0.0;
	for (unsigned int l=0; l<(N+1); l++) {
		double r = r0 + l*dr;
		unsigned int m = (Nt-1) + l*(Nt+1);
		integral1 += dr*r*phi(m)*pow(2.0/(r1-r0),0.5)*sin(momtm*r);
		integral2 += dr*r*(phi(m+1)-phi(m))*pow(2.0/(r1-r0),0.5)*sin(momtm*r)/dt;
	}
	Asqrd = pow(integral1,2.0) + pow(integral2,2.0)/freqSqrd;
	linErgContm += 2.0*PI*Asqrd*freqSqrd;
	linNumContm += 2.0*PI*Asqrd*pow(freqSqrd,0.5);
}
	
unsigned int N_print = 100, Nt_print = 100;
vec tVec(Nt_print*N_print), rVec(Nt_print*N_print);
double dtPrint = T/(Nt_print-1.0);
double dxPrint = (r1-r0)/(N_print-1.0);
for (unsigned int t=0;t<Nt_print;t++) {
	for (unsigned int r=0; r<N_print; r++) {
		unsigned int j= t + r*Nt_print;
		tVec(j) = t*dtPrint;
		rVec(j) = r0 + r*dxPrint;
	}
}

string datOutStr = (string)("data/sphaleronLinNum_L_"+numberToString<double>(r1)+".dat");
Filename datOut = datOutStr;
Filename pngOut = datOut;
pngOut.Suffix = ".png";
pngOut.Directory = "pics";

// linNum
save(datOut,so_simple,linNum);
po_simple.output = pngOut;
plot(datOut,po_simple);

// linErg
datOut.ID = "sphaleronLinErg";
pngOut.ID = datOut.ID;
save(datOut,so_simple,linErg);
po_simple.output = pngOut;
plot(datOut,po_simple);

// linErgField
datOut.ID = "sphaleronLinErgField";
pngOut.ID = datOut.ID;
save(datOut,so_simple,linErgField);
po_simple.output = pngOut;
plot(datOut,po_simple);

// nonLinErg
datOut.ID = "sphaleronNonLinErg";
pngOut.ID = datOut.ID;
save(datOut,so_simple,nonLinErg);
po_simple.output = pngOut;
plot(datOut,po_simple);

// erg
datOut.ID = "sphaleronErg";
pngOut.ID = datOut.ID;
save(datOut,so_simple,erg);
po_simple.output = pngOut;
plot(datOut,po_simple);

// phi
datOut.ID = "sphaleronEvo";
pngOut.ID = datOut.ID;

Parameters ps_print;
ps.Pot = 3;
ps.N = 300;
ps.Nb = 300;
ps.Tb = T;
ps.L = r1;

SaveOptions so_p;
so_p.printType = SaveOptions::ascii;
so_p.printMessage = true;
so_p.vectorType = SaveOptions::realB;
so_p.extras = SaveOptions::coords;
so_p.column = 0;
so_p.zeroModes = 0;
so_p.paramsIn = ps;
so_p.paramsOut = ps_print;

PlotOptions po_p;
po_p.printMessage = true;
po_p.gp = "gp/repi.gp";
po_p.style = "points";
po_p.output = pngOut;

save(datOut,so_p,phi);
plot(datOut,po_p);

printf("\nerg(Nt-1) = %8.4f\n",linErgField(Nt-1));
printf("linErgField(Nt-1) = %8.4f\n",linErgField(Nt-1));
printf("nonLinErg(Nt-1) = %8.4f\n",nonLinErg(Nt-1));
printf("linNum(Nt-1) = %8.4f\n",linNum(Nt-1));
printf("linErg(Nt-1) = %8.4f\n",linErg(Nt-1));
printf("linNumContm(Nt-1) = %8.4f\n",linNumContm);
printf("linErgContm(Nt-1) = %8.4f\n\n",linErgContm);

return 0;
}

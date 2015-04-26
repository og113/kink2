/*----------------------------------------------------------------------------------------------------------------------------
	static
		program to solve boundary value problem for static soliton
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
		1 - loading options, closenesses, argv inputs
		2 - getting inputs, checks etc.
		3 - assigning potential fucntions
		4 - defining quantites
		5 - calculating thin wall phi
		6 - beginning newton raphson loop
----------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------*/


int main(int argc, char** argv)
{
/*----------------------------------------------------------------------------------------------------------------------------
	1. loading options, argv inputs
		- loading options
		- argv inputs
		- defining timenumber
----------------------------------------------------------------------------------------------------------------------------*/

// loading opts
Options opts;
opts.load("optionsS");
//opts.print();

// defining timenumber
string timenumber = "";//currentDateTime();

// getting argv inputs
if (argc==2) timenumber = argv[1];
else if (argc % 2 && argc>1) {
	for (unsigned int j=0; j<(int)(argc/2); j++) {
		string id = argv[2*j+1];
		if (id[0]=='-') id = id.substr(1);
		if (id.compare("loops")==0) opts.loops = stringToNumber<uint>(argv[2*j+2]);
		else if (id.compare("printChoice")==0) opts.printChoice = argv[2*j+2];
		else {
			cerr << "input " << id << " unrecognized" << endl;
			return 1;
		}
	}
}
else if (argc != 1) {
	cerr << "must provide an even number of inputs in format '-name value':" << endl;
	for (int j=1; j<argc; j++) cerr << argv[j] << " ";
	cerr << endl;
	return 1;
}

/*----------------------------------------------------------------------------------------------------------------------------
	2. getting inputs, checks etc.
		- loading inputs
		- printing timenumber and parameters
		- clock
		- checks
----------------------------------------------------------------------------------------------------------------------------*/

// loading inputs
Parameters ps;
ps.load("inputsS");
if (ps.pot==3) {
	cerr << "no soliton for potential 3, retry with pot=1 or 2" << endl;
	return 1;
}
	
//printing timenumber and parameters
printf("%12s%12s\n","timenumber: ",timenumber.c_str());
ps.print();

//defining a time and starting the clock
clock_t time;
time = clock();

// loading closenesses
Closenesses closenesses;
closenesses.load("closenessesS");

// declaring Checks
Check checkMass("mass",closenesses.Action);
Check checkSoln("solution",closenesses.Soln);
Check checkSolnMax("solution max",closenesses.SolnMax);
Check checkDelta("delta",closenesses.Delta);
Check checkInv("matrix inversion",closenesses.Inv*ps.N*ps.NT);
Check checkProfile("phi input calculation",closenesses.Profile);

/*----------------------------------------------------------------------------------------------------------------------------
	3. assigning potential functions
		- assigning potential functions
		- assigning preliminary parameter structs
----------------------------------------------------------------------------------------------------------------------------*/

// assigning potential functions
Potential<double> V, dV, ddV;
if (ps.pot==1) {
	V((Potential<double>::PotentialType)&V1<comp>,ps);
	dV((Potential<double>::PotentialType)&dV1<comp>,ps);
	ddV((Potential<double>::PotentialType)&ddV1<comp>,ps);
}
else if (ps.pot==2) {
	V((Potential<double>::PotentialType)&V2<comp>,ps);
	dV((Potential<double>::PotentialType)&dV2<comp>,ps);
	ddV((Potential<double>::PotentialType)&ddV2<comp>,ps);
}
else {
	cerr << "pot option not available, pot = " << ps.pot << endl;
	return 1;
}

// assigning preliminary parameter structs
params_for_V paramsV0  = {ps.epsilon0, ps.A};
//params_for_V paramsV  = {ps.epsilon, ps.A};

/*----------------------------------------------------------------------------------------------------------------------------
	4. defining quantities
		- erg, linErg etc
		- p, minusDS, DDS
----------------------------------------------------------------------------------------------------------------------------*/

// Mass, posZero
double Mass;
double posZero;

//defining some quantities used to stop the Newton-Raphson loop when action stops varying
double Mass_last = 2.0/3.0;
uint runs_count = 0;
uint min_runs = 3;

//initializing phi (=p), DDS and minusDS
vec p(ps.N+1);
p = Eigen::VectorXd::Zero(ps.N+1);
spMat DDS(ps.N+1,ps.N+1);
vec minusDS(ps.N+1);

/*----------------------------------------------------------------------------------------------------------------------------
	5. calculating thin wall phi
		- finding phi profile between minima
		- assigning phi
		- fixing boundary conditions
		- printing input phi
		- posZero
----------------------------------------------------------------------------------------------------------------------------*/
	

//finding phi profile between minima
uint profileSize = ps.N; //more than the minimum
vector<double> phiProfile(profileSize);
vector<double> rhoProfile(profileSize);
double alphaL = -opts.alpha*ps.R, alphaR = opts.alpha*ps.R;
if (1.0<opts.alpha) {
	cerr << "R is too small. Not possible to give thinwall input. It should be >> " << opts.alpha*ps.R;
	return 1;
}
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
for (uint j=0; j<ps.N; j++) {
	double x = -ps.L/2.0+ps.a*(double)j;
	if (x>alphaR) {
		p(j) = ps.minima[0];
	}
	else if (x<alphaL) {
		p(j) = ps.minima[1];
	}
	else if (ps.pot==2) {
		vector<double> rhoPos (profileSize,x);
		for (uint k=0; k<profileSize; k++) {
			rhoPos[k] -= rhoProfile[k];
		}
		uint minLoc = smallestLoc(rhoPos);
        p(j) = phiProfile[minLoc];
	}
	else {
		p(j) = (ps.minima[1]+ps.minima[0])/2.0\
					+ (ps.minima[0]-ps.minima[1])*tanh(x/2.0)/2.0;
	}
	if (j>0) {
		if ((p(j)>0 && p(j-1)<0) || (p(j)<0 && p(j-1)>0))
			posZero = x*abs(p(j-1))/abs(p(j)-p(j-1)) + (x-ps.a)*abs(p(j))/abs(p(j)-p(j-1));
	}
}
p(ps.N) = 0.5; // lagrange multiplier for zero mode

// printing input phi
SaveOptions so_simple;
so_simple.paramsIn = ps; so_simple.paramsOut = ps;
so_simple.vectorType = SaveOptions::simple;
so_simple.extras = SaveOptions::none;
so_simple.printMessage = true;
Filename earlyFile = (string)("data/"+timenumber+"staticE.dat");
save(earlyFile,so_simple,p);

// plotting input phi
PlotOptions po_simple;
po_simple.column = 1;
po_simple.style = "linespoints";
po_simple.printMessage = true;
Filename earlyPlotFile = (string)("data/"+timenumber+"staticE.png");
po_simple.output = earlyPlotFile;
plot(earlyFile,po_simple);

// printing posZero
printf("   posZero: %12.4g\n",posZero);

/*----------------------------------------------------------------------------------------------------------------------------
	6. beginning newton-raphson loop
		- chiX, for fixing zero mode
		- reserving memory for DDS
		- initializing erg etc to zero
----------------------------------------------------------------------------------------------------------------------------*/

while(runs_count<min_runs || !checkSoln.good() || !checkSolnMax.good()) {
	runs_count++;

	//defining the zero mode, and finding posZero
	vec chiX(ps.N);
	chiX = Eigen::VectorXd::Zero(ps.N);
	for (uint j=1; j<(ps.N-1); j++){
	
		chiX(j) = p(j+1)-p(j-1); 
		
		if (runs_count>1) {
			if ((p(j)>0 && p(j-1)<0) || (p(j)<0 && p(j-1)>0)) {
				double x = -ps.L/2.0+ps.a*(double)j;
				posZero = x*abs(p(j-1))/abs(p(j)-p(j-1)) + (x-ps.a)*abs(p(j))/abs(p(j)-p(j-1));
			}
		}  
	}

	// allocating memory for DS, DDS
	minusDS = Eigen::VectorXd::Zero(ps.N+1); //initializing to zero
	DDS.setZero(); //just making sure
	Eigen::VectorXi DDS_to_reserve(ps.N+1);//number of non-zero elements per column
	DDS_to_reserve = Eigen::VectorXi::Constant(ps.N+1,4);
	DDS_to_reserve(0) = 2; //these need to be changed when boundary conditions need to be more compicated
	DDS_to_reserve(1) = 2;
	DDS_to_reserve(ps.N-2) = 2;
	DDS_to_reserve(ps.N-1) = 2;
	DDS_to_reserve(ps.N) = ps.N;
	DDS.reserve(DDS_to_reserve);

	//initializing to zero
	Mass = 0.0;

/*----------------------------------------------------------------------------------------------------------------------------
	9. assigning minusDS, DDS etc
		- beginning loop over lattice points
		- fixing zero mode
		- j=0
		- j=N-1
		- bulk
----------------------------------------------------------------------------------------------------------------------------*/

	// beginning loop over lattice points
	for (uint j = 0; j<ps.N; j++) {	
	
		double Dx = ((j==0 || j==(ps.N-1))? ps.a/2.0: ps.a);
		minusDS(ps.N) 		+= -Dx*p(j)*chiX(j);
		minusDS(j) 			+= -Dx*p(ps.N)*chiX(j);
		DDS.insert(j,ps.N)	= Dx*chiX(j);
		DDS.insert(ps.N,j)	= Dx*chiX(j);
		
		if (j==0) {
			double dx = ps.a;
			Mass += pow(p(j+1)-p(j),2.0)/dx;
			DDS.insert(j,j) = 1.0;
		}
		else if (j==(ps.N-1)) {
			DDS.insert(j,j) = 1.0;
		}
		else {
			double dx = ps.a;
			Mass += pow(p(j+1)-p(j),2.0)/dx;
			minusDS(j) 			+= +p(j+1)/dx + p(j-1)/dx - 2.0*p(j)/dx + dV(p(j))*dx;
			DDS.insert(j,j) 	= 2.0/dx - ddV(p(j))*dx;
			DDS.insert(j,j+1) 	= -1.0/dx;
			DDS.insert(j,j-1) 	= -1.0/dx;
		}
    }
    
/*----------------------------------------------------------------------------------------------------------------------------
	7. solving for delta
		- defining delta
		- analyzing pattern
		- factorizing
		- solving
		- check on inversion
		- p' = p + delta
----------------------------------------------------------------------------------------------------------------------------*/

	//solving for delta in DDS*delta=minusDS, where p' = p + delta		
	vec delta(ps.N+1);
	delta = Eigen::VectorXd::Zero(ps.N+1);
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
	vec diff(ps.N+1);
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
	
/*----------------------------------------------------------------------------------------------------------------------------
	8. printing early
		- p
		- minusDS
		- DDS
		- delta
		- chiX
----------------------------------------------------------------------------------------------------------------------------*/
		
	//printing early if desired	
	if ((opts.printChoice).compare("n")!=0) {
		Filename basic = (string)("data/"+timenumber+"staticE_run_"+numberToString<uint>(runs_count)+".dat");
		so_simple.printMessage = false;
		if ((opts.printChoice).compare("p")==0 || (opts.printChoice).compare("e")==0) {
			Filename pEFile = basic;
			pEFile.ID = "staticpE";
			save(pEFile,so_simple,p);
		}
		if ((opts.printChoice).compare("v")==0 || (opts.printChoice).compare("e")==0) {
			Filename vEFile = basic;
			vEFile.ID = "staticminusDSE";
			save(vEFile,so_simple,minusDS);
		}
		if ((opts.printChoice).compare("m")==0 || (opts.printChoice).compare("e")==0) {
			Filename mEFile = basic;
			mEFile.ID = "staticDDSE";
			save(mEFile,so_simple,DDS);
		}
		if ((opts.printChoice).compare("d")==0 || (opts.printChoice).compare("e")==0) {
			Filename dEFile = basic;
			dEFile.ID = "staticdeltaE";
			save(dEFile,so_simple,delta);
		}
		if ((opts.printChoice).compare("z")==0 || (opts.printChoice).compare("e")==0) {
			Filename zEFile = basic;
			zEFile.ID = "staticchiXE";
			save(zEFile,so_simple,chiX);
		}
	}
		
/*----------------------------------------------------------------------------------------------------------------------------
	9. convergence issues
		- checkReg
		- evaluating norms
		- evaluating checks to stop n-r loop
		- printing convergence tests
		- end of n-r loop
		- checkDT
----------------------------------------------------------------------------------------------------------------------------*/
	
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
	checkMass.add(absDiff(Mass,Mass_last));
	Mass_last = Mass;
	checkSoln.add(normDS/normP);
	checkSolnMax.add(maxDS);
	checkDelta.add(normDelta/normP);
		
	// printing tests to see convergence
	if (runs_count==1) {
		printf("%16s%16s%16s%16s%16s%16s\n","runsCount","MassTest","solTest","solMTest","deltaTest","posZero");
	}
	printf("%16i%16g%16g%16g%16g%16g\n",runs_count,checkMass.back(),checkSoln.back()\
								,checkSolnMax.back(),checkDelta.back(),posZero);
}

/*----------------------------------------------------------------------------------------------------------------------------
10. printing output
	- stopping clock
	- printing results to terminal
	- printing results to file
	- printing (and plotting if vectors):
		- p
		- DDS
		- minusDS
	- return
	
----------------------------------------------------------------------------------------------------------------------------*/

//stopping clock
time = clock() - time;
double realtime = time/1000000.0;

//printing results to terminal
printf("\n");
printf("%8s%8s%8s%12s%14s%14s%14s\n","runs","time","N","L","dE","posZero","Mass");
printf("%8i%8.1g%8i%12g%14.4g%14.4g%14.4g\n",runs_count,realtime,ps.N,ps.L,ps.dE,posZero,Mass);
printf("\n");
printf("%60s\n","----------------------------------------------------------------------------------------------------");

//printing results to file
FILE * staticfile;
staticfile = fopen("data/static.dat","a");
fprintf(staticfile,"%16s%8i%12g%12g%14.4g%14.4g%14.4g\n",timenumber.c_str()\
			,ps.N,ps.L,ps.dE,posZero,Mass,checkSoln.back());
fclose(staticfile);

bool printEverything = false;

string prefix = "data/"+timenumber;
string suffix = ".dat";
so_simple.printMessage = true;

//printing output phi on Euclidean time part
Filename pFile = (string)(prefix+"staticp"+suffix);
save(pFile,so_simple,p);
Filename plotFile = pFile;
plotFile.Directory = "pics";
plotFile.Suffix = ".png";
po_simple.output = plotFile;
plot(pFile,po_simple);

//printing output DDS
pFile.ID = "staticDDS";
save(pFile,so_simple,DDS);

if (printEverything) {
	//printing output minusDS
	pFile.ID = "staticminusDS";
	save(pFile,so_simple,minusDS);
}

if (!checkDelta.good()) {
	return 1;
}

return 0;
}

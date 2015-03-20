	/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for functions to save, load and plot
 -------------------------------------------------------------------------------------------------------------------------*/
 
 #include <fstream>
 #include <iostream>
 #include <sstream>
 #include <iomanip>
 #include <cmath>
 #include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "simple.h" // for countLines
#include "parameters.h"
#include "lattice.h"
#include "print.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. SaveOptions
	2. save
	3. load
	4. plot
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	2. save
		- vec
			- static simpleSaveVec
			- static saveVecB
			- static saveVec
		- cVec
			- static savecVecSimple
			- static savecVecB
			- static savecVec
		- mat
		- cMat
		- spMat
-------------------------------------------------------------------------------------------------------------------------*/

// save vec - simpleSaveVec
static void saveVecSimple(const string& f, const SaveOptions& opts, const vec& v) {
	fstream F;
	F.open(f.c_str(), ios::out);
	F.precision(16);
	F << left;
	uint length = v.size();
	if (opts.vectorType==SaveOptions::complex && !length%2) length = (uint)(length/2);
	for (uint j=0; j<length; j++) {
		if (opts.extras==SaveOptions::loc) 			F << setw(25) << j;
		if (opts.vectorType==SaveOptions::complex)	F << setw(25) << v(2*j) << setw(25) << v(2*j+1);
		else										F << setw(25) << v(j);
													F << endl;
	}
	F.close();
}

// save vec - saveVecB
static void saveVecB (const string& f, const SaveOptions& opts, const vec& v) {
	Parameters pin = opts.paramsIn, pout = opts.paramsOut;
	vec vo;
	if ((pin.N!=pout.N && pout.N!=0) || (pin.Nb!=pout.Nb && pout.Nb!=0)) {
		if (opts.vectorType==SaveOptions::realB)
			vo = interpolateReal(v,pin,pout);
		else if (opts.vectorType==SaveOptions::complexB)
			vo = interpolate(v,pin,pout);
		else {
			cerr << "save error: print vectorType option(" << opts.vectorType << ") not possible" << endl;
			return;
		}
	}
	else {
		vo = v;
	}
	fstream F;
	F.open(f.c_str(), ios::out);
	uint x0 = intCoord(0,1,pout.Nb);
	F.precision(16);
	F << left;
	for (lint j=0; j<pout.N*pout.Nb; j++) {
		uint x = intCoord(j,1,pout.Nb);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		switch(opts.extras) {
			case SaveOptions::none:		break;
			case SaveOptions::loc: 		F << setw(25) << j;
										break;
			case SaveOptions::coords:	F << setw(25) << real(coordB(j,0,pout)) << setw(25) << imag(coordB(j,0,pout));
										F << setw(25) << real(coordB(j,1,pout)); 
										break;
			default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
										break;
		}
		switch(opts.vectorType) {
			case SaveOptions::realB:	F << setw(25) << vo(j) << endl;
										break;
			case SaveOptions::complexB:	F << setw(25) << vo(2*j) << setw(25) << vo(2*j+1)  << endl;
										break;
			default:					cerr << "save error: print vectorType option(" << opts.vectorType << ") not possible" << endl;
										return;
										break;
		}
	}
	if (vo.size()>2*pout.N*pout.Nb) {
		F << endl;
		for (unsigned int k=0; k<(vo.size()-2*pout.N*pout.Nb);k++) {
			F << setw(25) << vo(2*pout.N*pout.Nb+k) << endl;
		}
	}
	F.close();
}

// save vec - saveVec
static void saveVec(const string& f, const SaveOptions& opts, const vec& v) {
	Parameters pin = opts.paramsIn, pout = opts.paramsOut;
	vec vo;
	if ((pin.N!=pout.N && pout.N!=0) || (pin.NT!=pout.NT && pout.NT!=0)) {
		if (SaveOptions::real)
			vo = interpolateReal(v,pin,pout);
		else if (SaveOptions::complex)
			vo = interpolate(v,pin,pout);
	}
	else {
		vo = v;
	}
	fstream F;
	F.open(f.c_str(), ios::out);
	unsigned int x0 = intCoord(0,1,pout.NT);
	F.precision(16);
	F << left;
	for (unsigned long int j=0; j<pout.N*pout.NT; j++) {
		unsigned int x = intCoord(j,1,pout.NT);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		switch(opts.extras) {
			case SaveOptions::none:		break;
			case SaveOptions::loc: 		F << setw(25) << j;
										break;
			case SaveOptions::coords:	F << setw(25) << real(coord(j,0,pout)) << setw(25) << imag(coord(j,0,pout));
										F << setw(25) << real(coord(j,1,pout)); 
										break;
			default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
										break;
		}
		switch(opts.vectorType) {
			case SaveOptions::realB:	F << setw(25) << vo(j) << endl;
										break;
			case SaveOptions::complexB:	F << setw(25) << vo(2*j) << setw(25) << vo(2*j+1)  << endl;
										break;
			default:					cerr << "save error: print vectorType option(" << opts.vectorType << ") not possible" << endl;
										break;
		}
	}
	if (vo.size()>2*pout.N*pout.NT) {
		F << endl;
		for (unsigned int j=0; j<(vo.size()-2*pout.N*pout.NT);j++) {
			F << setw(25) << vo(2*pout.N*pout.NT+j) << endl;
		}
	}
	F.close();
}

// save vec
void save(const string& f, const SaveOptions& opts, const vec& v) {
	switch(opts.vectorType) {
		case SaveOptions::simple:	saveVecSimple(f, opts, v);
									break;
		case SaveOptions::real:		saveVec(f,opts,v);
									break;
		case SaveOptions::complex:	saveVec(f,opts,v);
									break;
		case SaveOptions::realB:	saveVecB(f,opts,v);
									break;
		case SaveOptions::complexB:	saveVecB(f,opts,v);
									break;
		default:					cerr << "save error: print vectorType option(" << opts.vectorType << ") not possible" << endl;
									break;
	}
	if (opts.printMessage) {
		printf("%12s%30s\n","saved: ",f.c_str());
	}
}

// save cVec - simpleSaveVec
static void savecVecSimple(const string& f, const SaveOptions& opts, const cVec& v) {
	fstream F;
	F.open((f).c_str(), ios::out);
	F.precision(16);
	F << left;
	uint length = v.size();
	for (uint j=0; j<length; j++) {
		if (opts.extras==SaveOptions::loc) 	F << setw(25) << j;
											F << setw(25) << real(v(j)) << setw(25) << imag(v(j)) << endl;
	}
	F.close();
}

// save cVec - saveVecB
static void savecVecB (const string& f, const SaveOptions& opts, const cVec& v) {
	Parameters pin = opts.paramsIn, pout = opts.paramsOut;
	cVec vo;
	if ((pin.N!=pout.N && pout.N!=0) || (pin.Nb!=pout.Nb && pout.Nb!=0)) {
		vo = interpolate(v,pin,pout);
	}
	else {
		vo = v;
	}
	fstream F;
	F.open(f.c_str(), ios::out);
	uint x0 = intCoord(0,1,pout.Nb);
	F.precision(16);
	F << left;
	for (lint j=0; j<pout.N*pout.Nb; j++) {
		uint x = intCoord(j,1,pout.Nb);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		switch(opts.extras) {
			case SaveOptions::none:		break;
			case SaveOptions::loc: 		F << setw(25) << j;
										break;
			case SaveOptions::coords:	F << setw(25) << real(coordB(j,0,pout)) << setw(25) << imag(coordB(j,0,pout));
										F << setw(25) << real(coordB(j,1,pout)); 
										break;
			default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
										break;
		}
		F << setw(25) << real(vo(j)) << setw(25) << imag(vo(j))  << endl;
	}
	if (vo.size()>pout.N*pout.Nb) {
		F << endl;
		for (unsigned int k=0; k<(vo.size()-pout.N*pout.Nb);k++) {
			F << setw(25) << vo(pout.N*pout.Nb+k) << endl;
		}
	}
	F.close();
}

// save cVec - saveVec
static void savecVec(const string& f, const SaveOptions& opts, const cVec& v) {
	Parameters pin = opts.paramsIn, pout = opts.paramsOut;
	cVec vo;
	if ((pin.N!=pout.N && pout.N!=0) || (pin.NT!=pout.NT && pout.NT!=0)) {
		vo = interpolate(v,pin,pout);
	}
	else {
		vo = v;
	}
	fstream F;
	F.open(f.c_str(), ios::out);
	uint x0 = intCoord(0,1,pout.NT);
	F.precision(16);
	F << left;
	for (lint j=0; j<pout.N*pout.NT; j++) {
		uint x = intCoord(j,1,pout.NT);
		if (x!=x0) { //this is put in for gnuplot
			F << endl;
			x0 = x;
		}
		switch(opts.extras) {
			case SaveOptions::none:		break;
			case SaveOptions::loc: 		F << setw(25) << j;
										break;
			case SaveOptions::coords:	F << setw(25) << real(coord(j,0,pout)) << setw(25) << imag(coord(j,0,pout));
										F << setw(25) << real(coord(j,1,pout)); 
										break;
			default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
										break;
		}
		F << setw(25) << real(vo(j)) << setw(25) << imag(vo(j))  << endl;
	}
	if (vo.size()>pout.N*pout.NT) {
		F << endl;
		for (uint j=0; j<(vo.size()-pout.N*pout.NT);j++) {
			F << setw(25) << vo(pout.N*pout.NT+j) << endl;
		}
	}
	F.close();
}

// save cVec
void save(const string& f, const SaveOptions& opts, const cVec& v) {
	switch(opts.vectorType) {
		case SaveOptions::simple:	savecVecSimple(f, opts, v);
									break;
		case SaveOptions::real:		savecVec(f,opts,v);
									break;
		case SaveOptions::complex:	savecVec(f,opts,v);
									break;
		case SaveOptions::realB:	savecVecB(f,opts,v);
									break;
		case SaveOptions::complexB:	savecVecB(f,opts,v);
									break;
		default:					cerr << "save error: print vectorType option(" << opts.vectorType << ") not possible" << endl;
									break;
	}
	if (opts.printMessage) {
		printf("%12s%30s\n","saved: ",f.c_str());
	}
}

// save mat
void save(const string& f, const SaveOptions& opts, const mat& m) {
	fstream F;
	F.open(f.c_str(), ios::out);
	F << left;
	F.precision(16);
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			F << setw(25) << j << setw(25) << k;
			F << setw(25) << m(j,k) << endl;
		}
	}
	F.close();
	if (opts.printMessage) {
		printf("%12s%30s\n","saved: ",f.c_str());
	}
}

// save cMat
void save(const string& f, const SaveOptions& opts, const cMat& m) {
	fstream F;
	F.open(f.c_str(), ios::out);
	F << left;
	F.precision(16);
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			F << setw(25) << j << setw(25) << k;
			F << setw(25) << real(m(j,k)) << setw(25) << imag(m(j,k)) << endl;
		}
	}
	F.close();
	if (opts.printMessage) {
		printf("%12s%30s\n","saved: ",f.c_str());
	}
}

// save spMat
void save(const string& f, const SaveOptions& opts, const spMat& m) {
	fstream F;
	F.open(f.c_str(), ios::out);
	F << left;
	F.precision(16);
	for (int l=0; l<m.outerSize(); ++l) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(m,l); it; ++it) {
			F << setw(25) << it.row()+1 << setw(25) << it.col()+1 << setw(25) << it.value() << endl;
		}
	}
	F.close();
	if (opts.printMessage) {
		printf("%12s%30s\n","saved: ",f.c_str());
	}
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. load
		- vec
		- cVec
		- mat
		- cMat
		- spMat
-------------------------------------------------------------------------------------------------------------------------*/

// load vec
void load(const string& f, const SaveOptions& opts, vec& v) {
	uint col = opts.column;
	if (col==0) {
		switch(opts.extras) {
			case SaveOptions::none:		col = 1;
										break;
			case SaveOptions::loc:		col = 2;
										break;
			case SaveOptions::coords:	col = 4;
										break;
			default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
										break;
		}
	}
	uint fileLength = countLines(f);
	uint vLength;
	switch(opts.vectorType) {
		case SaveOptions::simple:	vLength = fileLength;
									break;
		case SaveOptions::real:		vLength = fileLength;
									break;
		case SaveOptions::complex:	vLength = 2*(fileLength-opts.zeroModes);
									break;
		case SaveOptions::realB:	vLength = fileLength;
									break;
		case SaveOptions::complexB:	vLength = 2*(fileLength-opts.zeroModes);
									break;
		default:					cerr << "save error: print vectorType option(" << opts.vectorType << ") not possible" << endl;
									break;
	}	
	vec vf = Eigen::VectorXd::Zero(vLength);
	fstream F;
	F.open(f.c_str(), ios::in);
	string line, temp;
	unsigned int j=0;
	while (getline(F, line)) {
		if (!line.empty()) {
			istringstream ss(line);
			if (col>1) for (unsigned int l=0; l<(col-1); l++) ss >> temp;
			if (opts.vectorType==SaveOptions::complex || opts.vectorType==SaveOptions::complexB) {
				if (j>=(fileLength-opts.zeroModes)) 	ss >> vf(j);
				else 									ss >> vf(2*j) >> vf(2*j+1);
				
			}
			else ss >> vf(j);
			j++;
		}
	}
	F.close();
	
	uint N1 = (opts.paramsIn).N, N2 = (opts.paramsOut).N;
	if (opts.vectorType!=SaveOptions::simple) {
		uint NT1 = (opts.paramsIn).NT, NT2 = (opts.paramsOut).NT;
		uint Nb1 = (opts.paramsIn).Nb, Nb2  = (opts.paramsOut).Nb;
		if (opts.vectorType==SaveOptions::complex && N1>0 && N2>0 && NT1>0 && NT2>0 && (N1!=N2 || NT1!=NT2)) {
			v = interpolate(vf,opts.paramsIn,opts.paramsOut);
		}
		else if (opts.vectorType==SaveOptions::real && N1>0 && N2>0 && NT1>0 && NT2>0 && (N1!=N2 || NT1!=NT2)) {
			v = interpolateReal(vf,opts.paramsIn,opts.paramsOut);
		}
		else if (opts.vectorType==SaveOptions::complexB && N1>0 && N2>0 && Nb1>0 && Nb2>0&& (N1!=N2 || Nb1!=Nb2)) {
			v = interpolate(vf,opts.paramsIn,opts.paramsOut);
		}
		else if (opts.vectorType==SaveOptions::realB && N1>0 && N2>0 && Nb1>0 && Nb2>0 && (N1!=N2 || Nb1!=Nb2)) {
			v = interpolateReal(vf,opts.paramsIn,opts.paramsOut);
		}
		else {
			v = vf;
		}
	}
	else {
		if (N1!=N2 && N2>0) {
			v = interpolate1d(vf,vf.size(),N2);
		}
		else {
			v = vf;
		}
	}
	if (opts.printMessage) {
		printf("%12s%30s\n","loaded: ",f.c_str());
	}
}

// load cVec
void load(const string& f, const SaveOptions& opts, cVec& v) {
	uint col = opts.column;
	if (col==0) {
		switch(opts.extras) {
			case SaveOptions::none:		col = 1;
										break;
			case SaveOptions::loc:		col = 2;
										break;
			case SaveOptions::coords:	col = 4;
										break;
			default:					cerr << "save error: print extras option(" << opts.extras << ") not possible" << endl;
										break;
		}
	}
	uint fileLength = countLines(f);
	uint vLength = fileLength;
	cVec vf = Eigen::VectorXcd::Zero(vLength);
	fstream F;
	F.open(f.c_str(), ios::in);
	string line, temp;
	uint j=0;
	double realPart, imagPart;
	while (getline(F, line)) {
		if (!line.empty()) {
			istringstream ss(line);
			if (col>1) for (unsigned int l=0; l<(col-1); l++) ss >> temp;
			ss >> realPart >> imagPart;
			vf(j) = realPart + comp(0.0,1.0)*imagPart;
			j++;
		}
	}
	F.close();
	
	uint N1 = (opts.paramsIn).N, N2 = (opts.paramsOut).N;
	if (opts.vectorType!=SaveOptions::simple) {
		uint NT1 = (opts.paramsIn).NT, NT2 = (opts.paramsOut).NT;
		uint Nb1 = (opts.paramsIn).Nb, Nb2  = (opts.paramsOut).Nb;
		if (opts.vectorType==SaveOptions::complex || opts.vectorType==SaveOptions::real) {
			if ((N1!=N2 || NT1!=NT2) && N1>0 && N2>0 && NT1>0 && NT2>0) {
				v = interpolate(vf,opts.paramsIn,opts.paramsOut);
			}
			else {
				v = vf;
			}
		}
		else if (opts.vectorType==SaveOptions::complexB || opts.vectorType==SaveOptions::realB) {
			if ((N1!=N2 || Nb1!=Nb2) && N1>0 && N2>0 && Nb1>0 && Nb2>0) {
				v = interpolate(vf,opts.paramsIn,opts.paramsOut);
			}
			else {
				v = vf;
			}
		}
		else {
			v = vf;
		}
	}
	else {
		if (N1!=N2 && N1>0 && N2>0) {
			v = interpolate1d(vf,N1,N2);
		}
		else {
			v = vf;
		}
	}
	if (opts.printMessage) {
		printf("%12s%30s\n","loaded: ",f.c_str());
	}
}

// load mat - assumes square matrix
void load(const string& f, const SaveOptions& opts, mat& m) {
	uint fileLength = countLines(f);
	uint matLength = (uint)sqrt(fileLength);
	fstream F;
	F.open(f.c_str(), ios::in);
	string line;
	uint j, k;
	double v;
	m = Eigen::MatrixXd::Zero(matLength,matLength);
	while (getline(F, line)) {
		if (!line.empty()) {
			istringstream ss(line);
			ss >> j >> k >> v;
			if (j>matLength || k>matLength) cerr << "load error: matrix index > sqrt(fileLength)" << endl;
			m(j,k) = v;
		}
	}
	F.close();
	if (opts.printMessage) {
		printf("%12s%30s\n","loaded: ",f.c_str());
	}
}

// load cMat
void load(const string& f, const SaveOptions& opts, cMat& m) {
	uint fileLength = countLines(f);
	uint matLength = (uint)sqrt(fileLength);
	fstream F;
	F.open(f.c_str(), ios::in);
	string line;
	uint j, k;
	double realPart, imagPart;
	m = Eigen::MatrixXcd::Zero(matLength,matLength);
	while (getline(F, line)) {
		if (!line.empty()) {
			istringstream ss(line);
			ss >> j >> k >> realPart >> imagPart;
			if (j>matLength || k>matLength) cerr << "load error: matrix index > sqrt(fileLength)" << endl;
			m(j,k) = realPart + comp(0.0,1.0)*imagPart;;
		}
	}
	F.close();
	if (opts.printMessage) {
		printf("%12s%30s\n","loaded: ",f.c_str());
	}
}

// load spMat
void load(const string& f , const SaveOptions& opts, spMat& m) {
	uint fileLength = countLines(f);
	Eigen::VectorXi to_reserve(fileLength); //an overestimate
	to_reserve.setZero(fileLength);
	fstream F;
	F.open(f.c_str(), ios::in);
	string line;
	unsigned int nnz = 0, length = 0, count = 0, row, column;
	double value;	
	Eigen::VectorXi rowVec(fileLength), columnVec(fileLength);
	vec valueVec(fileLength);
	while (getline(F, line)) {
		if (!line.empty()) {
			istringstream ss(line);
			ss >> row >> column >> value;
			if (abs(value)>MIN_NUMBER) {
				rowVec(count) = row-1;
				columnVec(count) = column-1;
				valueVec(count) = value;
				to_reserve(row-1) = to_reserve(row-1) + 1;
				nnz++;
				count++;
				if (row>length) {
					length = row;
				}
			}
		}
	}
	if (nnz==0) {
		cerr << "loadSpMat failed, no data in file: " << f << endl;
		}
	to_reserve.conservativeResize(length);
	rowVec.conservativeResize(count);
	columnVec.conservativeResize(count);
	valueVec.conservativeResize(count);
	spMat M(length,length);
	M.setZero(); //just making sure
	M.reserve(to_reserve);
	for (unsigned int l=0;l<count;l++) {
		M.insert(rowVec(l),columnVec(l)) = valueVec(l);
	}
	M.makeCompressed();
	m = M;
	if (opts.printMessage) {
		printf("%12s%30s\n","loaded: ",f.c_str());
	}
}
/*-------------------------------------------------------------------------------------------------------------------------
	4. plot
-------------------------------------------------------------------------------------------------------------------------*/

// plot
void plot(const string& f, const PlotOptions& opts) {
	string style = ((opts.style).empty()? "linespoints": opts.style);
	string output = ((opts.output).empty()? "pics/pic.png": opts.output);
	uint col = (opts.column==0? 1: opts.column);
	string columns = (opts.column2==0? numberToString<uint>(col): (numberToString<uint>(col)+":"+numberToString<uint>(opts.column2)));
	if ((opts.gp).empty()) {		
		string commandOpenStr = "gnuplot -persistent";
		const char * commandOpen = commandOpenStr.c_str();
		FILE * gnuplotPipe = popen (commandOpen,"w");
		string command1Str = "set term png size 1600,800";
		string command2Str = "set output \""+output+"\"";
		string command3Str = "plot \"" + f + "\" using " + columns + " with " + style;
		//string command4Str = "pause -1";
		fprintf(gnuplotPipe, "%s \n",command1Str.c_str());
		fprintf(gnuplotPipe, "%s \n",command2Str.c_str());
		fprintf(gnuplotPipe, "%s \n",command3Str.c_str());
		//fprintf(gnuplotPipe, "%s \n", command4Str.c_str());
		pclose(gnuplotPipe);
	}
	else {
		string commandStr = "gnuplot -e \"inFile='"+f+"'; outFile='"+output+"'; "+" style='"+style+"'\" "+opts.gp;
		FILE * gnuplotPipe = popen(commandStr.c_str(),"w");
		fprintf(gnuplotPipe,"%s\n"," ");
		pclose(gnuplotPipe);
	}
	if (opts.printMessage) {
		printf("%12s%30s\n","plotted: ",output.c_str());
	}
}

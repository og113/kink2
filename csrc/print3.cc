/*-------------------------------------------------------------------------------------------------------------------------
 	definitions for functions to save, load and plot
 -------------------------------------------------------------------------------------------------------------------------*/
 
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense>
#include "folder.h"
#include "simple.h" // for countLines
#include "print.h"

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. save
	2. load
	3. explicit instantiation
	
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. save
-------------------------------------------------------------------------------------------------------------------------*/

// save - saveVectorAscii
template <class T>
void saveVectorAscii(const string& f,  const T& v) {
	ofstream os;
	os.open(f.c_str());
	if (os.good()) {
		os << left << setprecision(16);
		for (uint j=0; j<v.size(); j++) {
			os << setw(25) << v[j] << endl;
		}
		os.close();
	}
	else {
		cerr << "saveVectorAscii error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// save - saveVectorAsciiAppend, appends as column
template <class T>
void saveVectorAsciiAppend(const string& f,  const T& v) {
	uint lengthOs = v.size();
	uint lengthIs = countLines(f);
	ifstream is;
	is.open(f.c_str(),ios::in);
	if (!is.good()) {
		cerr << "saveVectorAsciiAppend error: stream not good for " << f << endl;
		is.close();
		return;
	}
	ofstream os;
	Filename file = f;
	string tempFile = file.Directory+"tempAppend";
	os.open(tempFile.c_str(),ios::out); // may be an error here with multiple nodes due to always sending to same temp file
	if (!os.good()) {
		cerr << "saveVectorAsciiAppend error: stream not good for " << tempFile << endl;
		return;
	}
	os << left << setprecision(16);
	if (lengthOs!=lengthIs) {
		cerr << "saveVectorAsciiAppend error: length of vector, "<< lengthOs\
				 << ", to append not equal to file length, "<< lengthIs << endl;
		return;
	}
	else {
		string lineIn;
		for (uint j=0; j<lengthOs; j++){
			getline(is,lineIn);
			os << lineIn << setw(25) << v[j] << endl;
		}
		is.close();
		os.close();
		copyFile(tempFile,f);
	}
}

// save - saveVectorCsvAppend, appends as row
template <class T>
void saveVectorCsvAppend(const string& f,  const T& v) {
	ofstream os;
	os.open(f.c_str(),ios::app);
	if (os.good()) {
		os << setprecision(16);
		for (uint j=0; j<v.size(); j++) {
			os << v[j];
			if (j<(v.size()-1))
				os << ",";
		}
		os << endl;
		os.close();
	}
	else {
		cerr << "saveVectorCsvAppend error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// save - saveVectorBinary
template <class T>
void saveVectorBinary(const string& f,  const T& v) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	const number* r;
	if (os.good()) {
		for (uint j=0; j<v.size(); j++) {
			r = &v[j];
			os.write(reinterpret_cast<const char*>(r),sizeof(number));
		}
		os.close();
	}
	else {
		cerr << "saveVectorBinary error: cannot write to " << f << endl;
		os.close();
		return;
	}
}

// save - saveVectorBinaryAppend
template <class T>
void saveVectorBinaryAppend(const string& f,  const T& v) {
	ofstream os;
	os.open(f.c_str(),ios::binary | ios::app);
	const number* r;
	if (os.good()) {
		for (uint j=0; j<v.size(); j++) {
			r = &v[j];
			os.write(reinterpret_cast<const char*>(r),sizeof(T));
		}
		os.close();
	}
	else {
		cerr << "saveVectorBinaryAppend error: cannot write to " << f << endl;
		os.close();
		return;
	}
}


// save mat - binary
void saveMatrixBinary(const string& f, const Eigen::MatrixXd& m) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	const double* d;
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			d = &m(j,k);
			os.write(reinterpret_cast<const char*>(d),sizeof(double));
		}
	}
	os.close();
}

// save mat - ascii
void saveMatrixAscii(const string& f, const Eigen::MatrixXd& m) {
	ofstream F;
	F.open(f.c_str());
	if (!F.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	F << left << setprecision(16);
	for (uint j=0; j<m.rows(); j++) {
		for (uint k=0; k<m.cols(); k++) {
			F << setw(24) << m(j,k);
		}
		F << endl;
	}
	F.close();
}

// saveSparseMatrixBinary
/*void saveSparseMatrixBinary(const string& f, const Eigen::SparseMatrix<double>& m) {
	ofstream os;
	os.open(f.c_str(),ios::binary);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	const uint* row, col;
	const double* d;
	for (int l=0; l<m.outerSize(); ++l) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(m,l); it; ++it) {
			row = &(it.row());
			col = &(it.col());
			d = &m(j,k);
			os.write(reinterpret_cast<const char*>(row),sizeof(uint));
			os.write(reinterpret_cast<const char*>(col),sizeof(uint));
			os.write(reinterpret_cast<const char*>(d),sizeof(double));
		}
	}
	os.close();
}*/

// saveSparseMatrixAscii
void saveSparseMatrixAscii(const string& f, const Eigen::SparseMatrix<double>& m) {
	ofstream os;
	os.open(f.c_str(), ios::out);
	if (!os.good()) {
		cerr << "save error: stream not good for " << f << endl;
		return;
	}
	os << left;
	os.precision(16);
	for (int l=0; l<m.outerSize(); ++l) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(m,l); it; ++it) {
			os << setw(25) << it.row() << setw(25) << it.col(); // USED TO HAVE +1 ON it.row() AND it.col()
			os << setw(25) << it.value() << endl;
		}
	}
	os.close();
}


/*-------------------------------------------------------------------------------------------------------------------------
	2. load
-------------------------------------------------------------------------------------------------------------------------*/

// loadVectorAscii
template <class T>
void loadVectorAscii(const string& f, T& v) {
	uint lines = countLines(f);
	v.resize(lines);
	ifstream is;
	is.open(f.c_str());
	if (is.good()) {
		for (uint j=0; j<lines; j++) {
			is >> v[j];
		}
		is.close();
	}
	else {
		cerr << "loadVectorAscii error: cannot read from " << f << endl;
		is.close();
		return;
	}
}

// loadVectorBinary
template <class T>
void loadVectorBinary(const string& f, T& v) {
	number t;
	uint lines = countDoubles(f);
	v.resize(lines);
	ifstream is;
	is.open(f.c_str(),ios::binary);
	if (is.good()) {
		for (uint j=0; j<lines; j++) {
			is.read(reinterpret_cast<char*>(&t),sizeof(number));
			v[j] = t;
		}
		is.close();
	}
	else {
		cerr << "loadVectorBinary error: cannot read from " << f << endl;
		is.close();
		return;
	}
}

// loadVectorAsciiColumn
template <class T>
void loadVectorAsciiColumn(const string& f, T& v, const uint& col) {
	uint cols = countColumns(f);
	if (col>cols) {
		cerr << "loadVectorAsciiColumn error: col(" << col << ")>cols(" << cols << ")" << endl;
		return;
	}
	uint lines = countLines(f);
	v.resize(lines);
	ifstream is;
	is.open(f.c_str());
	if (is.good()) {
		number dross;
		string line;
		for (uint j=0; j<lines; j++) {
			getline(is,line);
			stringstream ss(line);
			for (uint k=0; k<(col-1); k++)
				ss >> dross;
			ss >> v[j];
		}
		is.close();
	}
	else {
		cerr << "loadVectorAsciiColumn error: cannot read from " << f << endl;
		is.close();
		return;
	}
}

// loadVectorCsvAppend
template <class T>
void loadVectorCsvAppend(const string& f,  T& v) {
	uint columns = countColumns(f,',');
	v.resize(columns);
	ifstream is;
	is.open(f.c_str());
	string line;
	if (is.good()) {
	
		is.seekg(0,ios_base::end);      	//Start at end of file
        char ch = ' ';                      //Init ch not equal to '\n'
        while(ch != '\n'){
            is.seekg(-2,ios_base::cur); 	//Two steps back, this means we
                                            //will NOT check the last character
            if((int)is.tellg() <= 0){       //If passed the start of the file,
                is.seekg(0);                //this is the start of the line
                break;
            }
            is.get(ch);                     //Check the next character
        }
        
		getline(is,line);
		stringstream ss(line);
		for (uint k=0; k<columns; k++)
			getline(ss,v[k],',');
			
		is.close();
	}
	else {
		cerr << "loadVectorCsvAppend error: cannot read from " << f << endl;
		is.close();
		return;
	}
}

// loadMatrixAsci - only works for square matrices
void loadMatrixAscii(const string& f, Eigen::MatrixXd& m) {

	uint rowsF = countLines(f), rows;
	uint columnsF = countColumns(f), columns;
	
	if (columnsF==1) {
		rows = (uint)sqrt(rowsF);
		columns = (uint)sqrt(rowsF);
	}
	else {
		rows = rowsF;
		columns = columnsF;
	}
	m = Eigen::MatrixXd::Zero(rows,columns);
	
	fstream is;
	is.open(f.c_str(), ios::in);
	string line;
	uint j = 0, k = 0;
	double v;
	while (getline(is, line)) {
		if (!line.empty()) {
			istringstream ss(line);
			if (columnsF==1) {
				ss >> v;
				m(j,k) = v;
				k++;
				if (k==columns) {
					k = 0;
					j++;
				}
			}
			else {
				while (ss >> v) {
					m(j,k) = v;
					k++;
					if (k==columns) {
						k = 0;
						j++;
					}
				}
			}
		}
	}
	is.close();

}

// load mat  - binary - only works for square matrices
void loadMatrixBinary(const string& f, Eigen::MatrixXd& m) {
	ifstream is;
	is.open(f.c_str(),ios::binary);
	if (!is.good()) {
		cerr << "cannot read from " << f << endl;
	}
	uint pos = is.tellg();
	int lines = -1;
	double d;
	while (!is.eof()) {
		is.read(reinterpret_cast<char*>(&d),sizeof(double));
		lines++;
	}
	if (lines==-1) {
		cerr << "load error: no lines in file " << f << endl;
		return;
	}
	is.clear();
	is.seekg(pos);
	uint rows = (uint)sqrt(lines);
	if (abs((double)rows-sqrt(lines))>MIN_NUMBER*1.0e2) {
		cerr << "load mat error: matrix in " << f << " not square" << endl; 
	}
	m = Eigen::MatrixXd::Zero(rows,rows);
	for (uint j=0; j<rows; j++) {
		for (uint k=0; k<rows; k++) {
			is.read(reinterpret_cast<char*>(&d),sizeof(double));
			m(j,k) = d;
		}
	}
	is.close();
}

// loadSparseMatrixBinary
//void loadSparseMatrixBinary(const string& f, const Eigen::SparseMatrix<double>& m);
// PROBABLY AUGHT TO SAVE AND LOAD FROM TRIPLETS SO THE sizeof(...) BIT WORKS

// loadSparseMatrixAscii
void loadSparseMatrixAscii(const string& f, Eigen::SparseMatrix<double>& m) {
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
				rowVec(count) = row;
				columnVec(count) = column;
				valueVec(count) = value;
				to_reserve(row) = to_reserve(row) + 1;
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
}

/*-------------------------------------------------------------------------------------------------------------------------
	3. explicit instantiation
-------------------------------------------------------------------------------------------------------------------------*/

// save
template void saveVectorBinary< vector<number> >(const string& f, const vector<number>& v);
template void saveVectorBinaryAppend< vector<number> >(const string& f, const vector<number>& v);
template void saveVectorAscii< vector<number> >(const string& f, const vector<number>& v);
template void saveVectorAsciiAppend< vector<number> >(const string& f, const vector<number>& v);
template void saveVectorCsvAppend< vector<number> >(const string& f, const vector<number>& v);
template void saveVectorCsvAppend< vector<string> >(const string& f, const vector<string>& v);
template void saveVectorBinary< Eigen::VectorXd >(const string& f, const Eigen::VectorXd& v);
template void saveVectorBinaryAppend< Eigen::VectorXd >(const string& f, const Eigen::VectorXd& v);
template void saveVectorAscii< Eigen::VectorXd >(const string& f, const Eigen::VectorXd& v);
template void saveVectorAsciiAppend< Eigen::VectorXd >(const string& f, const Eigen::VectorXd& v);
template void saveVectorCsvAppend< Eigen::VectorXd >(const string& f, const Eigen::VectorXd& v);

// load
template void loadVectorBinary< vector<number> >(const string& f, vector<number>& v);
template void loadVectorAscii< vector<number> >(const string& f, vector<number>& v);
template void loadVectorAsciiColumn< vector<number> >(const string& f, vector<number>& v, const uint& col);
template void loadVectorBinary< Eigen::VectorXd >(const string& f, Eigen::VectorXd& v);
template void loadVectorAscii< Eigen::VectorXd >(const string& f, Eigen::VectorXd& v);
template void loadVectorAsciiColumn< Eigen::VectorXd >(const string& f, Eigen::VectorXd& v, const uint& col);
template void loadVectorCsvAppend< vector<string> >(const string& f, vector<string>& v);

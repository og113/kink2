/*-------------------------------------------------------------------------------------------------------------------------
 	declarations for functions to save, load and plot
 -------------------------------------------------------------------------------------------------------------------------*/
  
#ifndef __PRINT3_H_INCLUDED__
#define __PRINT3_H_INCLUDED__

#include "simple.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;

/*-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------
CONTENTS
	1. save
	2. load
	
n.b. for these functions to work, as templates, the class T be a vector of numbers and must have indexing via [], .size() and .resize(uint) functions.
-------------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------------------------------------------
	1. save
-------------------------------------------------------------------------------------------------------------------------*/

// save - saveVectorAscii
template <class T>
void saveVectorAscii(const string& f,  const T& v);

// save - saveVectorAsciiAppend
template <class T>
void saveVectorAsciiAppend(const string& f,  const T& v);

// save - saveVectorCsvAppend
template <class T>
void saveVectorCsvAppend(const string& f,  const T& v);

// save - saveVectorBinary
template <class T>
void saveVectorBinary(const string& f,  const T& v);

// save - saveVectorBinaryAppend
template <class T>
void saveVectorBinaryAppend(const string& f,  const T& v);

// save - saveMatrixBinary
void saveMatrixBinary(const string& f, const Eigen::MatrixXd& m);

// saveMatrixAscii
void saveMatrixAscii(const string& f, const Eigen::MatrixXd& m);

// saveSparseMatrixBinary
//void saveSparseMatrixAscii(const string& f, const Eigen::SparseMatrix<double>& m);

// saveSparseMatrixAscii
void saveSparseMatrixAscii(const string& f, const Eigen::SparseMatrix<double>& m);

/*-------------------------------------------------------------------------------------------------------------------------
	2. load
-------------------------------------------------------------------------------------------------------------------------*/

// loadVectorBinary
template <class T>
void loadVectorBinary(const string& f, T& v);

// loadVectorAscii
template <class T>
void loadVectorAscii(const string& f, T& v);

// loadVectorAsciiColumn
template <class T>
void loadVectorAsciiColumn(const string& f, T& v, const uint& col);

// loadVectorCsvAppend
template <class T>
void loadVectorCsvAppend(const string& f,  T& v);

// loadMatrixAscii
void loadMatrixAscii(const string& f, Eigen::MatrixXd& m);

// loadMatrixBinary
void loadMatrixBinary(const string& f, Eigen::MatrixXd& m);

// loadSparseMatrixBinary
//void loadSparseMatrixAscii(const string& f, Eigen::SparseMatrix<double>& m);

// loadSparseMatrixAscii
void loadSparseMatrixAscii(const string& f, Eigen::SparseMatrix<double>& m);

#endif // __PRINT3_H_INCLUDED__

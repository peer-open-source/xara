/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Description: This file contains the class definition for Matrix.
// Matrix is a concrete class implementing the matrix abstraction.
// Matrix class is used to provide the abstraction for the most
// general type of matrix, that of an unsymmetric full matrix.
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#ifndef Matrix_h
#define Matrix_h 
#define NO_STATIC_WORK
#include <assert.h>
#include <cstddef>
using std::size_t;

class Vector;
class ID;
class Message;
class OPS_Stream;
namespace OpenSees {
  template<int n, typename T> struct VectorND;
  template<int, int, typename T> struct MatrixND;
};
// struct Vector3D;
using Vector3D = OpenSees::VectorND<3,double>;

class Matrix
{
  public:
    // constructors and destructor
    Matrix();	
    Matrix(int nrows, int ncols);
    Matrix(double *data, int nrows, int ncols);    
    Matrix(const Matrix &M);
#if !defined(NO_CXX11_MOVE)
    Matrix( Matrix &&M);    
#endif
    ~Matrix();

    // utility methods
    int setData(double *newData, int nRows, int nCols);
    template <int nr, int nc>
    inline int setData(OpenSees::MatrixND<nr,nc,double> &M) {
      if (!fromFree && data != nullptr)
        delete[] data;

      fromFree = 1; // Cannot delete data
      data = &M.values[0][0];
      numRows = nr;
      numCols = nc;
      dataSize= nr*nc;
      return 0;
    }

    inline int noRows() const;
    inline int noCols() const;

    void Zero();
    int resize(int numRow, int numCol);
    Vector diagonal() const;

    int  Assemble(const Matrix &,const ID &rows, const ID &cols, double fact = 1.0);  
    // methods added by Remo
    int  Assemble(const Matrix &V, int init_row, int init_col, double fact = 1.0);

    template <int nr, int nc> int  Assemble(const OpenSees::MatrixND<nr,nc,double> &V, int init_row, int init_col, double fact = 1.0);
    int  Assemble(const Vector &V, int init_row, int init_col, double fact = 1.0);
    int  AssembleTranspose(const Matrix &V, int init_row, int init_col, double fact = 1.0);
    int  AssembleTranspose(const Vector &V, int init_row, int init_col, double fact = 1.0);
    int  Extract(const Matrix &V, int init_row, int init_col, double fact = 1.0);

    
    int Solve(const Vector &V, Vector &res) const;
    int Solve(const Vector &V, Vector &res);
    int Solve(const Matrix &M, Matrix &res);
    int Invert(Matrix &res) const;
    int Invert();

    int addMatrix(const Matrix &other, double factOther);
    int addMatrix(double factThis, const Matrix &other, double factOther);
    int addMatrixTranspose(double factThis, const Matrix &other, double factOther);
    int addMatrixProduct(double factThis, const Matrix &A, const Matrix &B, double factOther); // AB
    int addMatrixTransposeProduct(double factThis, const Matrix &A, const Matrix &B, double factOther); // A'B
    int addMatrixTripleProduct(double factThis, const Matrix &A, const Matrix &B, double factOther); // A'BA
    int addMatrixTripleProduct(double factThis, const Matrix &A, const Matrix &B, const Matrix &C, double otherFact); //A'BC
    
    //
    // Inline operations
    //

    template<int NR, int NC>
    void addMatrix(const OpenSees::MatrixND<NR, NC, double>& M, double fact);

    template<class VecT>
    void addTensorProduct(const VecT& V, const VecT& W);
    template<class VecT>
    void addTensorProduct(const VecT& V, const VecT& W, double factThis);

    template<class VecT> void addSpin(const VecT& V);
    template<class VecT> void addSpin(const VecT& V, double mult);
    template<class VecT>
    void addSpinAtRow(const VecT& V, size_t row_index);
    template<class VecT>
    void addSpinAtRow(const VecT& V, size_t vector_index, size_t matrix_row_index);
    template<class VecT>
    void addSpinAtRow(const VecT& V, double mult, size_t row_index);
    template<class VecT>
    void addSpinAtRow(const VecT& V, double mult, size_t vector_index, size_t matrix_row_index);

    // overloaded operators 
    inline double &operator()(int row, int col);
    inline double  operator()(int row, int col) const;
    Matrix         operator()(const ID &rows, const ID & cols) const; 
    Matrix        &operator=(const Matrix &M);
    Matrix        &operator=(Matrix &&M);

    //
    // Matrix operations which will preserve the derived type and
    // which can be implemented efficiently without many constructor calls.
    //

    // Matrix-scalar operations
    Matrix &operator+=(double fact);
    Matrix &operator-=(double fact);
    Matrix &operator*=(double fact);
    Matrix &operator/=(double fact); 

    // Matrix operations which generate a new Matrix. They are not the
    // most efficient to use, as constructors must be called twice. They
    // however are usefull for matlab like expressions involving Matrices.

    // Matrix-scalar operations
    Matrix operator+(double fact) const;
    Matrix operator-(double fact) const;
    Matrix operator*(double fact) const;
    Matrix operator/(double fact) const;
    
    // Matrix-vector operations
    Vector operator^(const Vector &V) const;
    Vector operator*(const Vector &V) const;
    Vector operator*(const Vector3D &V) const;

    
    // matrix-matrix operations
    Matrix operator+(const Matrix &M) const;
    Matrix operator-(const Matrix &M) const;
    Matrix operator*(const Matrix &M) const;
//     Matrix operator/(const Matrix &M) const;    
    Matrix operator^(const Matrix &M) const;
    Matrix &operator+=(const Matrix &M);
    Matrix &operator-=(const Matrix &M);

    // methods to read/write to/from the matrix
    void Output(OPS_Stream &s) const;
    //    void Input(istream &s);

    friend OPS_Stream &operator<<(OPS_Stream &s, const Matrix &M);
    //    friend istream &operator>>(istream &s, Matrix &M);    
    friend Matrix operator*(double a, const Matrix &M);
    
    
    friend class Vector;
    template<int n, typename> friend struct OpenSees::VectorND;
    friend class Message;
    friend class UDP_Socket;
    friend class TCP_Socket;
    friend class TCP_SocketSSL;
    friend class TCP_SocketNoDelay;
    friend class MPI_Channel;
    friend class MySqlDatastore;
    friend class BerkeleyDbDatastore;

  protected:

  private:
    static double MATRIX_NOT_VALID_ENTRY;
#ifdef NO_STATIC_WORK
    double *matrixWork = nullptr;
    int *intWork = nullptr;
    int sizeDoubleWork = 400;
    int sizeIntWork = 20;
#else
    static double *matrixWork;
    static int *intWork;
    static int sizeDoubleWork;
    static int sizeIntWork;
#endif

    int numRows;
    int numCols;
    int dataSize;
    double *data;
    int fromFree;
};


/********* INLINED MATRIX FUNCTIONS ***********/
inline int 
Matrix::noRows() const 
{
  return numRows;
}

inline int 
Matrix::noCols() const 
{
  return numCols;
}


inline double &
Matrix::operator()(int row, int col)
{ 
#ifdef _G3DEBUG
  if ((row < 0) || (row >= numRows)) {
    opserr << "Matrix::operator() - row " << row << " our of range [0, " <<  numRows-1 << endln;
    return data[0];
  } else if ((col < 0) || (col >= numCols)) {
    opserr << "Matrix::operator() - row " << col << " our of range [0, " <<  numCols-1 << endln;
    return MATRIX_NOT_VALID_ENTRY;
  }
#endif
  return data[col*numRows + row];
}


inline double 
Matrix::operator()(int row, int col) const
{ 
#ifdef _G3DEBUG
  if ((row < 0) || (row >= numRows)) {
    opserr << "Matrix::operator() - row " << row << " our of range [0, " <<  numRows-1 << endln;
    return data[0];
  } else if ((col < 0) || (col >= numCols)) {
    opserr << "Matrix::operator() - row " << col << " our of range [0, " <<  numCols-1 << endln;
    return MATRIX_NOT_VALID_ENTRY;
  }
#endif
  return data[col*numRows + row];
}

template<int NR, int NC>
void Matrix::addMatrix(const OpenSees::MatrixND<NR, NC, double>& M, double fact)
{
  for (int i = 0; i< NR; i++)
    for (int j = 0; j< NC; j++)
      (*this)(i,j) += fact*M(i,j);
}

template <int nr, int nc>
int
Matrix::Assemble(const OpenSees::MatrixND<nr,nc,double> &M, int init_row, int init_col, double fact) 
{
  
  [[maybe_unused]] int final_row = init_row + nr - 1;
  [[maybe_unused]] int final_col = init_col + nc - 1;

  assert((init_row >= 0) && (final_row < numRows) && (init_col >= 0) && (final_col < numCols));

  for (int i=0; i<nc; i++) {
     const int pos_Cols = init_col + i;
     for (int j=0; j<nr; j++) {
        const int pos_Rows = init_row + j; 
        (*this)(pos_Rows,pos_Cols) += M(j,i)*fact;
     }
  }  

  return 0;
}

/**
* Computes the Spin of the input vector V, and saves the result into the output matrix S,
* at the specified row index.
* Note: no check is made on the size of the input-output arguments.
* @param V the input vector (assumed size: >= 3)
* @param S the output matrix (assumed size: >= 3x3)
* @param row_index the index of the first row in the output matrix where the spin has to be saved
*/
template< class TVec>
inline void Matrix::addSpinAtRow(const TVec& V, size_t row_index)
{
    size_t i0 = row_index;
    size_t i1 = 1 + row_index;
    size_t i2 = 2 + row_index;
    double v0 = V(i0);
    double v1 = V(i1);
    double v2 = V(i2);
                              (*this)(i0, 1) += -v2;     (*this)(i0, 2) +=  v1;
    (*this)(i1, 0) +=  v2;                               (*this)(i1, 2) += -v0;
    (*this)(i2, 0) += -v1;    (*this)(i2, 1) +=  v0;                          ;
}

/**
* Computes the Spin of the input vector V, from the specified index, and saves the result into the output matrix S,
* at the specified row index.
* Note: no check is made on the size of the input-output arguments.
* @param V the input vector (assumed size: >= 3)
* @param S the output matrix (assumed size: >= 3x3)
* @param vector_index the index of the first component of the input vector to be used to compute the spin
* @param row_index the index of the first row in the output matrix where the spin has to be saved
*/
template< class TVec>
inline void Matrix::addSpinAtRow(const TVec& V, size_t vector_index, size_t matrix_row_index)
{
    size_t i0 = matrix_row_index;
    size_t i1 = 1 + matrix_row_index;
    size_t i2 = 2 + matrix_row_index;
    double v0 = V(vector_index);
    double v1 = V(vector_index + 1);
    double v2 = V(vector_index + 2);
                              (*this)(i0, 1) += -v2;    (*this)(i0, 2) +=  v1;
    (*this)(i1, 0) +=  v2;                              (*this)(i1, 2) += -v0;
    (*this)(i2, 0) += -v1;    (*this)(i2, 1) +=  v0;
}

/**
* Computes the Spin of the input vector V, and saves the result into the output matrix S.
* This version uses a multiplier for the output values.
* Note: no check is made on the size of the input-output arguments.
* @param V the input vector (assumed size: >= 3)
* @param S the output matrix (assumed size: >= 3x3)
* @param mult the multiplier for the output values
*/
template< class VecT>
inline void Matrix::addSpin(const VecT& v, double mult)
{
   const double v0 = mult*v[0],
                v1 = mult*v[1],
                v2 = mult*v[2];

                          (*this)(0, 1) -=  v2;     (*this)(0, 2) +=  v1;
  (*this)(1, 0) +=  v2;                             (*this)(1, 2) -=  v0;
  (*this)(2, 0) -=  v1;   (*this)(2, 1) +=  v0;
}
/**
* Computes the Spin of the input vector V, and saves the result into the output matrix S.
* Note: no check is made on the size of the input-output arguments.
* @param V the input vector (assumed size: >= 3)
* @param S the output matrix (assumed size: >= 3x3)
*/
template< class TVec>
inline void Matrix::addSpin(const TVec& V)
{
                                  (*this)(0, 1) -=  V[2];       (*this)(0, 2) +=  V[1];
    (*this)(1, 0) += V[2];                                      (*this)(1, 2) -=  V[0];
    (*this)(2, 0) -= V[1];        (*this)(2, 1) +=  V[0];
}


/**
* Computes the Spin of the input vector V, and saves the result into the output matrix S,
* at the specified row index.
* This version uses a multiplier for the output values.
* Note: no check is made on the size of the input-output arguments.
* @param V the input vector (assumed size: >= 3)
* @param S the output matrix (assumed size: >= 3x3)
* @param mult the multiplier for the output values
* @param row_index the index of the first row in the output matrix where the spin has to be saved
*/
template< class TVec>
inline void Matrix::addSpinAtRow(const TVec& V, double mult, size_t row_index)
{
    size_t i0 = row_index;
    size_t i1 = 1 + row_index;
    size_t i2 = 2 + row_index;
    const double v0 = mult * V(i0);
    const double v1 = mult * V(i1);
    const double v2 = mult * V(i2);
    (*this)(i0, 0) = 0.0;    (*this)(i0, 1) = -v2;    (*this)(i0, 2) =  v1;
    (*this)(i1, 0) =  v2;    (*this)(i1, 1) = 0.0;    (*this)(i1, 2) = -v0;
    (*this)(i2, 0) = -v1;    (*this)(i2, 1) =  v0;    (*this)(i2, 2) = 0.0;
}

/**
* Computes the Spin of the input vector V, from the specified index, and saves the result into the output matrix S,
* at the specified row index.
* This version uses a multiplier for the output values.
* Note: no check is made on the size of the input-output arguments.
* @param V the input vector (assumed size: >= 3)
* @param S the output matrix (assumed size: >= 3x3)
* @param mult the multiplier for the output values
* @param vector_index the index of the first component of the input vector to be used to compute the spin
* @param row_index the index of the first row in the output matrix where the spin has to be saved
*/
template< class TVec>
inline void Matrix::addSpinAtRow(const TVec& V, double mult, size_t vector_index, size_t matrix_row_index)
{
    size_t i0 = matrix_row_index;
    size_t i1 = 1 + matrix_row_index;
    size_t i2 = 2 + matrix_row_index;
    const double v0 = mult * V(vector_index);
    const double v1 = mult * V(vector_index + 1);
    const double v2 = mult * V(vector_index + 2);

                               (*this)(i0, 1) += -v2;    (*this)(i0, 2) += v1;
    (*this)(i1, 0) +=  v2;                               (*this)(i1, 2) += -v0;
    (*this)(i2, 0) += -v1;     (*this)(i2, 1) +=  v0;
}

#endif


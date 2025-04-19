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
// Description: This file contains the class implementation for Matrix.
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#include "Matrix.h"
#include "Vector.h"
#include "Vector3D.h"
#include "routines/cmx.h"
#include "blasdecl.h"
#include "ID.h"

#include <stdlib.h>
#include <stdexcept>
#include <OPS_Stream.h>

#define MATRIX_VERY_LARGE_VALUE 1.0e213

#include <math.h>
#include <assert.h>

#ifndef NO_STATIC_WORK
# define MATRIX_WORK_AREA 400
# define INT_WORK_AREA 20
  int Matrix::sizeDoubleWork = MATRIX_WORK_AREA;
  int Matrix::sizeIntWork = INT_WORK_AREA;
  double *Matrix::matrixWork = nullptr;
  int    *Matrix::intWork    = nullptr;
#endif

//#define MATRIX_BLAS
//#define NO_WORK


//
// CONSTRUCTORS
//

Matrix::Matrix()
:numRows(0), numCols(0), dataSize(0), data(0), fromFree(0)
{
  // allocate work areas if the first
  if (matrixWork == nullptr) {
    matrixWork = new double[sizeDoubleWork];
    intWork    = new int[sizeIntWork];
  }
}


Matrix::Matrix(int nRows,int nCols)
:numRows(nRows), numCols(nCols), dataSize(0), data(0), fromFree(0)
{
//assert(nRows > 0);
//assert(nCols > 0);

  // allocate work areas if the first matrix
  if (matrixWork == nullptr) {
    matrixWork = new double[sizeDoubleWork];
    intWork    = new int[sizeIntWork];
  }

  dataSize = numRows * numCols;
  data = nullptr;

  if (dataSize > 0)
    data = new double[dataSize]{};
}

Matrix::Matrix(double *theData, int row, int col) 
:numRows(row),numCols(col),dataSize(row*col),data(theData),fromFree(1)
{
//assert(row > 0);
//assert(col > 0);

  // allocate work areas if the first matrix
  if (matrixWork == nullptr) {
    matrixWork = new double[sizeDoubleWork];
    intWork    = new int[sizeIntWork];
  }

}


Matrix::Matrix(const Matrix &other)
: numRows(0), numCols(0), dataSize(0), data(0), fromFree(0)
{
  // allocate work areas if the first matrix
  if (matrixWork == nullptr) {
    matrixWork = new double[sizeDoubleWork];
    intWork    = new int[sizeIntWork];
  }

  numRows  = other.numRows;
  numCols  = other.numCols;
  dataSize = other.dataSize;

  if (dataSize != 0) {
    data = new double[dataSize];
    // copy the data
    double *dataPtr = data;
    double *otherDataPtr = other.data;
    for (int i=0; i<dataSize; i++)
      *dataPtr++ = *otherDataPtr++;
  }
}

// Move constructor
#if !defined(NO_CXX11_MOVE)
Matrix::Matrix(Matrix &&other)
: numRows(other.numRows), numCols(other.numCols), 
  dataSize(other.dataSize), data(other.data),
  fromFree(other.fromFree)
{
  other.numRows  = 0;
  other.numCols  = 0;
  other.dataSize = 0;
  other.data     = nullptr;
  other.fromFree = 1;
  throw std::runtime_error("error");
}
#endif

//
// DESTRUCTOR
//

Matrix::~Matrix()
{
  if (data != nullptr) {
    if (fromFree == 0 && dataSize > 0){
      delete [] data; 
      data = nullptr;
    }
  }
#ifdef NO_STATIC_WORK
  if (matrixWork != nullptr)
    delete [] matrixWork;

  if (intWork != nullptr)
    delete [] intWork;
#endif
}
    

//
// METHODS - Zero, Assemble, Solve
//
int
Matrix::setData(double *theData, int row, int col) 
{
//assert(row > 0);
//assert(col > 0);

  // delete the old if allocated
  if (data != nullptr)
    if (fromFree == 0) {
      delete [] data; 
      data = 0;
    }
  numRows  = row;
  numCols  = col;
  dataSize = row*col;
  data     = theData;
  fromFree = 1;
  
  return 0;
}

void
Matrix::Zero(void)
{
  double *dataPtr = data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ = 0;
}


int
Matrix::resize(int rows, int cols) {

  int newSize = rows*cols;
  assert(newSize >= 0);

  if (newSize > dataSize) {

    // free the old space
    if (data != nullptr)
      if (fromFree == 0){
        delete [] data; 
        data = 0;
      }

    fromFree = 0;
    // create new space
    data = new double[newSize];
    dataSize = newSize;
    numRows = rows;
    numCols = cols;
  }

  // just reset the cols and rows - save two memory calls at expense of holding 
  // onto extra memory
  else {
    numRows = rows;
    numCols = cols;
  }

  return 0;
}


int
Matrix::Assemble(const Matrix &V, const ID &rows, const ID &cols, double fact) 
{
  int res = 0;
  int pos_Rows, pos_Cols;
  for (int i=0; i<cols.Size(); i++) {
    pos_Cols = cols(i);

    for (int j=0; j<rows.Size(); j++) {
      pos_Rows = rows(j);
      
      assert((pos_Cols >= 0)      && (pos_Rows >= 0) && (pos_Rows < numRows) &&
             (pos_Cols < numCols) && (i < V.numCols) && (j < V.numRows));
      (*this)(pos_Rows,pos_Cols) += V(j,i)*fact;
    }
  }

  return res;
}

int
Matrix::Solve(const Vector &b, Vector &x)
{
    int n = numRows;
    assert(numRows == numCols);
    assert(numRows == x.Size());
    assert(numRows == b.Size());

    // check work area can hold all the data
    if (dataSize > sizeDoubleWork) {
      if (matrixWork != nullptr) {
        delete [] matrixWork;
        matrixWork = nullptr;
      }
      matrixWork = new double[dataSize];
      sizeDoubleWork = dataSize;
    }

    // check work area can hold all the data
    if (n > sizeIntWork) {

      if (intWork != nullptr) {
        delete [] intWork;
        intWork = nullptr;
      }
      intWork = new int[n];
      sizeIntWork = n;  
    }
 
    // copy the data
    for (int i=0; i<dataSize; i++)
      matrixWork[i] = data[i];

    // set x equal to b
    x = b;

    int nrhs = 1;
    int ldA = n;
    int ldB = n;
    int info;
    double *Aptr = matrixWork;
    double *Xptr = x.theData;
    int *iPIV = intWork;

    DGESV(&n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);

    return -abs(info);
}

int
Matrix::Solve(const Vector &b, Vector &x) const
{
    int n = numRows;
    assert(numRows == numCols);
    assert(numRows == x.Size());
    assert(numRows == b.Size());

#if defined(NO_STATIC_WORK)
    /* static */ double *matrixWork = nullptr;
    /* static */ int *intWork = nullptr;
    /* static */ int sizeDoubleWork = 0;
    /* static */ int sizeIntWork    = 0;
#endif

    // check work area can hold all the data
    if (dataSize > sizeDoubleWork) {
      if (matrixWork != nullptr) {
        delete [] matrixWork;
        matrixWork = nullptr;
      }
      matrixWork = new double[dataSize];
      sizeDoubleWork = dataSize;
    }
    if (n > sizeIntWork) {
      if (intWork != nullptr) {
        delete [] intWork;
        intWork = nullptr;
      }
      intWork = new int[n];
      sizeIntWork = n;  
    }
 
    // copy the data
    int i;
    for (i=0; i<dataSize; i++)
      matrixWork[i] = data[i];

    // set x equal to b
    x = b;

    int nrhs = 1;
    int ldA = n;
    int ldB = n;
    int info;
    double *Aptr = matrixWork;
    double *Xptr = x.theData;
    int *iPIV = intWork;

    DGESV(&n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);

    delete [] intWork;
    delete [] matrixWork;

    return -abs(info);
}


int
Matrix::Solve(const Matrix &b, Matrix &x) // const
{

    int n = numRows;
    int nrhs = x.numCols;
    assert(numRows == numCols);
    assert(n == x.numRows);
    assert(n == b.numRows);
    assert(x.numCols == b.numCols);

    // check work area can hold all the data
    if (dataSize > sizeDoubleWork) {
      if (matrixWork != 0) {
        delete [] matrixWork;
        matrixWork = 0;
      }
      matrixWork = new double[dataSize];
      sizeDoubleWork = dataSize;
    }
    if (n > sizeIntWork) {
      if (intWork != nullptr) {
        delete [] intWork;
        intWork = 0;
      }
      intWork = new int[n];
      sizeIntWork = n;
    }

    // copy the data
    int i;
    for (i=0; i<dataSize; i++)
      matrixWork[i] = data[i];

    x = b;

    int ldA = n;
    int ldB = n;
    int info;
    double *Aptr = matrixWork;
    double *Xptr = x.data;
    
    int *iPIV = intWork;
    
    info = -1;

#ifdef _WIN32
    DGESV(&n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
#else
    dgesv_(&n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);

    /*
    // further correction if required
    double Bptr[n*n];
    for (int i=0; i<n*n; i++) Bptr[i] = b.data[i];
    double *origData = data;
    double Ferr[n];
    double Berr[n];
    double newWork[3*n];
    int newIwork[n];
    
    dgerfs_("N",&n,&n,origData,&ldA,Aptr,&n,iPIV,Bptr,&ldB,Xptr,&ldB,
            Ferr, Berr, newWork, newIwork, &info);
    */
#endif
    return -abs(info);
}

int
Matrix::Invert(Matrix &theInverse) const
{
    assert(numRows == numCols);
    assert(numRows == theInverse.numRows);
    theInverse = *this;
    return theInverse.Invert();
}

int
Matrix::Invert()
{

  int info;
  int n = numRows;
  assert(numRows == numCols);
  switch (numRows) {
    case 2:
      cmx_inv2(data, data, &info);
      break;
    case 3:
      cmx_inv3(data, data, &info);
      break;
    case 4:
      cmx_inv4(data, data, &info);
      break;
    case 5:
      cmx_inv5(data, data, &info);
      break;
    case 6:
      cmx_inv6(data, data, &info);
      break;

    default:
  
      // check work area can hold all the data
      if (dataSize > sizeDoubleWork) {
        if (matrixWork != nullptr) {
          delete [] matrixWork;
          matrixWork = 0;
        }
        matrixWork = new double[dataSize];
        sizeDoubleWork = dataSize;
      }

      // check work area can hold all the data
      if (n > sizeIntWork) {
        if (intWork != 0) {
          delete [] intWork;
          intWork = nullptr;
        }
        intWork = new int[n];
        sizeIntWork = n;  
      }
#if 0
      // copy the data 
      for (int i=0; i<dataSize; i++)
        matrixWork[i] = data[i];
#endif
      // theInverse = *this;
      int ldA = n;
      double *Wptr = matrixWork;
      double *Aptr = data;
      int workSize = sizeDoubleWork;
      
      int *iPIV = intWork;

      DGETRF(&n,&n,Aptr,&ldA,iPIV,&info);
      if (info != 0) 
        return -abs(info);
      DGETRI(&n,Aptr,&ldA,iPIV,Wptr,&workSize,&info);

  }
  return -abs(info);
}

int
Matrix::addMatrix(const Matrix &other, double factOther)
{
    assert(other.numRows == numRows);
    assert(other.numCols == numCols);

    if (factOther == 0.0)
      return 0;


    // want: this += other * factOther
    if (factOther == 1.0) {
      double *dataPtr = data;
      double *otherDataPtr = other.data;                    
      for (int i=0; i<dataSize; i++)
        *dataPtr++ += *otherDataPtr++;
    } else {
      double *dataPtr = data;
      double *otherDataPtr = other.data;                    
      for (int i=0; i<dataSize; i++)
        *dataPtr++ += *otherDataPtr++ * factOther;
    }
    // successfull
    return 0;
}

#include <stdexcept>

int
Matrix::addMatrix(double factThis, const Matrix &other, double factOther)
{
    assert(other.numRows == numRows);
    assert(other.numCols == numCols);

    if (factThis == 1.0 && factOther == 0.0)
      return 0;

    if (other.data == nullptr)
      throw std::runtime_error("error");

    if (factThis == 1.0) {
      // want: this += other * factOther
      if (factOther == 1.0) {
        double *dataPtr = data;
        double *otherDataPtr = other.data;                    
        for (int i=0; i<dataSize; i++)
          *dataPtr++ += *otherDataPtr++;
      } else {
        double *dataPtr = data;
        double *otherDataPtr = other.data;                    
        for (int i=0; i<dataSize; i++)
          *dataPtr++ += *otherDataPtr++ * factOther;
      }
    }
    else if (factThis == 0.0) {
      // want: this = other * factOther
      if (factOther == 1.0) {
        double *dataPtr = data;
        double *otherDataPtr = other.data;                    
        for (int i=0; i<dataSize; i++)
          *dataPtr++ = *otherDataPtr++;
      } else {
        double *dataPtr = data;
        double *otherDataPtr = other.data;                    
        for (int i=0; i<dataSize; i++)
          *dataPtr++ = *otherDataPtr++ * factOther;
      }
    }
    else {
      // want: this = this * thisFact + other * factOther
      if (factOther == 1.0) {
        double *dataPtr = data;
        double *otherDataPtr = other.data;                    
        for (int i=0; i<dataSize; i++) {
          double value = *dataPtr * factThis + *otherDataPtr++;
          *dataPtr++ = value;
        }
      } else {
        double *dataPtr = data;
        double *otherDataPtr = other.data;                    
        for (int i=0; i<dataSize; i++) {
          double value = *dataPtr * factThis + *otherDataPtr++ * factOther;
          *dataPtr++ = value;
        }
      }
    } 

    // successfull
    return 0;
}


int
Matrix::addMatrixTranspose(double factThis, const Matrix &other, double factOther)
{
    assert(other.numRows == numCols);
    assert(other.numCols == numRows);

    if (factThis == 1.0 && factOther == 0.0)
      return 0;

    if (factThis == 1.0) {

        // want: this += other^T * factOther
      if (factOther == 1.0) {
        double *dataPtr = data;
        for (int j=0; j<numCols; j++) {
          for (int i=0; i<numRows; i++)
                *dataPtr++ += (other.data)[j+i*numCols];
        }
      } else {
          double *dataPtr = data;
          for (int j=0; j<numCols; j++) {
            for (int i=0; i<numRows; i++)
                  *dataPtr++ += (other.data)[j+i*numCols] * factOther;
          }
        }
      } 

    else if (factThis == 0.0) {

      // want: this = other^T * factOther
      if (factOther == 1.0) {
        double *dataPtr = data;
        for (int j=0; j<numCols; j++) {
          for (int i=0; i<numRows; i++)
                *dataPtr++ = (other.data)[j+i*numCols];
        }
      } else {
        double *dataPtr = data;
        for (int j=0; j<numCols; j++) {
          for (int i=0; i<numRows; i++)
                *dataPtr++ = (other.data)[j+i*numCols] * factOther;
        }
      }
    } 

    else {

      // want: this = this * thisFact + other^T * factOther
      if (factOther == 1.0) {
        double *dataPtr = data;
        for (int j=0; j<numCols; j++) {
          for (int i=0; i<numRows; i++) {
            double value = *dataPtr * factThis + (other.data)[j+i*numCols];
                *dataPtr++ = value;
          }
        }
      } else {
        double *dataPtr = data;
        for (int j=0; j<numCols; j++) {
          for (int i=0; i<numRows; i++) {
                double value = *dataPtr * factThis + (other.data)[j+i*numCols] * factOther;
                *dataPtr++ = value;
          }
        }
      }
    } 

    // successfull
    return 0;
}


int
Matrix::addMatrixProduct(double thisFact, 
                         const Matrix &B, 
                         const Matrix &C, 
                         double otherFact)
{
    assert(B.numRows == numRows);
    assert(C.numCols == numCols);
    assert(B.numCols == C.numRows);

    if (thisFact == 1.0 && otherFact == 0.0)
      return 0;

#ifdef MATRIX_BLAS
    else if (numRows >  6) {
      int m = numRows,
          n = C.numCols,
          k = C.numRows;
      DGEMM("N", "N", &m, &n, &k,&otherFact, B.data, &  m,
                                             C.data, &  k,
                                 &thisFact,    data, &  m);
      return 0;
    }
#endif

    // NOTE: looping as per blas3 dgemm_: j,k,i
    else if (thisFact == 1.0) {

      // want: this += B * C  otherFact
      int numColB = B.numCols;
      double *ckjPtr  = &(C.data)[0];
      for (int j=0; j<numCols; j++) {
        double *aijPtrA = &data[j*numRows];
        for (int k=0; k<numColB; k++) {
          double tmp = *ckjPtr++ * otherFact;
          double *aijPtr = aijPtrA;
          double *bikPtr = &(B.data)[k*numRows];
          for (int i=0; i<numRows; i++)
            *aijPtr++ += *bikPtr++ * tmp;
        }
      }
    }

    else if (thisFact == 0.0) {
      // want: this = B * C  otherFact
      double *dataPtr = data;
      for (int i=0; i<dataSize; i++)
          *dataPtr++ = 0.0;
      int numColB = B.numCols;
      double *ckjPtr  = &(C.data)[0];
      for (int j=0; j<numCols; j++) {
        double *aijPtrA = &data[j*numRows];
        for (int k=0; k<numColB; k++) {
          double tmp = *ckjPtr++ * otherFact;
          double *aijPtr = aijPtrA;
          double *bikPtr = &(B.data)[k*numRows];
          for (int i=0; i<numRows; i++)
            *aijPtr++ += *bikPtr++ * tmp;
        }
      }
    } 

    else {
      // want: this = B * C  otherFact
      double *dataPtr = data;
      for (int i=0; i<dataSize; i++)
          *dataPtr++ *= thisFact;
      int numColB = B.numCols;
      double *ckjPtr  = &(C.data)[0];
      for (int j=0; j<numCols; j++) {
        double *aijPtrA = &data[j*numRows];
        for (int k=0; k<numColB; k++) {
          double tmp = *ckjPtr++ * otherFact;
          double *aijPtr = aijPtrA;
          double *bikPtr = &(B.data)[k*numRows];
          for (int i=0; i<numRows; i++)
            *aijPtr++ += *bikPtr++ * tmp;
        }
      }
    } 
    return 0;
}

int
Matrix::addMatrixTransposeProduct(double thisFact, 
                                  const Matrix &B, 
                                  const Matrix &C, 
                                  double otherFact)
{
  assert((B.numCols == numRows) && 
         (C.numCols == numCols) && 
         (B.numRows == C.numRows));

  if (thisFact == 1.0 && otherFact == 0.0)
    return 0;

#ifdef MATRIX_BLAS
  else if (numRows >  6) {
    int m = numRows,
        // n = C.numCols,
        n = numCols,
        k = C.numRows;
    DGEMM("T", "N", &m, &n, &k,&otherFact, B.data, &  k,
                                           C.data, &  k,
                               &thisFact,    data, &  m);
    return 0;
  }
#endif


  if (thisFact == 1.0) {
    int numMults = C.numRows;
    double *aijPtr = data;
    for (int j=0; j<numCols; j++) {
      for (int i=0; i<numRows; i++) {
        double *bkiPtr  = &(B.data)[i*numMults];
        double *cjkPtr  = &(C.data)[j*numMults];
        double sum = 0.0;
        for (int k=0; k<numMults; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr++ += sum * otherFact;
      }
    } 
  } else if (thisFact == 0.0) {
    int numMults = C.numRows;
    double *aijPtr = data;
    for (int j=0; j<numCols; j++) {
      for (int i=0; i<numRows; i++) {
        double *bkiPtr  = &(B.data)[i*numMults];
        double *cjkPtr  = &(C.data)[j*numMults];
        double sum = 0.0;
        for (int k=0; k<numMults; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr++ = sum * otherFact;
      }
    } 
  } else {
    int numMults = C.numRows;
    double *aijPtr = data;
    for (int j=0; j<numCols; j++) {
      for (int i=0; i<numRows; i++) {
        double *bkiPtr  = &(B.data)[i*numMults];
        double *cjkPtr  = &(C.data)[j*numMults];
        double sum = 0.0;
        for (int k=0; k<numMults; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr = *aijPtr * thisFact + sum * otherFact;
        aijPtr++;
      }
    } 
  }

  return 0;
}


// to perform this += T' * B * T
int
Matrix::addMatrixTripleProduct(double thisFact, 
                               const  Matrix &T, 
                               const  Matrix &B, 
                               double otherFact)
{
  if (thisFact == 1.0 && otherFact == 0.0)
    return 0;

  // check work area can hold the temporary matrix
  int dimB = B.numCols;
  int sizeWork = dimB * numCols;

  if (sizeWork > sizeDoubleWork) {
    // TODO
    this->addMatrix(thisFact, T^B*T, otherFact);
    return 0;
  }
  else {
    int m = B.numRows,
        n = T.numCols,
        k = B.numCols;
      //k = T.numRows;
    double zero = 0.0,
           one  = 1.0;

    DGEMM ("N", "N", &m      , &n      , &k,&one      , B.data, &B.numRows, // m
                                                        T.data, &T.numRows, // k
                                            &zero,  matrixWork, &m);

    DGEMM ("T", "N", &numRows, &numCols, &k,&otherFact, T.data, &T.numRows,
                                                    matrixWork, &m, // k
                                            &thisFact,    data, &numRows);
    return 0;
  }

  // zero out the work area
  double *matrixWorkPtr = matrixWork;
  for (int l=0; l<sizeWork; l++)
    *matrixWorkPtr++ = 0.0;
  
  // now form B * T * fact store in matrixWork == A area
  // NOTE: looping as per blas3 DGEMM : j,k,i

  double *tkjPtr  = &(T.data)[0];
  for (int j=0; j<numCols; j++) {
    double *aijPtrA = &matrixWork[j*dimB];
    for (int k=0; k<dimB; k++) {
      double tmp = *tkjPtr++ * otherFact;
      double *aijPtr = aijPtrA;
      double *bikPtr = &(B.data)[k*dimB];
      for (int i=0; i<dimB; i++) 
        *aijPtr++ += *bikPtr++ * tmp;
    }
  }

  // now form T' * matrixWork
  // NOTE: looping as per blas3 DGEMM : j,i,k
  if (thisFact == 1.0) {
    double *dataPtr = &data[0];
    for (int j=0; j< numCols; j++) {
      double *workkjPtrA = &matrixWork[j*dimB];
      for (int i=0; i<numRows; i++) {
        double *ckiPtr = &(T.data)[i*dimB];
        double *workkjPtr = workkjPtrA;
        double aij = 0.0;
        for (int k=0; k< dimB; k++)
          aij += *ckiPtr++ * *workkjPtr++;
        *dataPtr++ += aij;
      }
    }
  } else if (thisFact == 0.0) {
    double *dataPtr = &data[0];
    for (int j=0; j< numCols; j++) {
      double *workkjPtrA = &matrixWork[j*dimB];
      for (int i=0; i<numRows; i++) {
        double *ckiPtr = &(T.data)[i*dimB];
        double *workkjPtr = workkjPtrA;
        double aij = 0.0;
        for (int k=0; k< dimB; k++)
          aij += *ckiPtr++ * *workkjPtr++;
        *dataPtr++ = aij;
      }
    }

  } else {
    double *dataPtr = &data[0];
    for (int j=0; j< numCols; j++) {
      double *workkjPtrA = &matrixWork[j*dimB];
      for (int i=0; i<numRows; i++) {
        double *ckiPtr = &(T.data)[i*dimB];
        double *workkjPtr = workkjPtrA;
        double aij = 0.0;
        for (int k=0; k< dimB; k++)
          aij += *ckiPtr++ * *workkjPtr++;
        double value = *dataPtr * thisFact + aij;
        *dataPtr++ = value;
      }
    }
  }

  return 0;
}


//
// to perform this += At * B * C
//
int
Matrix::addMatrixTripleProduct(double thisFact, 
                               const Matrix &A, 
                               const Matrix &B,
                               const Matrix &C,
                               double otherFact)
{
    if (thisFact == 1.0 && otherFact == 0.0)
      return 0;

    // check work area can hold the temporary matrix
    int sizeWork = B.numRows * numCols;
#ifndef NO_WORK
    if (sizeWork > sizeDoubleWork) {
#endif
      this->addMatrix(thisFact, A^B*C, otherFact);
      return 0;
#ifndef NO_WORK
    }

    // zero out the work area
    double *matrixWorkPtr = matrixWork;
    for (int l=0; l<sizeWork; l++)
      *matrixWorkPtr++ = 0.0;

    // now form B * C * fact store in matrixWork == A area
    // NOTE: looping as per blas3 DGEMM : j,k,i

    int rowsB = B.numRows;
    double *ckjPtr  = &(C.data)[0];
    for (int j=0; j<numCols; j++) {
      double *aijPtrA = &matrixWork[j*rowsB];
      for (int k=0; k<rowsB; k++) {
        double tmp = *ckjPtr++ * otherFact;
        double *aijPtr = aijPtrA;
        double *bikPtr = &(B.data)[k*rowsB];
        for (int i=0; i<rowsB; i++) 
          *aijPtr++ += *bikPtr++ * tmp;
      }
    }

    // now form A' * matrixWork
    // NOTE: looping as per blas3 DGEMM : j,i,k
    int dimB = rowsB;
    if (thisFact == 1.0) {
      double *dataPtr = &data[0];
      for (int j=0; j< numCols; j++) {
        double *workkjPtrA = &matrixWork[j*dimB];
        for (int i=0; i<numRows; i++) {
          double *akiPtr = &(A.data)[i*dimB];
          double *workkjPtr = workkjPtrA;
          double aij = 0.0;
          for (int k=0; k< dimB; k++)
            aij += *akiPtr++ * *workkjPtr++;
          *dataPtr++ += aij;
        }
      }
    } else if (thisFact == 0.0) {
      double *dataPtr = &data[0];
      for (int j=0; j< numCols; j++) {
        double *workkjPtrA = &matrixWork[j*dimB];
        for (int i=0; i<numRows; i++) {
          double *akiPtr = &(A.data)[i*dimB];
          double *workkjPtr = workkjPtrA;
          double aij = 0.0;
          for (int k=0; k< dimB; k++)
            aij += *akiPtr++ * *workkjPtr++;
          *dataPtr++ = aij;
        }
      }

    } else {
      double *dataPtr = &data[0];
      for (int j=0; j< numCols; j++) {
        double *workkjPtrA = &matrixWork[j*dimB];
        for (int i=0; i<numRows; i++) {
          double *akiPtr = &(A.data)[i*dimB];
          double *workkjPtr = workkjPtrA;
          double aij = 0.0;
          for (int k=0; k< dimB; k++)
            aij += *akiPtr++ * *workkjPtr++;
          double value = *dataPtr * thisFact + aij;
          *dataPtr++ = value;
        }
      }
    }

    return 0;
#endif
}


int
Matrix::Assemble(const Matrix &V, int init_row, int init_col, double fact) 
{
  int pos_Rows, pos_Cols;
  int res = 0;
  
  int VnumRows = V.numRows;
  int VnumCols = V.numCols;
  
  [[maybe_unused]] int final_row = init_row + VnumRows - 1;
  [[maybe_unused]] int final_col = init_col + VnumCols - 1;
  
  assert((init_row >= 0) && (final_row < numRows) && (init_col >= 0) && (final_col < numCols));

  for (int i=0; i<VnumCols; i++) {
     pos_Cols = init_col + i;
     for (int j=0; j<VnumRows; j++) {
        pos_Rows = init_row + j;
   
        (*this)(pos_Rows,pos_Cols) += V(j,i)*fact;
     }
  }  

  return res;
}


int
Matrix::Assemble(const Vector &V, int init_row, int init_col, double fact) 
{

  const int VnumRows = V.sz;
  const int VnumCols = 1;
 
  {
    [[maybe_unused]] int final_row = init_row + VnumRows - 1;
    [[maybe_unused]] int final_col = init_col + VnumCols - 1;

    assert((init_row >= 0) && (final_row < numRows) && (init_col >= 0) && (final_col < numCols));
  }

  int pos_Rows, pos_Cols;
  int res = 0;
  
  for (int i=0; i<VnumCols; i++) {
     pos_Cols = init_col + i;
     for (int j=0; j<VnumRows; j++) {
        pos_Rows = init_row + j;
   
        (*this)(pos_Rows,pos_Cols) += V(j)*fact;
     }
  }

  return res;
}


int
Matrix::AssembleTranspose(const Matrix &V, int init_row, int init_col, double fact) 
{
  int VnumRows = V.numRows;
  int VnumCols = V.numCols;
  
  {
    [[maybe_unused]] int final_row = init_row + VnumCols - 1;
    [[maybe_unused]] int final_col = init_col + VnumRows - 1; 
    assert((init_row >= 0) && (final_row < numRows) && (init_col >= 0) && (final_col < numCols));
  }

  int pos_Rows, pos_Cols;
  int res = 0;

  for (int i=0; i<VnumRows; i++) {
     pos_Cols = init_col + i;
     for (int j=0; j<VnumCols; j++) {
        pos_Rows = init_row + j; 
        (*this)(pos_Rows,pos_Cols) += V(i,j)*fact;
     }
  }

  return res;
}


int
Matrix::AssembleTranspose(const Vector &V, int init_row, int init_col, double fact) 
{
  int VnumRows = V.sz;
  int VnumCols = 1;
  { 
    [[maybe_unused]] int final_row = init_row + VnumCols - 1;
    [[maybe_unused]] int final_col = init_col + VnumRows - 1;
    
    assert((init_row >= 0) && (final_row < numRows) && (init_col >= 0) && (final_col < numCols));
  }

  int pos_Rows, pos_Cols;
  int res = 0;
  
  for (int i=0; i<VnumRows; i++) {
     pos_Cols = init_col + i;
     for (int j=0; j<VnumCols; j++) {
        pos_Rows = init_row + j;   
        (*this)(pos_Rows,pos_Cols) += V(i)*fact;
     }
  }

  return res;
}


int
Matrix::Extract(const Matrix &V, int init_row, int init_col, double fact) 
{
  [[maybe_unused]] int final_row = init_row + numRows - 1;
  [[maybe_unused]] int final_col = init_col + numCols - 1;  
  assert((init_row >= 0) && (final_row < V.numRows) && 
         (init_col >= 0) && (final_col < V.numCols));

  int res = 0;
  int pos_Rows, pos_Cols;
  for (int i=0; i<numCols; i++) {
     pos_Cols = init_col + i;
     for (int j=0; j<numRows; j++) {
        pos_Rows = init_row + j;
   
        (*this)(j,i) = V(pos_Rows,pos_Cols)*fact;
     }
  }

  return res;
}


//
// OVERLOADED OPERATOR () to CONSTRUCT A NEW MATRIX
//

Matrix
Matrix::operator()(const ID &rows, const ID & cols) const
{
    int nRows, nCols;
    nRows = rows.Size();
    nCols = cols.Size();
    Matrix result(nRows,nCols);
    double *dataPtr = result.data;
    for (int i=0; i<nCols; i++)
        for (int j=0; j<nRows; j++)
            *dataPtr++ = (*this)(rows(j),cols(i));

    return result;
}
                
// Matrix &operator=(const Matrix  &V):
//      the assignment operator, This is assigned to be a copy of V. if sizes
//      are not compatable this.data [] is deleted. The data pointers will not
//      point to the same area in mem after the assignment.
//
Matrix &
Matrix::operator=(const Matrix &other)
{
  // first check we are not trying other = other
  if (this == &other) 
    return *this;

  if ((numCols != other.numCols) || (numRows != other.numRows)) {
#ifdef _G3DEBUG
      opserr << "Matrix::operator=() - matrix dimensions do not match\n";
#endif

      if (this->data != 0) {
          delete [] this->data;
          this->data = 0;
      }

      int theSize = other.numCols*other.numRows;

      data = new double[theSize];

      this->dataSize = theSize;
      this->numCols  = other.numCols;
      this->numRows  = other.numRows;
  }

  // now copy the data
  double *dataPtr = data;
  double *otherDataPtr = other.data;                    
  for (int i=0; i<dataSize; i++)
      *dataPtr++ = *otherDataPtr++;
  
  return *this;
}


// Move assignment
//
#if !defined(NO_CXX11_MOVE)
Matrix &
Matrix::operator=( Matrix &&other)
{
  // first check we are not trying other = other
  if (this == &other) 
    return *this;

  if (this->data != 0 && fromFree == 0){
    delete [] this->data;
    this->data = 0;
  }
        
  this->data = other.data;
  this->dataSize = other.numCols*other.numRows;
  this->numCols  = other.numCols;
  this->numRows  = other.numRows;
  this->fromFree = other.fromFree;
  other.data     = 0;
  other.dataSize = 0;
  other.numCols  = 0;
  other.numRows  = 0;
  other.fromFree = 1;

  return *this;
}
#endif


// virtual Matrix &operator+=(double fact);
// virtual Matrix &operator-=(double fact);
// virtual Matrix &operator*=(double fact);
// virtual Matrix &operator/=(double fact); 
//        The above methods all modify the current matrix. If in
//        derived matrices data kept in data and of sizeData no redef necessary.

Matrix &
Matrix::operator+=(double fact)
{
  // check if quick return
  if (fact == 0.0)
    return *this;

  double *dataPtr = data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ += fact;
  
  return *this;
}


Matrix &
Matrix::operator-=(double fact)
{
  // check if quick return
  if (fact == 0.0)
    return *this;
  
  double *dataPtr = data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ -= fact;

  return *this;
}


Matrix &
Matrix::operator*=(double fact)
{
  // check if quick return
  if (fact == 1.0)
    return *this;
  
  double *dataPtr = data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ *= fact;
  
  return *this;
}

Matrix &
Matrix::operator/=(double fact)
{
    // check if quick return
    if (fact == 1.0)
        return *this;

    else {
      double val = 1.0/fact;

      double *dataPtr = data;
      for (int i=0; i<dataSize; i++)
        *dataPtr++ *= val;

      return *this;
    }
}


//
//    virtual Matrix operator+(double fact);
//    virtual Matrix operator-(double fact);
//    virtual Matrix operator*(double fact);
//    virtual Matrix operator/(double fact);
//        The above methods all return a new full general matrix.
//

Matrix
Matrix::operator+(double fact) const
{
    Matrix result(*this);
    result += fact;
    return result;
}

Matrix
Matrix::operator-(double fact) const
{
    Matrix result(*this);
    result -= fact;
    return result;
}

Matrix
Matrix::operator*(double fact) const
{
    Matrix result(*this);
    result *= fact;
    return result;
}

Matrix
Matrix::operator/(double fact) const
{
    Matrix result(*this);
    result /= fact;
    return result;
}


//
// MATRIX_VECTOR OPERATIONS
//
Vector
Matrix::operator*(const Vector3D &V) const
{
    Vector result(numRows);

    double *dataPtr = data;
    for (int i=0; i<numCols; i++)
      for (int j=0; j<numRows; j++)
        result(j) += *dataPtr++ * V[i];

    return result;
}

Vector
Matrix::operator*(const Vector &V) const
{
    Vector result(numRows);
#ifdef MATRIX_BLAS
    result.addMatrixVector(0.0, *this, V, 1.0);
    return result;
#else

    double *dataPtr = data;
    for (int i=0; i<numCols; i++)
      for (int j=0; j<numRows; j++)
        result(j) += *dataPtr++ * V(i);

    return result;
#endif
}

Vector
Matrix::operator^(const Vector &V) const
{
    assert(V.Size() == numRows);

    Vector result(numCols);
#ifdef MATRIX_BLAS
    result.addMatrixTransposeVector(0.0, *this, V, 1.0);
    return result;
#else

    double *dataPtr = data;
    for (int i=0; i<numCols; i++)
      for (int j=0; j<numRows; j++)
        result(i) += *dataPtr++ * V(j);

    return result;
#endif
}


//
// MATRIX - MATRIX OPERATIONS
//
Matrix
Matrix::operator+(const Matrix &M) const
{
    Matrix result(*this);
    result.addMatrix(M,1.0);    
    return result;
}
            
Matrix
Matrix::operator-(const Matrix &M) const
{
    Matrix result(*this);
    result.addMatrix(M,-1.0);    
    return result;
}
            
    
Matrix
Matrix::operator*(const Matrix &M) const
{
    Matrix result(numRows,M.numCols);
    result.addMatrixProduct(0.0, *this, M, 1.0);
    return result;
}



// Matrix operator^(const Matrix &M) const
//        We overload the * operator to perform matrix^t-matrix multiplication.
//        reults = (*this)transposed * M.

Matrix
Matrix::operator^(const Matrix &M) const
{
  Matrix result(numCols,M.numCols);
#ifdef MATRIX_BLAS
  result.addMatrixTransposeProduct(0.0, *this, M, 1.0);
#else

  assert(numRows == M.numRows && result.numRows == numCols);

    double *resDataPtr = result.data;            

    int innerDim = numRows;
    int nCols = result.numCols;
    for (int i=0; i<nCols; i++) {
      double *aDataPtr = data;
      double *bStartColDataPtr = &(M.data[i*innerDim]);
      for (int j=0; j<numCols; j++) {
        double *bDataPtr = bStartColDataPtr;
        double sum = 0.0;
        for (int k=0; k<innerDim; k++) {
          sum += *aDataPtr++ * *bDataPtr++;
        }
        *resDataPtr++ = sum;
      }
    }
#endif
    return result;
}
    

Matrix &
Matrix::operator+=(const Matrix &M)
{
  assert(numRows == M.numRows && numCols == M.numCols);

#ifdef MATRIX_BLAS
  addMatrix( M, 1.0);
#else
  double *dataPtr = data;
  double *otherData = M.data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ += *otherData++;
#endif 
  return *this;
}

Matrix &
Matrix::operator-=(const Matrix &M)
{
#ifdef MATRIX_BLAS
  addMatrix(M, -1.0);
#else

  assert(numRows == M.numRows && numCols == M.numCols);

  double *dataPtr = data;
  double *otherData = M.data;
  for (int i=0; i<dataSize; i++)
    *dataPtr++ -= *otherData++;
#endif 
  return *this;
}


//
// Input/Output Methods
//

void 
Matrix::Output(OPS_Stream &s) const
{
    for (int i=0; i<noRows(); i++) {
        for (int j=0; j<noCols(); j++)
            s <<  (*this)(i,j) << " ";
        s << "\n";
    }
}


//
// friend stream functions for input and output
//

OPS_Stream &operator<<(OPS_Stream &s, const Matrix &V)
{
    s << "\n";
    V.Output(s);
    s << "\n";        
    return s;
}


Matrix operator*(double a, const Matrix &V)
{
  return V * a;
}

Vector Matrix::diagonal() const
{
  
  assert(numRows == numCols);

  int size = numRows < numCols ? numRows : numCols;
  Vector diagonal(size);

  for (int i = 0; i < size; ++i)
    diagonal(i) = data[i*numRows + i];

  return diagonal;
}

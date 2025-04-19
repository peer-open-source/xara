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
// Description: This file contains the class implementation for Vector.
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#include "Vector.h"
#include "Matrix.h"
#include "ID.h"
#include <OPS_Stream.h>

#include <math.h>
#include <assert.h>
#include "blasdecl.h"

#if 0
#define VECTOR_BLAS
#endif

// Vector():
//        Standard constructor, sets size = 0;

Vector::Vector()
: sz(0), theData(nullptr), fromFree(0)
{

}

// Vector(int size):
//        Constructor used to allocate a vector of size size.

Vector::Vector(int size)
: sz(size), theData(nullptr), fromFree(0)
{
  assert(size >= 0);

  // get some space for the vector
  if (size > 0)
    theData = new double [size]{};
}

Vector::Vector(std::shared_ptr<double[]> data, int size)
: sz(size), theData(nullptr), fromFree(0)
{
  if (size > 0) {
    theData = new double [size];

    for (int i=0; i<sz; i++)
      theData[i] = data[i];
  }
}


Vector::Vector(double *data, int size)
: sz(size),theData(data),fromFree(1)
{
  assert(sz >= 0);
}
 

// Vector(const Vector&):
//        Constructor to init a vector from another.

Vector::Vector(const Vector &other)
: sz(other.sz),theData(0),fromFree(0)
{
  if (sz != 0) {
    theData = new double [other.sz];    
  }
  // copy the component data
  for (int i=0; i<sz; i++)
    theData[i] = other.theData[i];
}

// Vector(const Vector&):
//  Move constructor
#if !defined(NO_CXX11_MOVE)   
Vector::Vector(Vector &&other)
: sz(other.sz),theData(other.theData),fromFree(0)
{
  other.theData = nullptr;
  other.sz = 0;
} 
#endif



Vector::~Vector()
{
  if (fromFree == 0 && theData != nullptr)
    delete [] theData;
  theData = nullptr;
}


int 
Vector::setData(double *newData, int size)
{
  assert(size >  0);

  if (theData != nullptr && fromFree == 0) {
    delete [] theData;
    theData = nullptr;
  }
  sz = size;
  theData = newData;
  fromFree = 1;

  return 0;
}


Vector
Vector::view(int start, int end) const
{
  assert(end < sz);
  assert(start < sz);

  return Vector (&(this->theData[start]), end-start);
}
 

int 
Vector::resize(int newSize)
{
  // first check that newSize is valid
  assert(newSize >= 0);
  
  // otherwise if newSize is gretaer than oldSize free old space and get new space
  if (newSize > sz) {

    // delete the old array
    if (theData != 0 && fromFree == 0) {
      delete [] theData;
      theData = nullptr;
    }
    sz = 0;
    fromFree = 0;
    
    // create new memory
    // theData = (double *)malloc(newSize*sizeof(double));    
    theData = new double[newSize];

    sz = newSize;
  }  

  // just set the size to be newSize .. penalty of holding onto additional
  // memory .. but then save the free() and malloc() calls
  else 
    sz = newSize;    

  return 0;
}


// Assemble(Vector &x, ID &y, double fact ):
//     Method to assemble into object the Vector V using the ID l.
//     If ID(x) does not exist program writes error message if
//     VECTOR_CHECK defined, otherwise ignores it and goes on.

int 
Vector::Assemble(const Vector &V, const ID &l, double fact )
{
  int result = 0;
  int pos;
  for (int i=0; i<l.Size(); i++) {
    pos = l(i);
    
    if (pos < 0)
      ;
    else {
      assert((pos < sz) && (i < V.Size()));
      // assemble into vector
      theData[pos] += V.theData[i] *fact;
    }
  }

  return result;
}
    

int
Vector::Normalize() 
{
  double length = Norm();
  /*
  for (int i=0; i<sz; i++)
    length += theData[i] * theData[i];
  length = sqrt(length);
  */
  
  if (length == 0.0) 
    return -1;

  length = 1.0/length;
#ifdef VECTOR_BLAS
  const int incx = 1;
  dscal_(&sz, &length, theData, &incx);
#else
  for (int j=0; j<sz; j++)
    theData[j] *= length;
#endif
  return 0;
}

int
Vector::addVector(const Vector &other, double otherFact )
{
  // if sizes are compatable add
  assert(sz == other.sz);

  // check if quick return
  if (otherFact == 0.0)
    return 0; 

  else {
#ifdef VECTOR_BLAS
    const int incr = 1;
    daxpy_(&sz, &otherFact, other.theData, &incr, theData, &incr);
    return 0;
#else
    // want: this += other * otherFact
    double *dataPtr = theData;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
        *dataPtr++ += *otherDataPtr++;
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
        *dataPtr++ -= *otherDataPtr++;
    } else 
      for (int i=0; i<sz; i++) 
        *dataPtr++ += *otherDataPtr++ * otherFact;
#endif
  }

  // successfull
  return 0;
}
int
Vector::addVector(double thisFact, const Vector &other, double otherFact )
{
  // if sizes are compatable add
  assert(sz == other.sz);

  // check if quick return
  if (otherFact == 0.0 && thisFact == 1.0)
    return 0; 

  else if (thisFact == 1.0) {
#ifdef VECTOR_BLAS
    const int incr = 1;
    daxpy_(&sz, &otherFact, other.theData, &incr, theData, &incr);
    return 0;
#else
    // want: this += other * otherFact
    double *dataPtr = theData;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
        *dataPtr++ += *otherDataPtr++;
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
        *dataPtr++ -= *otherDataPtr++;
    } else 
      for (int i=0; i<sz; i++) 
        *dataPtr++ += *otherDataPtr++ * otherFact;
#endif
  }

  else if (thisFact == 0.0) {

    // want: this = other * otherFact
    double *dataPtr = theData;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
        *dataPtr++ = *otherDataPtr++;
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
        *dataPtr++ = -(*otherDataPtr++);
    } else 
      for (int i=0; i<sz; i++) 
        *dataPtr++ = *otherDataPtr++ * otherFact;
  }

  else {

    // want: this = this * thisFact + other * otherFact
    double *dataPtr = theData;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) {
        double value = *dataPtr * thisFact + *otherDataPtr++;
        *dataPtr++ = value;
      }
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) {
        double value = *dataPtr * thisFact - *otherDataPtr++;
        *dataPtr++ = value;
      }
    } else 
      for (int i=0; i<sz; i++) {
        double value = *dataPtr * thisFact + *otherDataPtr++ * otherFact;
        *dataPtr++ = value;
      }
  } 

  // successfull
  return 0;
}


int
Vector::addMatrixVector(double thisFact, const Matrix &m, const Vector &v, double otherFact )
{
  // check the sizes are compatable
  assert(sz == m.noRows());
  assert(m.noCols() == v.sz);

  // see if quick return
  if (thisFact == 1.0 && otherFact == 0.0)
    return 0;

#ifdef VECTOR_BLAS
  else if (v.sz > 10) {
    int incr = 1,
           i = m.numRows,
           n = m.numCols;
    return
      DGEMV("N", &i, &n,
            &otherFact,
            m.data, &i,
            v.theData, &incr,
            &thisFact,
            theData,   &incr);
  }
#endif

  else if (thisFact == 1.0) {

    // want: this += m * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++;
        for (int j=0; j<sz; j++)
          theData[j] += *matrixDataPtr++ * otherData;
      }
    } 
    else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++;
        for (int j=0; j<sz; j++)
          theData[j] -= *matrixDataPtr++ * otherData;
      }
    } 
    else { // have to do the multiplication
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++ * otherFact;
        for (int j=0; j<sz; j++)
          theData[j] += *matrixDataPtr++ * otherData;
      }
    }
  }

  else if (thisFact == 0.0) {
    
    // want: this = m * v * otherFact
    for (int i=0; i<sz; i++)
      theData[i] = 0.0;

    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++;
        for (int j=0; j<sz; j++)
          theData[j] += *matrixDataPtr++ * otherData;
      }
    } 
    else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++;
        for (int j=0; j<sz; j++)
          theData[j] -= *matrixDataPtr++ * otherData;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++ * otherFact;
        for (int j=0; j<sz; j++)
          theData[j] += *matrixDataPtr++ * otherData;
      }
    }
  }

  else {

    // want: this = this * thisFact + m * v * otherFact
    for (int i=0; i<sz; i++)
      theData[i] *= thisFact;

    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++;
        for (int j=0; j<sz; j++)
          theData[j] += *matrixDataPtr++ * otherData;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++;
        for (int j=0; j<sz; j++)
          theData[j] -= *matrixDataPtr++ * otherData;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
        double otherData = *otherDataPtr++ * otherFact;
        for (int j=0; j<sz; j++)
          theData[j] += *matrixDataPtr++ * otherData;
      }
    }
  }
  
  // successfull
  return 0;
}



int
Vector::addMatrixTransposeVector(double thisFact, 
                                 const Matrix &m, 
                                 const Vector &v, 
                                 double otherFact )
{
#ifdef _G3DEBUG
  // check the sizes are compatable
  if ((sz != m.noRows()) && (m.noRows() != v.sz)) {
    opserr << "Vector::addMatrixTransposeVector() - incompatable sizes\n";
    return -1;    
  }
#endif

  // see if quick return
  if (otherFact == 0.0 && thisFact == 1.0)
    return 0;

#ifdef VECTOR_BLAS
  else if (v.sz > 10) {
    int incr = 1,
           i = m.numRows,
           n = m.numCols;
    return
      DGEMV("T", &i, &n,
            &otherFact,
            m.data, &i,
            v.theData, &incr,
            &thisFact,
            theData,   &incr);
  }
#endif

  else if (thisFact == 1.0) {

    // want: this += m^t * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
        double *otherDataPtr = otherDataPtrA;
        double sum = 0.0;
        for (int j=0; j<otherSize; j++)
          sum += *matrixDataPtr++ * *otherDataPtr++;
        theData[i] += sum;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
        double *otherDataPtr = otherDataPtrA;
        double sum = 0.0;
        for (int j=0; j<otherSize; j++)
          sum += *matrixDataPtr++ * *otherDataPtr++;
        theData[i] -= sum;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
        double *otherDataPtr = otherDataPtrA;
        double sum = 0.0;
        for (int j=0; j<otherSize; j++)
          sum += *matrixDataPtr++ * *otherDataPtr++;
        theData[i] += sum * otherFact;
      }
    }
  }

  else if (thisFact == 0.0) {

    // want: this = m^t * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
        double *otherDataPtr = otherDataPtrA;
        double sum = 0.0;
        for (int j=0; j<otherSize; j++)
          sum += *matrixDataPtr++ * *otherDataPtr++;
        theData[i] = sum;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
        double *otherDataPtr = otherDataPtrA;
        double sum = 0.0;
        for (int j=0; j<otherSize; j++)
          sum += *matrixDataPtr++ * *otherDataPtr++;
        theData[i] = -sum;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
        double *otherDataPtr = otherDataPtrA;
        double sum = 0.0;
        for (int j=0; j<otherSize; j++)
          sum += *matrixDataPtr++ * *otherDataPtr++;
        theData[i] = sum * otherFact;
      }
    }
  } 

  else {

    // want: this = this * thisFact + m^t * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
        double *otherDataPtr = otherDataPtrA;
        double sum = 0.0;
        for (int j=0; j<otherSize; j++)
          sum += *matrixDataPtr++ * *otherDataPtr++;
        double value = theData[i] * thisFact + sum;
        theData[i] = value;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
        double *otherDataPtr = otherDataPtrA;
        double sum = 0.0;
        for (int j=0; j<otherSize; j++)
          sum += *matrixDataPtr++ * *otherDataPtr++;
        double value = theData[i] * thisFact - sum;
        theData[i] = value;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
        double *otherDataPtr = otherDataPtrA;
        double sum = 0.0;
        for (int j=0; j<otherSize; j++)
          sum += *matrixDataPtr++ * *otherDataPtr++;
        double value = theData[i] * thisFact + sum * otherFact;
        theData[i] = value;
      }
    }
}

  return 0;
}
        
        
// double Norm();
//  return the norm of  vector. (non-const as may save norm for later)

double
Vector::Norm() const
{
  double value = 0;
  for (int i=0; i<sz; i++) {
    double data = theData[i];
    value += data*data;
  }
  return sqrt(value);
}

double
Vector::pNorm(int p) const
{
  double value = 0;
  
  if (p > 0) {
    for (int i=0; i<sz; i++) {
      double data = fabs(theData[i]);
      value += pow(data,p);
    }
    return pow(value,1.0/p);
  } else {
    for (int i=0; i<sz; i++) {
      double data = fabs(theData[i]);
      value = (data>value) ? data : value;
    }
    return value;
  }
}


double &
Vector::operator[](int x) 
{
  // check if it is inside range [0,sz-1]
  assert(x >= 0);
  
  if (x >= sz) {
    // TODO: Is this expected?
    double *dataNew = new double[x+1];
    for (int i=0; i<sz; i++)
      dataNew[i] = theData[i];
    for (int j=sz; j<x; j++)
      dataNew[j] = 0.0;
    
    if (fromFree == 0)
      if (theData != 0){
        delete [] theData;
      theData = 0;
    }
    theData = dataNew;
    sz = x+1;
  }

  return theData[x];
}



// operator()(const ID &rows) const
//  return a vector whose components are the components of the
//        current vector located in positions given by the ID rows.
Vector 
Vector::operator()(const ID &rows) const
{
  // create a new Vector to be returned
  Vector result(rows.Size());

  // copy the appropraite contents from current to result     
  int pos;
  for (int i=0; i<rows.Size(); i++) {
    pos = rows(i);
    assert(pos >= 0 && pos < sz);
    result(i) = (*this)(pos);
  }
  return result;
}


// Vector &operator=(const Vector  &V):
//        the assignment operator, This is assigned to be a copy of V. if sizes
//        are not compatable this.theData [] is deleted. The data pointers will not
//        point to the same area in mem after the assignment.
//

Vector &
Vector::operator=(const Vector &V) 
{
  // assert(sz == V.sz); // TODO
  // first check we are not trying v = v
  if (this != &V) {

      if (sz != V.sz)  {
          // Check that we are not deleting an empty Vector
          if (this->theData != nullptr){
            delete [] this->theData;
            this->theData = nullptr;
          }
          this->sz = V.sz;
          
          // Check that we are not creating an empty Vector
          this->theData = (sz != 0) ? new double[sz] : nullptr;
      }

      // copy the data
      for (int i=0; i<sz; i++)
        theData[i] = V.theData[i];
  }

  return *this;
}

// Move assignment operator.  
#if !defined(NO_CXX11_MOVE)   
Vector &
Vector::operator=(Vector &&V) 
{
  // first check we are not trying v = v
  if (this != &V) {
    if (this->theData != nullptr) { 
      delete [] this->theData;
      this->theData = 0;
    }
    theData = V.theData;
    this->sz = V.sz;
    V.theData = 0;
    V.sz = 0;
  }
  return *this;
}
#endif


// Vector &operator+=(double fact):
//        The += operator adds fact to each element of the vector, data[i] = data[i]+fact.

Vector &Vector::operator+=(double fact)
{
  if (fact != 0.0)
    for (int i=0; i<sz; i++)
      theData[i] += fact;
  return *this;    
}


// Vector &operator-=(double fact)
//        The -= operator subtracts fact from each element of the vector, data[i] = data[i]-fact.

Vector &Vector::operator-=(double fact)
{
  if (fact != 0.0)
    for (int i=0; i<sz; i++)
      theData[i] -= fact;
  return *this;    
}



// Vector &operator*=(double fact):
//        The -= operator subtracts fact from each element of the vector, data[i] = data[i]-fact.

Vector &Vector::operator*=(double fact)
{
#ifdef VECTOR_BLAS
  const int incr = 1;
  dscal_(&sz, &fact, theData, &incr);
  return *this;
#else
  for (int i=0; i<sz; i++)
    theData[i] *= fact;
  return *this;
#endif
}


// Vector &operator/=(double fact):
//        The /= operator divides each element of the vector by fact, theData[i] = theData[i]/fact.
//        Program exits if divide-by-zero would occur with warning message.

Vector &Vector::operator/=(double fact)
{
#if 0
  if (fact == 0.0) { // instead of divide-by-zero error set to VECTOR_VERY_LARGE_VALUE
    for (int i=0; i<sz; i++)
      theData[i] = VECTOR_VERY_LARGE_VALUE;
  } else 
#endif
  {
    for (int i=0; i<sz; i++)
      theData[i] /= fact;
  }
  
  return *this;
}




// Vector operator+(double fact):
//        The + operator returns a Vector of the same size as current, whose components
//        are return(i) = theData[i]+fact;

Vector 
Vector::operator+(double fact) const
{
  Vector result(*this);
  result += fact;
  return result;
}



// Vector operator-(double fact):
//        The + operator returns a Vector of the same size as current, whose components
//        are return(i) = theData[i]-fact;

Vector 
Vector::operator-(double fact) const
{
    Vector result(*this);
    result -= fact;
    return result;
}



// Vector operator*(double fact):
//        The + operator returns a Vector of the same size as current, whose components
//        are return(i) = theData[i]*fact;

Vector 
Vector::operator*(double fact) const
{
    Vector result(*this);
    result *= fact;
    return result;
}


// Vector operator/(double fact):
//        The + operator returns a Vector of the same size as current, whose components
//        are return(i) = theData[i]/fact; Exits if divide-by-zero error.

Vector 
Vector::operator/(double fact) const
{
    assert(fact != 0.0);
    Vector result(*this);
    result /= fact;
    return result;
}



// Vector &operator+=(const Vector &V):
//        The += operator adds V's data to data, data[i]+=V(i). A check to see if
//        vectors are of same size is performed if VECTOR_CHECK is defined.

Vector &
Vector::operator+=(const Vector &other)
{
  assert(sz == other.sz);
  for (int i=0; i<sz; i++)
    theData[i] += other.theData[i];
  return *this;
}


// Vector &operator-=(const Vector &V):
//        The -= operator subtracts V's data from  data, data[i]+=V(i). A check 
//           to see if vectors are of same size is performed if VECTOR_CHECK is defined.
Vector &
Vector::operator-=(const Vector &other)
{
  assert(sz == other.sz);
  for (int i=0; i<sz; i++)
    theData[i] -= other.theData[i];
  return *this;    
}


// Vector operator+(const Vector &V):
//        The + operator checks the two vectors are of the same size if VECTOR_CHECK is defined.
//         Then returns a Vector whose components are the vector sum of current and V's data.

Vector 
Vector::operator+(const Vector &b) const
{
  assert(sz == b.sz);
  Vector result(*this);
  result += b;
  return result;
}


// Vector operator-(const Vector &V):
//        The - operator checks the two vectors are of the same size and then returns a Vector
//        whose components are the vector difference of current and V's data.

Vector 
Vector::operator-(const Vector &b) const
{
  assert(sz == b.sz);
  Vector result(*this);
  result -= b;
  return result;
}



// double operator^(const Vector &V) const;
//  perform (Vector)transposed * vector.
double
Vector::operator^(const Vector &V) const
{
  assert(sz == V.sz);

  double result = 0.0;
  double *dataThis = theData;
  double *dataV = V.theData;
    for (int i=0; i<sz; i++)
      result += *dataThis++ * *dataV++;

    return result;
}


// Vector operator/(const Matrix &M) const;    
//  return inv(M)*this
#if 1
Vector
Vector::operator/(const Matrix &M) const
{
  Vector res(M.noRows());
    
  if (M.noRows() != M.noCols()) { // if not square do least squares solution
    Matrix A(M^M);
    A.Solve(*this, res);    
  }
  else {
    M.Solve(*this, res);
  }
  return res;
}
#endif
    
        
// Vector operator==(const Vector &V):
//    The == operator checks the two vectors are of the same size if VECTOR_CHECK is defined.
//    Then returns 1 if all the components of the two vectors are equal and 0 otherwise.

int 
Vector::operator==(const Vector &V) const
{
  if (sz != V.sz)
    return 0;

  double *dataThis = theData;
  double *dataV = V.theData;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != *dataV++)
      return 0;

  return 1;
}

int 
Vector::operator==(double value) const
{
  double *dataThis = theData;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != value)
      return 0;

  return 1;
}


// Vector operator!=(const Vector &V):
//  The != operator checks the two vectors are of the same size if VECTOR_CHECK is defined.
//  Then returns 1 if any of the components of the two vectors are unequal and 0 otherwise.

int
Vector::operator!=(const Vector &V) const
{
  if (sz != V.sz)
    return 1;

  double *dataThis = theData;
  double *dataV = V.theData;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != *dataV++)
      return 1;

  return 0;
}

int 
Vector::operator!=(double value) const
{
  double *dataThis = theData;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != value)
      return 1;

  return 0;
}


// friend OPS_Stream &operator<<(OPS_Stream &s, const Vector &V)
//        A function is defined to allow user to print the vectors using OPS_Streams.

OPS_Stream &operator<<(OPS_Stream &s, const Vector &V)
{
  return s.write(V.theData, V.sz);
}

// friend istream &operator>>(istream &s, Vector &V)
//        A function is defined to allow user to input the data into a Vector which has already
//        been constructed with data, i.e. Vector(int) or Vector(const Vector &) constructors.

/*
istream &operator>>(istream &s, Vector &V)
{
  for (int i=0; i<V.Size(); i++) 
    s >> V(i);
  
    return s;
}
*/


Vector operator*(double a, const Vector &V)
{
  return V * a;
}


int
Vector::Assemble(const Vector &V, int init_pos, double fact) 
{
  assert((init_pos >= 0) && ((init_pos + V.sz - 1) < sz));

  int cur_pos = init_pos;
  for (int j=0; j<V.sz; j++) 
     (*this)(cur_pos++) += V(j)*fact;

  return 0;
}

int
Vector::Extract(const Vector &V, int init_pos, double fact) 
{

  assert((init_pos >= 0) && ((init_pos + sz - 1) < V.sz));

  int cur_pos   = init_pos;  
  for (int j=0; j<sz; j++) 
     (*this)(j) = V(cur_pos++)*fact;

  return 0;
}

Matrix
Vector::operator%(const Vector &V) const
{
  // if sizes are compatable add
#ifdef _G3DEBUG
  if (sz != V.sz) {
    // else sizes are incompatable, do nothing but warning
    opserr <<  "WARNING Vector::tensor multiplication operator % - incompatable Vector sizes\n";
    return -1;
  }
#endif
  
  // we want result=a*b'
  
  Matrix result(sz,sz);
  for (int i=0; i<sz; i++)
    for (int j=0; j<sz; j++)
      result(i,j)=theData[i]*V.theData[j];
  
  return result;
}

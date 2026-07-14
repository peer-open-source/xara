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
// Purpose: This file contains the class implementation for Information.
//
// Written: fmk 10/99
// Revised:
//
#include <Information.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

Information::Information() 
  :theType(UnknownType),
   theID(0), theVector(0), theMatrix(0), theString(0)
{
  // does nothing
}

Information::Information(int val) 
  :theType(IntType), theInt(val),
  theID(0), theVector(0), theMatrix(0), theString(0)
{
  // does nothing
}

Information::Information(double val) 
  :theType(DoubleType), theDouble(val),
  theID(0), theVector(0), theMatrix(0), theString(0)
{
  // does nothing
}

Information::Information(const ID &val) 
  :theType(IdType),
  theID(0), theVector(0), theMatrix(0), theString(0)
{
  // Make a copy
  theID = new ID(val);
}

Information::Information(const Vector &val) 
  :theType(VectorType),
  theID(0), theVector(0), theMatrix(0), theString(0)
{
  // Make a copy
  // theVector = new Vector(val);
  if (val.Size() <= small_vector_size) {
    theVector.setData(small_vector, val.Size());
  } else {
    theVector.resize(val.Size());
  }
  theVector = val;
}

Information::Information(const Matrix &val) 
  :theType(MatrixType),
   theID(0), theVector(0), theMatrix(0), theString(0)
{
  // Make a copy
  theMatrix = new Matrix(val);
}

Information::Information(const ID &val1, const Vector &val2) 
  :theType(IdType),
   theID(0), theVector(0), theMatrix(0), theString(0)
{
  // Make a copy
  theID = new ID(val1);
#if 0
  theVector = new Vector(val2);
#else
  theVector.resize(val2.Size());
  theVector = val2;
#endif
}


Information::~Information() 
{
  if (theID != 0)
    delete theID;
  
  // if (theVector != 0)
  //   delete theVector;
  
  if (theMatrix != 0)
    delete theMatrix;
  
}

int 
Information::setInt(int newInt)
{
  theInt = newInt;
  return 0;
}

int 
Information::setDouble(double newDouble)
{
  theDouble = newDouble;
  return 0;
}

int 
Information::setID(const ID &newID)
{
  if (theID != nullptr) {
    *theID = newID;
  } else {
    theID = new ID(newID);
  }
  return 0;
}

int 
Information::setVector(const Vector &newVector)
{
#if 0
  if (theVector != nullptr) {
    *theVector = newVector;
  } else {
    theVector = new Vector(newVector);
  }
#else
  if (newVector.Size() != theVector.Size()) {
    theVector.resize(newVector.Size());
  }
  theVector = newVector;
#endif
  return 0;
}

int
Information::setVector(int i, double value)
{
  if (i < 0 || i >= theVector.Size())
    return -1;

  theVector(i) = value;
  return 0;
}

int
Information::setVectorAt(const Vector &newVector, int startIndex)
{
  if (startIndex < 0 || startIndex + newVector.Size() > theVector.Size())
    return -1;


  for (int i = 0; i < newVector.Size(); ++i)
    theVector(startIndex + i) = newVector(i);

  return 0;
}


int 
Information::setMatrix(const Matrix &newMatrix)
{
  if (theMatrix != 0) {
    *theMatrix = newMatrix;
  }
  else {
    theMatrix = new Matrix(newMatrix);
  }

  return 0;
}

int 
Information::setMatrix(int i, int j, double value)
{
  if (theMatrix == nullptr)
    return -1;

  theMatrix->operator()(i,j) = value;
  return 0;
}

int 
Information::setString(const char *newString)
{
  int newLength = int(strlen(newString));

  if (theString != 0) {
    int oldLength = int(strlen(theString));
    if (oldLength >= newLength) 
      strcpy(theString, newString);
    else {
      delete [] theString;
      theString = new char[newLength+1];
      strcpy(theString, newString);      
    }
  } else {
      theString = new char[newLength+1];
      strcpy(theString, newString);      
  }

  return 0;
}

void 
Information::Print(OPS_Stream &s, int flag)
{
  if (theType == IntType)
    s << theInt << " ";
  else if (theType == DoubleType)
    s << theDouble << " ";
  else if (theType == IdType && theID != 0)
    for (int i=0; i<theID->Size(); i++)
      s << (*theID)(i) << " ";
  else if (theType == VectorType && theVector != 0)
    for (int i=0; i<theVector.Size(); i++)
      s << (theVector)(i) << " ";
  else if (theType == MatrixType && theMatrix != 0) {
    for (int i=0; i<theMatrix->noRows(); i++) {
      for (int j=0; j<theMatrix->noCols(); j++)
        s <<  (*theMatrix)(i,j) << " ";
      s << "\n";
    }
  }
  return;
}


void 
Information::Print(std::ofstream &s, int flag)
{
  if (theType == IntType)
    s << theInt << " ";
  else if (theType == DoubleType)
    s << theDouble << " ";
  else if (theType == IdType && theID != 0)
    for (int i=0; i<theID->Size(); i++)
      s << (*theID)(i) << " ";
  else if (theType == VectorType && theVector != 0)
    for (int i=0; i<theVector.Size(); i++)
      s << (theVector)(i) << " ";
  else if (theType == MatrixType && theMatrix != 0) {
    for (int i=0; i<theMatrix->noRows(); i++) {
      for (int j=0; j<theMatrix->noCols(); j++)
	s <<  (*theMatrix)(i,j) << " ";
      s << "\n";
    }
  }
  return;
}

const Vector &
Information::getData() 
{
  if (theType == IntType) {
#if 0
    if (theVector == 0) 
      theVector = new Vector(1);
    (*theVector)(0) = theInt;
#else
    if (theVector.Size() != 1) {
      theVector.setData(small_vector, 1);
    }
    theVector(0) = theInt;
#endif
  } else if (theType == DoubleType) {
    #if 0
    if (theVector == 0) 
      theVector = new Vector(1);
    (*theVector)(0) = theDouble;
    #else
    if (theVector.Size() != 1) {
      theVector.setData(small_vector, 1);
    }
    theVector(0) = theDouble;
    #endif
  }
  else if (theType == IdType && theID != 0) {
#if 0
    if (theVector == 0) 
      theVector = new Vector(theID->Size());
    for (int i=0; i<theID->Size(); i++)
      (*theVector)(i) =  (*theID)(i);
#else
    if (theVector.Size() != theID->Size()) {
      theVector.resize(theID->Size());
    }
    for (int i=0; i<theID->Size(); i++)
      theVector(i) =  (*theID)(i);
#endif
  } else if (theType == MatrixType && theMatrix != 0) {
    int noRows = theMatrix->noRows();
    int noCols = theMatrix->noCols();
#if 0
    if (theVector == 0) 
      theVector = new Vector(noRows * noCols);
    int count = 0;
    for (int i=0; i<noRows; i++)
      for (int j=0; j<noCols; j++)
        (*theVector)(count++) = (*theMatrix)(i,j);
#else
    if (theVector.Size() != noRows * noCols) {
      theVector.resize(noRows * noCols);
    }
    int count = 0;
    for (int i=0; i<noRows; i++)
      for (int j=0; j<noCols; j++)
        theVector(count++) = (*theMatrix)(i,j);
#endif
  }
  
  return theVector;
}

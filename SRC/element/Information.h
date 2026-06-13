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
// Description: This file contains the class definition for Information.
// Information is a class in which all data members are public, i.e. basically
// a struct.
//
// Written: fmk 
// Created: 10/99
// Revision: A
//
#ifndef Information_h
#define Information_h

#include <OPS_Globals.h>
#include <fstream>

using std::ofstream;
class ID;
class Matrix;
class Vector;

enum InfoType {
  UnknownType, 
  IntType, 
  DoubleType, 
  IdType, 
  VectorType, 
  MatrixType
};

#define VIRTUAL
class Information
{
  public:
    Information();
    Information(int val);
    Information(double val);
    Information(const ID &val);
    Information(const Vector &val);
    Information(const Matrix &val);
    Information(const ID &val1, const Vector &val2);
    
    VIRTUAL ~Information();
    
    VIRTUAL int setInt(int newInt);
    VIRTUAL int setDouble(double newDouble);
    VIRTUAL int setID(const ID &newID);
    VIRTUAL int setVector(const Vector &newVector);
    VIRTUAL int setMatrix(const Matrix &newMatrix);
    VIRTUAL int setString(const char *theString);
    
    VIRTUAL void Print(OPS_Stream &s, int flag = 0);
    VIRTUAL void Print(ofstream &s, int flag = 0);
    VIRTUAL const Vector &getData(void);

    // data that is stored in the information object
    InfoType	theType;    // information about data type
    int		theInt;     // an integer value
    double	theDouble;  // a double value
    ID		*theID;     // pointer to an ID object, created elsewhere
    Vector 	*theVector; // pointer to a Vector object, created elsewhere
    Matrix	*theMatrix; // pointer to a Matrix object, created elsewhere
    char        *theString; // pointer to string

  private:

};

#endif


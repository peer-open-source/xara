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
#include <Vector.h>

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
    
    ~Information();
    
    int setInt(int newInt);
    int setDouble(double newDouble);
    int setID(const ID &newID);
    int setVector(const Vector &);
    int setVector(int i, double value);
    int setVectorAt(const Vector &, int start);
    int setMatrix(const Matrix &);
    int setMatrix(int i, int j, double value);
    int setString(const char *);

    int setMessage(const char *msg) {message = msg; return 0;}
    const char* getMessage() {return message;}
    
    void Print(OPS_Stream &s, int flag = 0);
    void Print(std::ofstream &s, int flag = 0);
    const Vector &getData();

    // data that is stored in the information object
    InfoType	theType;    // information about data type
    int		  theInt;     // an integer value
    double	theDouble;  // a double value
    ID		  *theID;     // pointer to an ID object, created elsewhere
    Vector 	 theVector; //
    Matrix	*theMatrix; //
    char    *theString; // pointer to string

  private:
    // small vector to avoid dynamic allocation
    static constexpr int small_vector_size = 12;
    double small_vector[small_vector_size];
    const char* message = "";
};

#endif


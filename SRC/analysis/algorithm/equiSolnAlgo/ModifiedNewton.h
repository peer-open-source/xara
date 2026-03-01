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
// Description: This file contains the class definition for 
// ModifiedNewton. ModifiedNewton is a class which performs a modified 
// Newton-Raphson solution algorithm in solving the equations.
// 
// Written: fmk 
// Created: 11/96 
// Revision: A 
//
#ifndef ModifiedNewton_h
#define ModifiedNewton_h

#include <EquiSolnAlgo.h>
#include <Vector.h>

class ConvergenceTest;

class ModifiedNewton: public EquiSolnAlgo
{
public:
  ModifiedNewton(int tangent, double iFactor = 0.0, double cFactor = 1.0);
  ~ModifiedNewton();

  int solveCurrentStep();    
  int getNumIterations() const override;

  virtual int sendSelf(int commitTag, Channel &);
  virtual int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

  void Print(OPS_Stream &, int flag) const final;
    
  protected:
    
  private:
    int tangent;
    int numIterations;

    double iFactor;
    double cFactor;
};

#endif



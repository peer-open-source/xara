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

// $Revision: 1.6 $
// $Date: 2007-04-02 23:41:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/NewtonRaphson.h,v $


#ifndef NewtonRaphson_h
#define NewtonRaphson_h

// File: ~/OOP/analysis/algorithm/NewtonRaphson.h
//
// Written: fmk
// Created: 11/96
// Revision: A
//

// Description: This file contains the class definition for
// NewtonRaphson. NewtonRaphson is a class which performs a Newton-Raphson
// solution algorithm in solving the equations.
//
// What: "@(#)NewtonRaphson.h, revA"

#include <EquiSolnAlgo.h>
#include <Vector.h>
#include <utility/XaraTimer.h>


class NewtonRaphson: public EquiSolnAlgo
{
  public:
  NewtonRaphson();
  NewtonRaphson(IncrementalIntegrator::TangentFlagType prediction_tangent,
                IncrementalIntegrator::TangentFlagType correction_tangent,
                double iFact,
                double cFact);
  ~NewtonRaphson();

  int solveCurrentStep();

  virtual int sendSelf(int commitTag, Channel &);
  virtual int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &, int flag) const final;

  int getNumIterations() const override;

 protected:


 private:
  IncrementalIntegrator::TangentFlagType
    correction_tangent, prediction_tangent;

  int numIterations;

  double iFactor;
  double cFactor;

  enum class Steps {
    Residual,
    Tangent,
    Solve,
    Update
  };

  Timer<Steps> timer;
};

#endif

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
// NewtonHallM. NewtonHallM is a class which performs a modified Newton-Raphson-Hall
// solution algorithm in solving the equations
//
// Written: fmk 
// Created: 03/18
// Revision: A 
//
#pragma once
#include <EquiSolnAlgo.h>
#include <Vector.h>

class NewtonHallM: public EquiSolnAlgo
{
public:
  NewtonHallM();
  NewtonHallM(double initFactor, int method, double alpha, double c);    
  ~NewtonHallM();
  
  int solveCurrentStep() override;

  void Print(OPS_Stream &, int flag) const final;    
  
  int getNumIterations() const override;

private:
  int numIterations;

  int method;
  double alpha;
  double c;

  double iFactor;
  double cFactor;
};


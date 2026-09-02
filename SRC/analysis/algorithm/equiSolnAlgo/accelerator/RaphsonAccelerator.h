/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
// Written: MHS
// Created: April 2002
//
// Description: This file contains the class definition for 
// RaphsonAccelerator. 
//
#pragma once
#include "Accelerator.h"
#include <IncrementalIntegrator.h>

class RaphsonAccelerator: public Accelerator
{
 public:
  RaphsonAccelerator(int tangent, double iFactor, double cFactor);
  virtual ~RaphsonAccelerator();
  
  int newStep(LinearSOE &theSOE);
  int accelerate(Vector &v, LinearSOE &theSOE, 
		 IncrementalIntegrator &theIntegrator);
  int updateTangent(IncrementalIntegrator &theIntegrator, bool& updated);
  bool updateTangent() override {return true;}

  int getTangent() override {return theTangent;}

  void Print(OPS_Stream &, int flag) const final;


 private:
  // Flag indicating which tangent to form
  int theTangent;
  double iFactor;  // Initial factor
  double cFactor;  // Current factor

  // Total number of iterations for the time step
  int totalIter;
};

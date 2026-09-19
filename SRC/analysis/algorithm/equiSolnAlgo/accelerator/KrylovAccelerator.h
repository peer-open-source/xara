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

// Description: This file contains the class definition for 
// KrylovAccelerator. 
#pragma once
#include "Accelerator.h"
#include <IncrementalIntegrator.h>

class KrylovAccelerator : public Accelerator
{
 public:
  KrylovAccelerator(int maxDim = 3, int tangent = CURRENT_TANGENT);
  virtual ~KrylovAccelerator();
  
  int  newStep(const LinearSOE &);
  int  accelerate(Vector &v, LinearSOE&, IncrementalIntegrator &);
  int  updateTangent(IncrementalIntegrator &, bool& factored);
  bool updateTangent();

  void Print(OPS_Stream &, int flag) const final;

 private:
  // Current dimension of Krylov subspace
  int dimension;
  
  // Size information
  int numEqns;
  int maxDimension;
  
  // Storage for update vectors
  Vector **v;
  
  // Storage for subspace vectors
  Vector **Av;
  
  // Array data sent to LAPACK subroutine
  double *AvData;
  double *rData;
  double *work;
  
  // Length of work array
  int lwork;
  // Which tangent to form at restart
  int theTangent;
};

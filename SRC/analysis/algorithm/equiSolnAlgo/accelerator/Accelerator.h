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
// Description: This file contains the class definition for 
// Accelerator. 
// Accelerator is an abstract base class. Its subclasses seek
// to find an improvement to the modified Newton prediction of
// the displacement increment.
//
// Written: MHS
// Created: April 2002
//
#pragma once
#include <Logging.h> // TODO
#include <IncrementalIntegrator.h>

class OPS_Stream;
class SolutionAlgorithm;
class LinearSOE;
class Vector;

class Accelerator
{
 public:
  Accelerator(int classTag);
  virtual ~Accelerator();
  
  // virtual functions
  virtual int newStep(LinearSOE &theSOE) = 0;
  virtual int accelerate(Vector &v, LinearSOE &, 
                         IncrementalIntegrator &) = 0;
  virtual int updateTangent(IncrementalIntegrator &theIntegrator, bool& factored);

  virtual bool updateTangent() {return false;}

  virtual IncrementalIntegrator::TangentFlagType getTangent() {return NO_TANGENT;}

  virtual void Print(OPS_Stream &s, int flag) const = 0;
  
};


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
#ifndef Accelerator_h
#define Accelerator_h

#include <MovableObject.h>
#include <OPS_Stream.h>
#include <Logging.h> // TODO
#include <IncrementalIntegrator.h>

class SolutionAlgorithm;
class LinearSOE;
class Vector;

class Accelerator : public MovableObject
{
 public:
  Accelerator(int classTag);
  virtual ~Accelerator();
  
  // virtual functions
  virtual int newStep(LinearSOE &theSOE) = 0;
  virtual int accelerate(Vector &v, LinearSOE &theSOE, 
                         IncrementalIntegrator &theIntegrator) = 0;
  virtual int updateTangent(IncrementalIntegrator &theIntegrator, bool& factored);
  virtual bool updateTangent(void) {return false;}

  virtual int getTangent(void) {return NO_TANGENT;}

  virtual void Print(OPS_Stream &s, int flag=0) = 0;
  
};

#endif

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
// Description: This file contains the class implementation for 
// PeriodicAccelerator. 

#include <PeriodicAccelerator.h>

#include <Vector.h>
#include <LinearSOE.h>
#include <IncrementalIntegrator.h>

#include <ID.h>

PeriodicAccelerator::PeriodicAccelerator(int iter, int tangent)
  :Accelerator(),
   iteration(0), totalIter(0), maxIter(iter), theTangent(tangent)
{
  if (maxIter < 1)
    maxIter = 1;
}

PeriodicAccelerator::~PeriodicAccelerator()
{

}

int 
PeriodicAccelerator::newStep(const LinearSOE &theSOE)
{
  totalIter = 0;

  // Reset iteration counter
  iteration = (theTangent == CURRENT_TANGENT) ? maxIter : 0;

  return 0;
}

int
PeriodicAccelerator::accelerate(Vector &vStar, LinearSOE &theSOE,
				IncrementalIntegrator &theIntegrator)
{
  iteration++;
  totalIter++;

  return 0; 
}

int
PeriodicAccelerator::updateTangent(IncrementalIntegrator &theIntegrator, bool& factored)
{
  /*
  if (theTangent == NO_TANGENT)
    return 0;

  else if (theTangent == SECOND_TANGENT) {
    if (totalIter == maxIter) {
      theIntegrator.formTangent(CURRENT_TANGENT);
      return 1;
    }
    else
      return 0;
  }

  else { // CURRENT_TANGENT or INITIAL_TANGENT
    if (iteration >= maxIter) {
      iteration = 0;
      theIntegrator.formTangent(theTangent);
      if (theTangent == CURRENT_TANGENT)
	return 1;
      else
	return 0;
    }
    else
      return 0;
  }
  */
  factored = false;
  if (iteration < maxIter)
    return 0;

  switch (theTangent) {
  case CURRENT_TANGENT:
    iteration = 0;
    theIntegrator.formTangent(CURRENT_TANGENT);
    return 1;
    break;
  case INITIAL_TANGENT:
    iteration = 0;
    theIntegrator.formTangent(INITIAL_TANGENT);
    return 0;
    break;
  case NO_TANGENT:
    iteration = 0;
    return 0;
    break;
  default:
    return 0;
  }
}

bool
PeriodicAccelerator::updateTangent()
{
  if (iteration > maxIter) {
    iteration = 0;
    return true;
  }
  else 
    return false;
}

void
PeriodicAccelerator::Print(OPS_Stream &s, int flag) const
{
  s << "PeriodicAccelerator" << "\n";
  s << "\tIterations till restart: " << maxIter << "\n";
}

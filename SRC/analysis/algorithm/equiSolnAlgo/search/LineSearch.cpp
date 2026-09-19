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
// Description: This file contains the class implementation for 
// LineSearch. 
//
// Written: fmk 
// Created: 11/01
//
#include <LineSearch.h>
#include <ConvergenceTest.h>
#include <IncrementalResidual.h>
#include <LinearSOE.h>
#include <Vector.h>


LineSearch::LineSearch()
: Gn(0), dXs(0)
{

}

LineSearch::~LineSearch()
{

}


int
LineSearch::apply(IncrementalResidual &theIntegrator, 
                  LinearSOE &theSOE, 
                  ConvergenceTest &theTest,
                  Vector &dX, 
                  Vector &Go)
{
  Gn.resize(Go.Size());
  dXs.resize(dX.Size());

  // initial value of s
  const double s0 = dX ^ Go;

  Gn.Zero();
  if (theIntegrator.formUnbalance(Gn) < 0)
    return -1;

  
  int status = 0;
  theTest.start(theSOE);
  if (theTest.test(Gn, dX) < 1) {
    const double s = dX ^ Gn;
    this->newStep(Go);
    status = this->search(s0, s, dX, Go, dXs, theIntegrator);
    if (status == 0) {
      // Search was successful, update dX with the new value
      dX = dXs;
    }
  }
  return status;
}


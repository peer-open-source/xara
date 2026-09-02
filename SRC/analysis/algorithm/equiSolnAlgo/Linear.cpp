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
// Linear. Linear is a class which performs a linear solution algorithm
// to solve the equations.
// 
// File: ~/OOP/analysis/algorithm/Linear.C 
// 
// Written: fmk 
// Created: Sun Sept 15 15:06:47: 1996 
//
#include <Linear.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>

// Constructor
Linear::Linear(int theTangent, int Fact)
  : EquiSolnAlgo(EquiALGORITHM_TAGS_Linear), incrTangent(theTangent), factorOnce(Fact)
{

}

// Destructor
Linear::~Linear()
{

}

// Perform the linear solution algorithm.

int 
Linear::solveCurrentStep()
{
  // set up some pointers and check they are valid
  // NOTE this could be taken away if we set Ptrs as protecetd in superclass

  LinearSOE  *theSOE = this->getLinearSOEptr();
  IncrementalIntegrator  *theIncIntegrator = this->getIncrementalIntegratorPtr(); 

  if ((theIncIntegrator ==0 ) || (theSOE == 0)) {
    return SolutionAlgorithm::BadAlgorithm;
  }

  if (factorOnce != 2) {
    if (theIncIntegrator->formTangent(incrTangent) < 0)
      return SolutionAlgorithm::BadFormTangent;

    if (factorOnce == 1)
      factorOnce = 2;
  }


  if (theIncIntegrator->formUnbalance() < 0)
    return SolutionAlgorithm::BadFormResidual;


  if (theSOE->solve() < 0)
    return SolutionAlgorithm::BadLinearSolve;


  if (theIncIntegrator->update(theSOE->getX()) < 0)
    return SolutionAlgorithm::BadStepUpdate;

  return 0;
}


void
Linear::Print(OPS_Stream &s, int flag) const
{
  s << "\t Linear algorithm";
}

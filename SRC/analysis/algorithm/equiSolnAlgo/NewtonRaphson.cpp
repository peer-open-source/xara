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
// NewtonRaphson. NewtonRaphson is a class which uses the
// Newton-Raphson solution algorithm
// to solve the equations. No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.
// 
// What: "@(#)NewtonRaphson.C, revA"
// 
// Written: fmk 
// Created: Sun Sept 15 15:06:47: 1996 
//
#include <NewtonRaphson.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>


NewtonRaphson::NewtonRaphson(IncrementalIntegrator::TangentFlagType prediction_tangent, 
                             IncrementalIntegrator::TangentFlagType correction_tangent,
                            double iFact, 
                            double cFact)
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
 prediction_tangent(prediction_tangent), 
 correction_tangent(correction_tangent),
 iFactor(iFact), cFactor(cFact)
{

}

NewtonRaphson::NewtonRaphson()
    :EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
    prediction_tangent(CURRENT_TANGENT), 
    correction_tangent(CURRENT_TANGENT),
    iFactor(0.), cFactor(1.)
{

}



NewtonRaphson::~NewtonRaphson()
{
  
}


int 
NewtonRaphson::solveCurrentStep()
{
  // set up some pointers and check they are valid
  // NOTE this could be taken away if we set Ptrs as protecetd in superclass
  IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
  LinearSOE  *theSOE = this->getLinearSOEptr();

  if  (  (theIntegrator== nullptr)
      || (theSOE       == nullptr)
      || (theTest      == nullptr)) {
      opserr << "WARNING NewtonRaphson::solveCurrentStep() - setLinks() has";
      opserr << " not been called - or no ConvergenceTest has been set\n";
      return SolutionAlgorithm::BadAlgorithm;
  }        

  //
  // 1 Form unbalance
  //
  if (theIntegrator->formUnbalance() < 0)
    return SolutionAlgorithm::BadFormResidual;


  // Its prbably good to pass theTest as an argument to solveCurrentStep.
  if (theTest->start(*theSOE) < 0)
    return SolutionAlgorithm::BadTestStart;


  int result = ConvergenceTest::Continue;

  numIterations = 0;

  //
  // 2 Corrections
  //
  do {
    //
    // 2.1 Form tangent
    //

    if (numIterations == 0) {
      SOLUTION_ALGORITHM_tangentFlag = prediction_tangent;
      if (theIntegrator->formTangent(prediction_tangent) < 0)
        return SolutionAlgorithm::BadFormTangent;
    }
    else {
      SOLUTION_ALGORITHM_tangentFlag = correction_tangent;
      if (theIntegrator->formTangent(correction_tangent, iFactor, cFactor) < 0)
        return SolutionAlgorithm::BadFormTangent;
    }

    //
    // 2.2 Solve for dx
    //
    if (theSOE->solve() < 0) 
      return SolutionAlgorithm::BadLinearSolve;

    if (theIntegrator->update(theSOE->getX()) < 0)
      return SolutionAlgorithm::BadStepUpdate;

    //
    // 2.3 Form updated residual
    //

    if (theIntegrator->formUnbalance() < 0)
      return SolutionAlgorithm::BadFormResidual;

    //
    // 2.4 Test on updated residual
    //
    result = theTest->test(*theSOE);
    numIterations++;
    this->record(numIterations);

  }  while (result == ConvergenceTest::Continue);

  if (result == ConvergenceTest::Failure)
    return SolutionAlgorithm::TestFailed;

  // if postive result, we are returning what the convergence test
  // returned which should be the number of iterations  
  return result;
}


int
NewtonRaphson::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(3);
  data(0) = double(correction_tangent);
  data(1) = double(prediction_tangent);
  data(2) = iFactor;
  data(3) = cFactor;
  return theChannel.sendVector(this->getDbTag(), cTag, data);
}

int
NewtonRaphson::recvSelf(int cTag, 
                        Channel &theChannel, 
                        FEM_ObjectBroker &theBroker)
{
  static Vector data(3);
  theChannel.recvVector(this->getDbTag(), cTag, data);
  correction_tangent = int(data(0));
  prediction_tangent = int(data(1));
  iFactor = data(2);
  cFactor = data(3);
  return 0;
}


void
NewtonRaphson::Print(OPS_Stream &s, int flag) const
{
  if (flag == 0) {
    s << "NewtonRaphson" << endln;
  }
}


int
NewtonRaphson::getNumIterations() const
{
  return numIterations;
}



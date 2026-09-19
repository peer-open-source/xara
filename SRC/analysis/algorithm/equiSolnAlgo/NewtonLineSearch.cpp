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
// Description: This file contains the implementation for NewtonLineSearch. 
//
// Written: fmk 
// Created: 11/96 
// Modified: Ed "C++" Love 10/00 to perform the line search
//
#include <NewtonLineSearch.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <ConvergenceTest.h>
#include <LineSearch.h>
#include <update/UpdateBFGS.h>
#include <ID.h>


NewtonLineSearch::NewtonLineSearch(LineSearch *theSearch, 
                                   IncrementalIntegrator::TangentFlagType prediction_tangent,
                                   IncrementalIntegrator::TangentFlagType correction_tangent) 
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonLineSearch),
 theLineSearch(theSearch),
 prediction_tangent(prediction_tangent),
 correction_tangent(correction_tangent)
{

}


NewtonLineSearch::~NewtonLineSearch()
{
  if (theLineSearch != nullptr)
    delete theLineSearch;
}


int 
NewtonLineSearch::solveCurrentStep()
{
  // set up some pointers and check they are valid
  IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
  LinearSOE  *theSOE = this->getLinearSOEptr();

  if ((theIntegrator == 0) || (theSOE == 0) || (theTest == 0)) {
    return SolutionAlgorithm::BadAlgorithm;
  }
  Go.resize(theSOE->getNumEqn());
  Gn.resize(theSOE->getNumEqn());
  dX.resize(theSOE->getNumEqn());
  dXs.resize(theSOE->getNumEqn());


  ConvergenceTest *theOtherTest = nullptr;
  theOtherTest = theTest->getCopy(10);

  //
  // 1 Form unbalance
  //
  Go.Zero();
  if (theIntegrator->formUnbalance(Go) < 0) 
    return SolutionAlgorithm::BadFormResidual;
  theSOE->setB(Go);

  //
  //
  //
  if (theTest->start(*theSOE) < 0)
    return SolutionAlgorithm::BadTestStart;

  int result = ConvergenceTest::Continue;
  int numIterations = 0;
  do {

    // 2.1 Form the tangent
    if (numIterations == 0) {
      SOLUTION_ALGORITHM_tangentFlag = prediction_tangent;
      if (theIntegrator->formTangent(prediction_tangent) < 0)
        return SolutionAlgorithm::BadFormTangent;
    }
    else if (correction_tangent == PREDICTOR_TANGENT) {
      // here we reuse the tangent formed at the first iteration, relying
      // on it being maintained in the LinearSOE.
      ;
    }
    else {
      SOLUTION_ALGORITHM_tangentFlag = correction_tangent;
      if (theIntegrator->formTangent(correction_tangent) < 0)
        return SolutionAlgorithm::BadFormTangent;
    }

    //
    // 2.2 Solve for dx
    //
    if (theSOE->solve(Go, dX) < 0)
      return SolutionAlgorithm::BadLinearSolve;      


    if (theIntegrator->update(dX) < 0)     
      return SolutionAlgorithm::BadStepUpdate;

    //
    // line search 
    //
    dXs = dX;
    if (theLineSearch != nullptr) {
      // initial value of s
      double so = dX ^ Go;

      Gn.Zero();
      if (theIntegrator->formUnbalance(Gn) < 0)
        return SolutionAlgorithm::BadFormResidual;

      // do a line search only if convergence criteria not met
      theOtherTest->start(*theSOE);
      result = theOtherTest->test(Gn, dX);

      if (result < 1) {

        // new value of s
        double su =  dX ^ Gn;
    
        int search_result = 0;
        theLineSearch->newStep(Go);
        search_result = theLineSearch->search(so, su, dX, Gn, dXs, *theIntegrator);

        if (search_result < 0)
          return search_result;
      }
    }

    this->record(0);

    result = theTest->test(Gn, dXs);
    numIterations++;
    if (result == ConvergenceTest::Continue)
      Go = Gn;

  } while (result == ConvergenceTest::Continue);

  if (result == ConvergenceTest::Failure)
    return SolutionAlgorithm::TestFailed;


  // note - if positive result we are returning what the convergence test returned
  // which should be the number of iterations
  return result;
}


void
NewtonLineSearch::Print(OPS_Stream &s, int flag) const
{
  if (flag == 0) 
    s << "NewtonLineSearch\n";

  if (theLineSearch != 0)
    theLineSearch->Print(s, flag);
}


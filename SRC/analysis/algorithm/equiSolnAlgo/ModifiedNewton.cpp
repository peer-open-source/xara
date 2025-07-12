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
// ModifiedNewton. ModifiedNewton is a class which uses the
// Newton-Raphson solution algorithm
// to solve the equations. No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.
// 
// Written: fmk 
// Created: 11/96 
//
#include <ModifiedNewton.h>
#include <AnalysisModel.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>


ModifiedNewton::ModifiedNewton(int theTangentToUse, double iFact, double cFact)
:EquiSolnAlgo(EquiALGORITHM_TAGS_ModifiedNewton),
 tangent(theTangentToUse), iFactor(iFact), cFactor(cFact)
{
  
}


ModifiedNewton::~ModifiedNewton()
{

}


int 
ModifiedNewton::solveCurrentStep()
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass
    AnalysisModel       *theAnalysisModel = this->getAnalysisModelPtr();
    IncrementalIntegrator *theIncIntegratorr = this->getIncrementalIntegratorPtr();
    LinearSOE                *theSOE = this->getLinearSOEptr();

    if ((theAnalysisModel == nullptr) 
        || (theIncIntegratorr == nullptr) 
        || (theSOE == nullptr)
        || (theTest == nullptr)) {
      opserr << "WARNING ModifiedNewton::solveCurrentStep() - setLinks() has";
      opserr << " not been called - or no ConvergenceTest has been set\n";
      return SolutionAlgorithm::BadAlgorithm;
    }

    if (theIncIntegratorr->formUnbalance() < 0) {
      opserr << "WARNING ModifiedNewton::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance()\n";        
      return SolutionAlgorithm::BadFormResidual;
    }        

    SOLUTION_ALGORITHM_tangentFlag = tangent;
    if (theIncIntegratorr->formTangent(tangent, iFactor, cFactor) < 0)
      return SolutionAlgorithm::BadFormTangent;


    // set itself as the ConvergenceTest objects EquiSolnAlgo
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      opserr << "ModifiedNewton::solveCurrentStep() -";
      opserr << "the ConvergenceTest object failed in start()\n";
      return SolutionAlgorithm::BadTestStart;
    }

    // repeat until convergence is obtained or reach max num iterations
    int result = -1;
    numIterations = 0;
    do {
      if (theSOE->solve() < 0)
        return SolutionAlgorithm::BadLinearSolve;
      
      if (theIncIntegratorr->update(theSOE->getX()) < 0)
        return SolutionAlgorithm::BadStepUpdate;

      if (theIncIntegratorr->formUnbalance() < 0)
        return SolutionAlgorithm::BadFormResidual;

      result = theTest->test();
      numIterations++;
      this->record(numIterations);

    }  while (result == ConvergenceTest::Continue);

    if (result == ConvergenceTest::Failure)
      return SolutionAlgorithm::TestFailed;

    return result;
}

int
ModifiedNewton::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(3);
  data(0) = tangent;
  data(1) = iFactor;
  data(2) = cFactor;
  return theChannel.sendVector(this->getDbTag(), cTag, data);
}

int
ModifiedNewton::recvSelf(int cTag, 
                        Channel &theChannel, 
                        FEM_ObjectBroker &theBroker)
{
  static Vector data(3);
  theChannel.recvVector(this->getDbTag(), cTag, data);
  tangent = data(0);
  iFactor = data(1);
  cFactor = data(2);
  return 0;
}

void
ModifiedNewton::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) {
        s << "ModifiedNewton";
    }
}

int
ModifiedNewton::getNumIterations()
{
  return numIterations;
}

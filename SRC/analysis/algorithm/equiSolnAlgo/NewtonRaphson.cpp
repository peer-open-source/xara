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
#include <AnalysisModel.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>

#if 1
#include <elementAPI.h>

void *
OPS_ADD_RUNTIME_VPV(OPS_NewtonRaphsonAlgorithm)
{
    int formTangent = CURRENT_TANGENT;
    double iFactor = 0;
    double cFactor = 1;

    while (OPS_GetNumRemainingInputArgs() > 0) {
      const char* type = OPS_GetString();
      if (strcmp(type,"-secant")==0 || 
          strcmp(type,"-Secant")==0) {
        formTangent = CURRENT_SECANT;
        iFactor = 0;
        cFactor = 1.0;

      } else if (strcmp(type,"-initial")==0 || 
                 strcmp(type,"-Initial")==0) {
        formTangent = INITIAL_TANGENT;
        iFactor = 1.;
        cFactor = 0;

      } else if (strcmp(type,"-intialThenCurrent")==0 || 
                 strcmp(type,"-intialCurrent")==0) {
        formTangent = INITIAL_THEN_CURRENT_TANGENT;
        iFactor = 0;
        cFactor = 1.0;

      } else if (strcmp(type,"-hall")==0 || 
                 strcmp(type,"-Hall")==0) {

        formTangent = HALL_TANGENT;
        iFactor = 0.1;
        cFactor = 0.9;
        if (OPS_GetNumRemainingInputArgs() == 2) {
          double data[2];
          int numData = 2;
          if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
            opserr << "WARNING invalid data reading 2 hall factors\n";
            return 0;
          }
          iFactor = data[0];
          cFactor = data[1];
        }
      }
    }

    return new NewtonRaphson(formTangent, iFactor, cFactor);
}
#endif 

// Constructor
NewtonRaphson::NewtonRaphson(int theTangentToUse, double iFact, double cFact)
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
 tangent(theTangentToUse), iFactor(iFact), cFactor(cFact)
{

}

NewtonRaphson::NewtonRaphson()
        :EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
        tangent(CURRENT_TANGENT), iFactor(0.), cFactor(1.)
{

}


NewtonRaphson::NewtonRaphson(ConvergenceTest &theT, int theTangentToUse, double iFact, double cFact)
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonRaphson),
 tangent(theTangentToUse), iFactor(iFact), cFactor(cFact)
{

}

// Destructor
NewtonRaphson::~NewtonRaphson()
{
  
}


int 
NewtonRaphson::solveCurrentStep()
{
    // set up some pointers and check they are valid
    // NOTE this could be taken away if we set Ptrs as protecetd in superclass
    AnalysisModel   *theAnaModel = this->getAnalysisModelPtr();
    IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
    LinearSOE  *theSOE = this->getLinearSOEptr();

    if  ( (theAnaModel  == nullptr) 
       || (theIntegrator== nullptr)
       || (theSOE       == nullptr)
       || (theTest      == nullptr)){
        opserr << "WARNING NewtonRaphson::solveCurrentStep() - setLinks() has";
        opserr << " not been called - or no ConvergenceTest has been set\n";
        return SolutionAlgorithm::BadAlgorithm;
    }        

    //
    // 1 Form unbalance
    //
    if (theIntegrator->formUnbalance() < 0)
      return SolutionAlgorithm::BadFormResidual;


    // Set self as the ConvergenceTest's EquiSolnAlgo;
    // TODO: Perhaps move this to the BasicAnalysisBuilder;
    // at least the call to start. Its prbably good to pass
    // theTest as an argument to solveCurrentStep.
    theTest->setEquiSolnAlgo(*this);
    if (theTest->start() < 0) {
      return SolutionAlgorithm::BadTestStart;
    }

    int result = ConvergenceTest::Continue;

    numIterations = 0;

    //
    // 2 Corrections
    //
    do {
      //
      // 2.1 Form tangent
      //
      if (tangent == INITIAL_THEN_CURRENT_TANGENT) {

        if (numIterations == 0) {

          SOLUTION_ALGORITHM_tangentFlag = INITIAL_TANGENT;
          if (theIntegrator->formTangent(INITIAL_TANGENT) < 0)
            return SolutionAlgorithm::BadFormTangent;

        } else {
          SOLUTION_ALGORITHM_tangentFlag = CURRENT_TANGENT;
          if (theIntegrator->formTangent(CURRENT_TANGENT) < 0)
            return SolutionAlgorithm::BadFormTangent;
        }

      } else {

        SOLUTION_ALGORITHM_tangentFlag = tangent;
        if (theIntegrator->formTangent(tangent, iFactor, cFactor) < 0)
          return SolutionAlgorithm::BadFormTangent;

      }

      //
      // 2.2 Solve for dx
      //
      if (theSOE->solve() < 0) 
        return SolutionAlgorithm::BadLinearSolve;

      //
      // 2.3 Form updated residual
      //
      if (theIntegrator->update(theSOE->getX()) < 0)
        return SolutionAlgorithm::BadStepUpdate;

      if (theIntegrator->formUnbalance() < 0)
        return SolutionAlgorithm::BadFormResidual;

      //
      // 2.4 Test on updated residual
      //
      result = theTest->test();
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
  data(0) = tangent;
  data(1) = iFactor;
  data(2) = cFactor;
  return theChannel.sendVector(this->getDbTag(), cTag, data);
}

int
NewtonRaphson::recvSelf(int cTag, 
                        Channel &theChannel, 
                        FEM_ObjectBroker &theBroker)
{
  static Vector data(3);
  theChannel.recvVector(this->getDbTag(), cTag, data);
  tangent = int(data(0));
  iFactor = data(1);
  cFactor = data(2);
  return 0;
}


void
NewtonRaphson::Print(OPS_Stream &s, int flag)
{
  if (flag == 0) {
    s << "NewtonRaphson" << endln;
  }
}


int
NewtonRaphson::getNumIterations(void)
{
  return numIterations;
}



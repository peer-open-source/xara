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
                                                                        
// Written: fmk 03/18

// Description: This file contains the class definition for 
// NewtonHallM. 

#include <NewtonHallM.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <ID.h>

#include <math.h>

#if 1
#include <elementAPI.h>
void *
OPS_ADD_RUNTIME_VPV(OPS_NewtonHallM)
{
  int method = 0;
  double iFactor = .1;
  double alpha = .01;
  double c = 100;
  double data[2];
  
  int numData = 1;
  if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
    opserr << "WARNING invalid data reading 2 hall factors\n";
    return 0;
  }
  iFactor = data[0];

  while (OPS_GetNumRemainingInputArgs() > 0) {
      const char* type = OPS_GetString();
      if(strcmp(type,"-exp")==0 || strcmp(type,"-Exp")==0) {
        numData = 1;
        if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
          opserr << "WARNING invalid data reading 2 hall factors\n";
          return 0;
        } else 
          alpha = data[0];
      } else if(strcmp(type,"-sigmoid")==0 || strcmp(type,"-Sigmoid")==0) {
        method = 1;
        int numData = 2;
        if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
          opserr << "WARNING invalid data reading 2 hall factors\n";
          return 0;
        } else {
          alpha = data[0];
          c = data[1];
        }
      } else if(strcmp(type,"-constant")==0 || strcmp(type,"-Constant")==0) {
        method = 2;
        int numData = 1;
        if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
          opserr << "WARNING invalid data reading 2 hall factors\n";
          return 0;
        } else {
          c = data[0];
        }
      }
    }

  return new NewtonHallM(iFactor, method, alpha, c);
}
#endif 

// Constructor
NewtonHallM::NewtonHallM(double initFactor, int mthd, double alphaFact, double cFact)
:EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonHallM),
 numIterations(0), method(mthd), alpha(alphaFact), c(cFact), iFactor(initFactor), cFactor(cFact)
{

}

NewtonHallM::NewtonHallM()
  :EquiSolnAlgo(EquiALGORITHM_TAGS_NewtonHallM),
   numIterations(0), method(0), alpha(.1), c(0), iFactor(0.1)
{

}

NewtonHallM::~NewtonHallM()
{
  

}


int 
NewtonHallM::solveCurrentStep()
{
  // set up some pointers and check they are valid
  // NOTE this could be taken away if we set Ptrs as protecetd in superclass
  IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
  //IncrementalIntegrator *theIntegratorSens=this->getIncrementalIntegratorPtr();//Abbas
  LinearSOE  *theSOE = this->getLinearSOEptr();

  if ((theIntegrator == 0) || (theSOE == 0)  || (theTest == 0)){
      return -5;
  }        

  if (theIntegrator->formUnbalance() < 0) {      
    return -2;
  }            

  if (theTest->start(*theSOE) < 0) {
    return SolutionAlgorithm::BadTestStart;
  }

  int result = -1;
  numIterations = 0;

  do {

    int tangent = HALL_TANGENT;
    SOLUTION_ALGORITHM_tangentFlag = tangent;

    double iFact, cFact;
    if (method == 0) {
      iFact = iFactor*exp(-alpha*numIterations);
      cFact = 1.0 - iFact;
    } else if (method ==1) {
      double iFact0 = 1.0/(1. + exp(-alpha*c));
      iFact = 1/(1 + exp(alpha*(numIterations-c)));
      iFact = iFactor*iFact/iFact0;
      cFact = 1.0 - iFact;
    } else {
      iFact = iFactor;
      cFact = cFactor;
    }
    
    if (theIntegrator->formTangent(tangent, iFact, cFact) < 0){
      opserr << "WARNING NewtonHallM::solveCurrentStep() -";
      opserr << "the Integrator failed in formTangent()\n";
      return -1;
    }                    
    
    if (theSOE->solve() < 0) {
      opserr << "WARNING NewtonHallM::solveCurrentStep() -";
      opserr << "the LinearSysOfEqn failed in solve()\n";        
      return -3;
    }            
    
    if (theIntegrator->update(theSOE->getX()) < 0) {
      opserr << "WARNING NewtonHallM::solveCurrentStep() -";
      opserr << "the Integrator failed in update()\n";        
      return -4;
    }                
    if (theIntegrator->formUnbalance() < 0) {
      opserr << "WARNING NewtonHallM::solveCurrentStep() -";
      opserr << "the Integrator failed in formUnbalance()\n";        
      return -2;
    }        
    
    result = theTest->test(*theSOE);
    numIterations++;
    this->record(numIterations);
    
    
  }  while (result == ConvergenceTest::Continue);
  
  if (result == ConvergenceTest::Failure) {
    opserr << "NewtnRaphson::solveCurrentStep() -";
    opserr << "the ConvergenceTest object failed in test()\n";
    return -3;
  }
  
  return result;
}


void
NewtonHallM::Print(OPS_Stream &s, int flag) const
{
  if (flag == 0) {
    s << "NewtonHallM" << "\n";
    if (method == 0)
      s << "  -exp method with alpha = " << alpha << "\n";
    else
      s << "  -sigmoid method with alpha: " << alpha << " c: " << c << "\n";
  }
}


int
NewtonHallM::getNumIterations() const
{
  return numIterations;
}



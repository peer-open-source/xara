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
// Description: This file contains the class definition for AcceleratedNewton.  
// AcceleratedNewton is a class which uses a Krylov 
// subspace accelerator on the modified Newton method. 
//
// The accelerator is described by Carlson and Miller. 
//
// [1] Carlson and Miller "Design and Application of a 1D GWMFE Code"
//     from SIAM Journal of Scientific Computing (Vol. 19, No. 3,
//     pp. 728-765, May 1998)
//
// Written: MHS
// Created: Oct 2001
//
#include <AcceleratedNewton.h>
#include <Accelerator.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <LineSearch.h>
#include <fstream>



AcceleratedNewton::AcceleratedNewton(Accelerator *theAccel,
                                     int incr_tangent)
  : EquiSolnAlgo(EquiALGORITHM_TAGS_AcceleratedNewton),
    tangent(incr_tangent),
    search(nullptr),
    theAccelerator(theAccel), vAccel(0), 
    numFactorizations(0), numIterations(0)
{
 
}


AcceleratedNewton::~AcceleratedNewton()
{
  if (theAccelerator != nullptr)
    delete theAccelerator;

  if (vAccel != 0)
    delete vAccel;
}


int 
AcceleratedNewton::solveCurrentStep()
{
  // set up some pointers and check they are valid
  // NOTE this could be taken away if we set Ptrs as protecetd in superclass
  IncrementalIntegrator *theIntegrator = this->getIncrementalIntegratorPtr();
  LinearSOE *theSOE = this->getLinearSOEptr();
  
  if ((theIntegrator == 0) || (theSOE == 0)  || (theTest == 0)) {
    return SolutionAlgorithm::BadAlgorithm;
  }        

  // Set up memory in the accelerator
  if (theAccelerator != nullptr)
    theAccelerator->newStep(*theSOE);

  int numEqns = theSOE->getNumEqn();

  if (vAccel == nullptr)
    vAccel = new Vector(numEqns);

  if (vAccel->Size() != numEqns) {
    vAccel->resize(numEqns);
  }

  //
  // 1. Form unbalance
  //
  // Evaluate system residual R(y_0)
  if (theIntegrator->formUnbalance() < 0) {
    return SolutionAlgorithm::BadFormResidual;
  }

  // Evaluate system Jacobian J = R'(y)|y_0
  if (theIntegrator->formTangent(tangent) < 0) {
    return SolutionAlgorithm::BadFormTangent;
  }
  // Count factorization of the first tangent
  numFactorizations++;


  if (theTest->start(*theSOE) < 0) {
    return SolutionAlgorithm::BadTestStart;
  }

  // Loop counter
  int k = 1;

  int result = ConvergenceTest::Continue;

  numIterations = 0;

  do {

    // Solve for displacement increment
    if (theSOE->solve() < 0) 
      return SolutionAlgorithm::BadLinearSolve;

    // Get the modified Newton increment
    *vAccel = theSOE->getX();

    // Accelerate the displacement increment
    if (theAccelerator != nullptr) {
      if (theAccelerator->accelerate(*vAccel, *theSOE, *theIntegrator) < 0) {
        opserr << "the Accelerator failed\n";
        return -1;
      }
    }

    // if (search != nullptr) {
    //   const double s0 = - (*vAccel ^ theSOE->getB());
    //   if (search->search(*vAccel, *theSOE, *theIntegrator) < 0) {
    //     return -1;// SolutionAlgorithm::BadLineSearch;
    //   }
    // }

    // Update system with accelerated displacement increment v_{k+1}
    if (theIntegrator->update(*vAccel) < 0) {
      return SolutionAlgorithm::BadStepUpdate;
    }        

    // Evaluate residual
    if (theIntegrator->formUnbalance() < 0) 
      return SolutionAlgorithm::BadFormResidual;

    numIterations++;

    // Check convergence criteria
    result = theTest->test(*theSOE);

    if (result == ConvergenceTest::Continue) {
      // Let the accelerator update the tangent if needed
      if (theAccelerator != nullptr) {
        bool did_update = false;
        int ret = theAccelerator->updateTangent(*theIntegrator, did_update);
        if (ret < 0) {
          opserr << "the Accelerator failed to update tangent\n";
          return -1;
        }
        if (did_update)
          numFactorizations++;
      }
    }
    this->record(k++);

  } while (result == ConvergenceTest::Continue);
 

  if (result == ConvergenceTest::Failure) 
    return SolutionAlgorithm::TestFailed;
  
  // note - if positive result we are returning what the convergence
  // test returned which should be the number of iterations
  return result;
}


void
AcceleratedNewton::Print(OPS_Stream &s, int flag) const
{
  s << "QuasiNewton" << "\n";
  LinearSOE *theSOE = this->getLinearSOEptr();
  s << "\tNumber of equations: " << theSOE->getNumEqn() << "\n";

  if (theAccelerator != 0)
    theAccelerator->Print(s,flag);
  else
    s << "\tNo accelerator --> Modified Newton" << "\n";
}


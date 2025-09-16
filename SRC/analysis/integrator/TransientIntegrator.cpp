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
// Description: This file contains the class definition for TransientIntegrator.
// TransientIntegrator is an algorithmic class for setting up the finite element
// equations for a static analysis and for Incrementing the nodal displacements
// with the values in the soln vector to the LinearSOE object. 
//
// Written: fmk 
// Created: Tue Sept 17 15:54:47: 1996
// Revision: A
//
#include <TransientIntegrator.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>

TransientIntegrator::TransientIntegrator(int clasTag)
 : IncrementalIntegrator(clasTag)
{

}

TransientIntegrator::~TransientIntegrator()
{

}

#if 1
int 
TransientIntegrator::formTangent(int statFlag, double iFact, double cFact)
{
  iFactor = iFact;
  cFactor = cFact;
  return this->formTangent(statFlag);
}
#endif

int 
TransientIntegrator::formTangent(int statFlag)
{
  statusFlag = statFlag;


  int result = 0;
  LinearSOE *theLinSOE = this->getLinearSOE();
  AnalysisModel *theModel = this->getAnalysisModel();
  if (theLinSOE == nullptr || theModel == nullptr) {
    opserr << "WARNING TransientIntegrator::formTangent() ";
    opserr << "no LinearSOE or AnalysisModel has been set\n";
    return -1;
  }

  // the loops to form and add the tangents are broken into two for 
  // efficiency when performing parallel computations
  
  theLinSOE->zeroA();

  // do modal damping
  bool inclModalMatrix = theModel->inclModalDampingMatrix();
  if (inclModalMatrix == true) {
    const Vector *modalValues = theModel->getModalDampingFactors();
    if (modalValues != 0) {
      this->addModalDampingMatrix(modalValues);
    }
  }


  // loop through the DOF_Groups and add the unbalance
  DOF_GrpIter &theDOFs = theModel->getDOFs();
  DOF_Group *dofPtr; 
  while ((dofPtr = theDOFs()) != nullptr) {
    if (theLinSOE->addA(dofPtr->getTangent(this),dofPtr->getID()) <0) {
      opserr << "TransientIntegrator::formTangent() - failed to addA:dof\n";
      result = -1;
    }
  }

  // loop through the FE_Elements getting them to add the tangent    
  FE_EleIter &theEles2 = theModel->getFEs();    
  FE_Element *elePtr;    
  while((elePtr = theEles2()) != nullptr) {
    if (theLinSOE->addA(elePtr->getTangent(this),elePtr->getID()) < 0) {
      opserr << "TransientIntegrator::formTangent() - failed to addA:ele\n";
      result = -2;
    }
  }
  return result;
}



int
TransientIntegrator::formUnbalance()
{
  // The dynamic residual is the sum of three sources:
  //   Pm: Modal damping force
  //   Pe: Element residual force
  //   Pn: Nodal unbalance force, including Rayleigh forces
  //
  LinearSOE *theLinSOE = this->getLinearSOE();
  AnalysisModel *theModel = this->getAnalysisModel();

  if (theModel == nullptr || theLinSOE == nullptr) {
    opserr << "WARNING IncrementalIntegrator::formUnbalance -";
    opserr << " no AnalysisModel or LinearSOE has been set\n";
    return -1;
  }

  theLinSOE->zeroB();

  // do modal damping
  const Vector *modalValues = theModel->getModalDampingFactors();
  if (modalValues != nullptr)
    this->addModalDampingForce(modalValues);

  if (this->formElementResidual() < 0)
    return -1;

  if (this->formNodalUnbalance() < 0)
    return -2;

  return 0;
}

int
TransientIntegrator::formEleResidual(FE_Element *theEle)
{
  theEle->zeroResidual();
  theEle->addRIncInertiaToResidual(1.0);
  return 0;
}    

int
TransientIntegrator::formNodUnbalance(DOF_Group *theDof)
{
  theDof->zeroUnbalance();
  theDof->addPIncInertiaToUnbalance();
  return 0;
}    


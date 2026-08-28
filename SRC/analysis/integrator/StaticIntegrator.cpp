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
// Description: This file contains the class definition for StaticIntegrator.
// StaticIntegrator is an algorithmic class for setting up the finite element
// equations for a static analysis and for Incrementing the nodal displacements
// with the values in the soln vector to the LinearSOE object. 
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#include <StaticIntegrator.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>

#include <DOF_Group.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <DOF_Group.h>


StaticIntegrator::StaticIntegrator(int clasTag)
 : IncrementalIntegrator(clasTag)
{
  // perhaps this should go in the constructor for StaticIntegrator;
  // The children should declare it
  this->setResidualType(ResidualType::StaticUnbalance);
}

StaticIntegrator::~StaticIntegrator()
{

}

int 
StaticIntegrator::formUnbalance()
{
  LinearSOE* theLinSOE = this->getLinearSOE();

  if (theLinSOE == nullptr)
    return -1;

  theLinSOE->zeroB();

  if (this->getAnalysisModel()->applyResidual(*this, *theLinSOE) < 0)
    return -1;

  return 0;
}



int
StaticIntegrator::formEleTangent(FE_Element *theEle)
{
  if (statusFlag == CURRENT_TANGENT) {
    theEle->zeroTangent();
    theEle->addKtToTang();

  } else if (statusFlag == INITIAL_TANGENT) {
    theEle->zeroTangent();
    theEle->addKiToTang();

  } else if (statusFlag == HALL_TANGENT) {
    theEle->zeroTangent();
    theEle->addKtToTang(cFactor);
    theEle->addKiToTang(iFactor);
  } 

  return 0;
}

int
StaticIntegrator::formEleResidual(FE_Element *theEle)
{
  theEle->zeroResidual();
  switch (this->getResidualType()) {
    case ResidualType::StaticUnbalance:
      theEle->addRtoResidual();
      break;
    case ResidualType::StaticSensitivity:
      theEle->addResistingForceSensitivity(this->getGradIndex());
      break;
    case ResidualType::TransientUnbalance:
      theEle->addRIncInertiaToResidual(1.0);
      break;
  }
  return 0;
}

int
StaticIntegrator::formNodTangent(DOF_Group *theDof)
{
  // For most static analysis, nodes dont contribute to tangent.
  theDof->zeroTangent();
//theDof->addKtoTang(1.0);
  return 0;
}    

int
StaticIntegrator::formNodUnbalance(DOF_Group *theDof)
{
  // only nodes unbalance need be added
  theDof->zeroUnbalance();
  theDof->addPtoUnbalance();
  return 0;
}    

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
// Description: This file contains the class definition for LoadControl.
// LoadControl is an algorithmic class for performing a static analysis
// using a load control integration scheme.
// 
// Written: fmk 
// Created: 07/98
//
#include <cmath>
#include <LoadControl.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Node.h> // getDOF_GroupPtr
#include <DOF_Group.h>
#include <LoadPattern.h>
#include <Domain.h>
#include <Parameter.h>
#include <ParameterIter.h>


LoadControl::LoadControl(double dLambda, int numIncr, double min, double max, double exponent)
  : StaticIntegrator(INTEGRATOR_TAGS_LoadControl),
    deltaLambda(dLambda), 
    specNumIncrStep(numIncr), 
    numIncrLastStep(numIncr),
    dLambdaMin(min), 
    dLambdaMax(max),
    pcontrol(dLambda, numIncr, min, max, exponent)
{
  // to avoid divide-by-zero error on first update() ensure numIncr != 0
  if (numIncr == 0) {
    opserr << "WARNING LoadControl::LoadControl() - numIncr set to 0, 1 assumed\n";
    specNumIncrStep = 1.0;
    numIncrLastStep = 1.0;
  }
}


LoadControl::~LoadControl()
{
    
}

int
LoadControl::domainChanged()
{
  AnalysisModel *theModel = this->getAnalysisModel();    
  if (theModel == nullptr) {
    return -1;
  }

  pcontrol.domainChanged(*this->getLinearSOE(), *theModel);

  return 0;
}

// increment
int 
LoadControl::newStep()
{

  AnalysisModel *theModel = this->getAnalysisModel();    
  if (theModel == nullptr) {
    return -1;
  }
  double current_time = theModel->getCurrentDomainTime();

  double deltaLambda;
  if (pcontrol.predict(*this->getLinearSOE(), deltaLambda) < 0) {
    opserr << "LoadControl::newStep() - predictor control failed\n";
    return -1;
  }


  current_time += deltaLambda;
  theModel->applyLoadDomain(current_time);

  return 0;
}


int
LoadControl::update(const Vector &deltaU)
{
  AnalysisModel *myModel = this->getAnalysisModel();
  LinearSOE *theSOE = this->getLinearSOE();
  if (myModel == nullptr || theSOE == nullptr)
    return -10;

  myModel->incrDisp(deltaU);
  if (myModel->updateDomain() < 0) {
    opserr << "LoadControl::update - model failed to update for new dU\n";
    return -1;
  }

  // Set deltaU for the convergence test
  theSOE->setX(deltaU);

  numIncrLastStep++;
  if (pcontrol.update(deltaU) < 0) {
    opserr << "LoadControl::update() - predictor control failed\n";
    return -1;
  }

  return 0;
}


int
LoadControl::commit()
{
  int result = StaticIntegrator::commit();
  if (result == 0)
    pcontrol.commit();
  return result;
}


int
LoadControl::revertToLastStep()
{
  pcontrol.revert();
  return 0;
}




int
LoadControl::formSensitivityRHS(int grad)
{
  // Set a couple of data members
  this->setResidualType(ResidualType::StaticSensitivity);
  this->setGradIndex(grad);


  // get model
  AnalysisModel* theAnalysisModel = this->getAnalysisModel();
  LinearSOE* theSOE = this->getLinearSOE();

  //
  // add residual
  //
  FE_Element *elePtr;
  FE_EleIter &theEles = theAnalysisModel->getFEs();
  while((elePtr = theEles()) != nullptr)
    theSOE->addB( elePtr->getResidual(this),  elePtr->getID());


  //
  // add dPext/dh contributions
  //
  static Vector oneDimVectorWithOne(1);
  oneDimVectorWithOne(0) = 1.0;
  static ID oneDimID(1);

  LoadPattern *loadPatternPtr;
  Domain *theDomain = theAnalysisModel->getDomainPtr();
  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

  while ((loadPatternPtr = thePatterns()) != nullptr) {
    const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(grad);
    int sizeRandomLoads = randomLoads.Size();
    if (sizeRandomLoads == 1) {
      ;  // No random loads in this load pattern
    }
    else {
      int numRandomLoads = (int)(sizeRandomLoads/2);
      for (int i=0; i<numRandomLoads*2; i=i+2) {
        int nodeNumber = (int)randomLoads(i);
        int dofNumber = (int)randomLoads(i+1);
        Node* aNode = theDomain->getNode(nodeNumber);
        DOF_Group* aDofGroup = aNode->getDOF_GroupPtr();
        const ID &anID = aDofGroup->getID();
        int relevantID = anID(dofNumber-1);
        oneDimID(0) = relevantID;
        theSOE->addB(oneDimVectorWithOne, oneDimID);
      }
    }
  }

  // reset residual type
  this->setResidualType(ResidualType::StaticUnbalance);

  return 0;
}

int
LoadControl::updateGradient(const Vector &v, int gradNum, int numGrads)
{
  // get model
  AnalysisModel* theAnalysisModel = this->getAnalysisModel();
  
  DOF_GrpIter &theDOFGrps = theAnalysisModel->getDOFs();

  DOF_Group   *dofPtr;
  
  while ( (dofPtr = theDOFGrps() ) != nullptr)
    dofPtr->saveDispSensitivity(v,gradNum,numGrads); 
  
  return 0;
}

int 
LoadControl::commitSensitivity(int gradNum, int numGrads)
{
  return 0;
}

// false for LC and true for DC
bool 
LoadControl::computeSensitivityAtEachIteration()
{
  return false;
}


int 
LoadControl::computeSensitivities()
{
  LinearSOE *theSOE = this->getLinearSOE();
  AnalysisModel *theModel = this->getAnalysisModel();

  // Zero out the old right-hand side of the SOE
  theSOE->zeroB();

  Domain *theDomain=theModel->getDomainPtr();
  ParameterIter &paramIter = theDomain->getParameters();      
  Parameter *theParam;
  // De-activate all parameters
  while ((theParam = paramIter()) != nullptr)
    theParam->activate(false);

  // compute sensitivity wrt each parameter
  int numGrads = theDomain->getNumParameters();
  paramIter = theDomain->getParameters();      
  while ((theParam = paramIter()) != nullptr) {

    // Activate this parameter
    theParam->activate(true);

    // Zero the RHS vector
    theSOE->zeroB();

    // Get the grad index for this parameter
    int gradIndex = theParam->getGradIndex();

    // Form the RHS
    this->formSensitivityRHS(gradIndex);
      
    // Solve for displacement sensitivity
    
    theSOE->solve();

    // Save sensitivity to nodes
    this->updateGradient(theSOE->getX(), gradIndex, numGrads);
    
    // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
    theModel->commitGradient(gradIndex, numGrads);
    this->commitSensitivity(gradIndex, numGrads);
    
    // De-activate this parameter for next sensitivity calc
    theParam->activate(false);
  }

  return 0;
}


void
LoadControl::Print(OPS_Stream &s, int flag)
{ 

  s << "LoadControl Integrator\n";
  s << "  expon: " << expon << "\n";
  s << "  specNumIncrStep: " << specNumIncrStep << "\n";
  s << "  numIncrLastStep: " << numIncrLastStep << "\n";
  s << "  dLambdaMin: " << dLambdaMin << "\n";
  s << "  dLambdaMax: " << dLambdaMax << "\n";

  AnalysisModel *theModel = this->getAnalysisModel();

  if (theModel != nullptr) {
    double currentLambda = theModel->getCurrentDomainTime();
    s << "\t LoadControl - currentLambda: " << currentLambda;
    s << "  deltaLambda: " << deltaLambda << "\n";
  }

}

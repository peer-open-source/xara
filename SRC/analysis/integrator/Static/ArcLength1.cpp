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
// Description: This file contains the class definition for ArcLength1.
// ArcLength1 is an algorithmic class for performing a static analysis
// using the arc length scheme, that is within a load step the following
// constraint is enforced: dU^TdU + alpha^2*dLambda^2 = ArcLength1^2
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and ArcLength1 is a control parameter.
//
// File: ~/analysis/integrator/ArcLength1.C
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
#include <ArcLength1.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <stdlib.h>
#include <elementAPI.h>

void *
OPS_ADD_RUNTIME_VPV(OPS_ArcLength1)
{
    double arcLength;
    double alpha;
    if (OPS_GetNumRemainingInputArgs() < 2) {
      opserr << "WARNING integrator ArcLength arcLength alpha \n";
      return 0;
    }

    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &arcLength) < 0) {
      opserr << "WARNING integrator ArcLength failed to read arc length\n";
      return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &alpha) < 0) {
      opserr << "WARNING integrator ArcLength failed to read alpha\n";
      return 0;
    }
    return new ArcLength1(arcLength,alpha); 
}

ArcLength1::ArcLength1(double arcLength, double alpha)
:StaticIntegrator(INTEGRATOR_TAGS_ArcLength1),
 arcLength2(arcLength*arcLength), alpha2(alpha*alpha),
 deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), 
 phat(0), deltaLambdaStep(0.0), currentLambda(0.0), 
 signLastDeltaLambdaStep(1)
{

}

ArcLength1::~ArcLength1()
{
  // delete any vector object created
  if (deltaUhat != nullptr)
    delete deltaUhat;
  if (deltaU != nullptr)
    delete deltaU;
  if (deltaUstep != nullptr)
    delete deltaUstep;
  if (deltaUbar != nullptr)
    delete deltaUbar;
  if (phat != nullptr)
    delete phat;
}


int
ArcLength1::newStep()
{
  // get pointers to AnalysisModel and LinearSOE
  AnalysisModel *theModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();
  assert(theModel != nullptr && theLinSOE != nullptr);

  // get the current load factor
  currentLambda = theModel->getCurrentDomainTime();

  if (deltaLambdaStep < 0)
    signLastDeltaLambdaStep = -1;
  else
    signLastDeltaLambdaStep = +1;

  // determine dUhat
  this->formTangent();
  theLinSOE->setB(*phat);
  theLinSOE->solve();
  (*deltaUhat) = theLinSOE->getX();
  Vector &dUhat = *deltaUhat;
  
  // determine delta lambda(1) == dlambda
  double dLambda = sqrt(arcLength2/((dUhat^dUhat) + alpha2));
  dLambda *= signLastDeltaLambdaStep; // base sign of load change
                                      // on what was happening last step
  deltaLambdaStep = dLambda;
  currentLambda += dLambda;

  // determine delta U(1) == dU
  (*deltaU) = dUhat;
  (*deltaU) *= dLambda;
  (*deltaUstep) = (*deltaU);

  // update model with delta lambda and delta U
  theModel->incrDisp(*deltaU);    
  theModel->applyLoadDomain(currentLambda);    
  theModel->updateDomain();

  return 0;
}

int
ArcLength1::update(const Vector &dU)
{
  AnalysisModel *theModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();    
  assert(theModel != nullptr && theLinSOE != nullptr);

  (*deltaUbar) = dU;

  // determine dUhat    
  theLinSOE->setB(*phat);
  theLinSOE->solve();
  (*deltaUhat) = theLinSOE->getX();    

  // determine delta lambda(i)
  double a =  (*deltaUstep)^(*deltaUbar);
  double b = ((*deltaUstep)^(*deltaUhat)) + alpha2*deltaLambdaStep;
  if (b == 0) {
    opserr << "ArcLength1::update() - zero denominator,";
    opserr << " alpha was set to 0.0 and zero reference load\n";
    return -1;
  }
  double dLambda = -a/b;

  // determine delta U(i)
  (*deltaU) = (*deltaUbar);    
  deltaU->addVector(1.0, *deltaUhat,dLambda);
  
  // update dU and dlambda
  (*deltaUstep) += *deltaU;
  deltaLambdaStep += dLambda;
  currentLambda += dLambda;

  // update the model
  theModel->incrDisp(*deltaU);    
  theModel->applyLoadDomain(currentLambda);
  theModel->updateDomain();
  
  // set the X soln in linearSOE to be deltaU for convergence Test
  theLinSOE->setX(*deltaU);

  return 0;
}



int 
ArcLength1::domainChanged()
{
  // we first create the Vectors needed
  AnalysisModel *theModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();    
  if (theModel == 0 || theLinSOE == 0) {
    opserr << "WARNING ArcLength1::update() ";
    opserr << "No AnalysisModel or LinearSOE has been set\n";
    return -1;
  }    
  int size = theModel->getNumEqn(); // ask model in case N+1 space

  if (deltaUhat == nullptr || deltaUhat->Size() != size) { // create new Vector
    if (deltaUhat != nullptr)
        delete deltaUhat;   // delete the old
    deltaUhat = new Vector(size);
  }

  if (deltaUbar == nullptr || deltaUbar->Size() != size) { // create new Vector
    if (deltaUbar != nullptr)
        delete deltaUbar;   // delete the old
    deltaUbar = new Vector(size);
  }

  
  if (deltaU == nullptr || deltaU->Size() != size) { // create new Vector
    if (deltaU != nullptr)
        delete deltaU;   // delete the old
    deltaU = new Vector(size);
  }

  if (deltaUstep == nullptr || deltaUstep->Size() != size) { 
    if (deltaUstep != nullptr)
        delete deltaUstep;  
    deltaUstep = new Vector(size);
  }

  if (phat == 0 || phat->Size() != size) { 
    if (phat != nullptr)
        delete phat;  
    phat = new Vector(size);
  }    

  // now we have to determine phat
  // do this by incrementing lambda by 1, applying load
  // and getting phat from unbalance.
  currentLambda = theModel->getCurrentDomainTime();
  currentLambda += 1.0;
  theModel->applyLoadDomain(currentLambda);    
  this->formUnbalance(); // NOTE: this assumes unbalance at last was 0
  (*phat) = theLinSOE->getB();
  currentLambda -= 1.0;
  theModel->setCurrentDomainTime(currentLambda);    
  
  return 0;
}


int
ArcLength1::sendSelf(int cTag, Channel &theChannel)
{
  Vector data(5);
  data(0) = arcLength2;
  data(1) = alpha2;
  data(2) = deltaLambdaStep;
  data(3) = currentLambda;
  data(4)  = signLastDeltaLambdaStep;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "ArcLength1::sendSelf() - failed to send the data\n";
      return -1;
  }
  return 0;
}


int
ArcLength1::recvSelf(int cTag,
                Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(5);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "ArcLength1::sendSelf() - failed to send the data\n";
    return -1;
  }      

  // set the data
  arcLength2 = data(0);
  alpha2 = data(1);
  deltaLambdaStep = data(2);
  currentLambda = data(3);
  signLastDeltaLambdaStep = data(4);
  return 0;
}

void
ArcLength1::Print(OPS_Stream &s, int flag)
{
  AnalysisModel *theModel = this->getAnalysisModel();
  if (theModel != nullptr) {
    double cLambda = theModel->getCurrentDomainTime();
    s << "\t ArcLength1 - currentLambda: " << cLambda;
    s << "  ArcLength1: " << sqrt(arcLength2) <<  "  alpha: ";
    s << sqrt(alpha2) << endln;
  } else 
    s << "\t ArcLength1 - no associated AnalysisModel\n";
}



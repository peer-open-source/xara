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
// Description: This file contains the class definition for ArcLength.
// ArcLength is an algorithmic class for perfroming a static analysis
// using the arc length scheme, that is within a load step the follwing
// constraint is enforced: dU^TdU + alpha^2*dLambda^2 = arcLength^2
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and arcLength is a control parameter.
//
// Written: fmk
// Created: 07/98
// Revision: A
//
#include <math.h>
#include <assert.h>
#include <ArcLength.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Domain.h>
#include <ID.h>
#include <FE_Element.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <Node.h>
#include <DOF_Group.h>
#include <EquiSolnAlgo.h>
#include <control/TrackSign.h>

template <typename T> static int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

ArcLength::ArcLength(double arcLength, 
                     double alpha,
                     double numIter, 
                     double expon_, 
                     bool use_det_,
                     ReferencePattern reference_type_)
:StaticIntegrator(INTEGRATOR_TAGS_ArcLength),
 arcLength(arcLength),
 alpha2(alpha*alpha),
 expon(expon_),
 use_det(use_det_),
 reference_type(reference_type_),
 numLastIter(numIter),
 numSpecIter(numIter),
 dUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0),
 phat(0), 
 deltaLambdaStep(0.0),
 dDeltaLambdaStepdh(0.0), 
#ifndef NO_ARC_LAMBDA
 currentLambda(0.0), 
#endif
 dLAMBDA(0.0),
 a(0.0),b(0.0),c(0.0),b24ac(0.0),
 //
 deltaUstep2(0),dDeltaUstepdh(0), 
 dLAMBDA2(0.0),dlambda1dh(0.0),dLAMBDAdh(0),Residual(0),sensU(0),sensitivityFlag(0),
 signLastDeltaLambdaStep(1), signLastDeterminant(1), 
 dUhatdh(0),dphatdh(0),dUIJdh(0),dlambdaJdh(0.0),gradNumber(0)
{

}

ArcLength::~ArcLength()
{
  if(dUhatdh !=0)
    delete dUhatdh;
  if(dphatdh !=0)
    delete dphatdh;
  if(dLAMBDAdh !=0)
    delete dLAMBDAdh;
  if (dUIJdh !=0)
    delete dUIJdh;
  if (dDeltaUstepdh !=0)
    delete dDeltaUstepdh;
  if (Residual !=0)
    delete Residual;
  if (sensU !=0)
    delete sensU;
}

int
ArcLength::newStep()
{
  // get pointers to AnalysisModel and LinearSOE
  AnalysisModel *theModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();
  assert(theModel != nullptr && theLinSOE != nullptr);



  //
  // Determine dUhat
  //
  this->formTangent();
  if (theLinSOE->solve(phat, dUhat) < 0)
    return -1;



  //
  // 3) determine delta lambda(1) == dlambda
  //
  // Determine sign of increment
  double signStep = 1;
  if (deltaLambdaStep < 0)
      signStep = -1;

  if (use_det) {
    double det = theLinSOE->getDeterminant();
    int sgndet = sgn(det); //(isnan(det) || det == 0.0)? -1 : sgn(det);

    if (sgndet != signLastDeterminant && !isnan(det))
      signStep *= -1;

    signLastDeterminant = sgndet;
  }
  // Determine the size of the increment
  arcLength     *= std::pow(numSpecIter/numLastIter, expon);
  double dLambda = std::sqrt((arcLength*arcLength)/((dUhat^dUhat) + alpha2));
  dLambda *= signStep; // base sign of load change on what was happening last step

  //
  // 4) Increment current load factor
  //
#ifdef NO_ARC_LAMBDA
  double 
#endif
  currentLambda = theModel->getCurrentDomainTime();
  deltaLambdaStep = dLambda;
  currentLambda  += dLambda;
  numLastIter = 0;

  //
  // determine delta U(1) == dU = dUhat * dlambda
  //
  deltaU     = dUhat;
  deltaU    *= dLambda;
  deltaUstep  = deltaU;

  // update model with delta lambda and delta U
  theModel->incrDisp(deltaU);    

  if (this->activateSensitivity() == true) {
    (*deltaUstep2) = deltaU;
    dLAMBDA         = dLambda;
    Domain *theDomain=theModel->getDomainPtr();
    ParameterIter &paramIter = theDomain->getParameters();
    Parameter *theParam;

    // De-activate all parameters
    while ((theParam = paramIter()) != nullptr)
      theParam->activate(false);

    // Now, compute sensitivity wrt each parameter
    // int numGrads = theDomain->getNumParameters();
    paramIter = theDomain->getParameters();
    while ((theParam = paramIter()) != nullptr) {
      // Activate this parameter
      theParam->activate(true);
      // Get the grad index for this parameter
      gradNumber = theParam->getGradIndex();

      this->formTangDispSensitivity(gradNumber);

      this->formdLambdaDh(gradNumber);

      dDeltaUstepdh->addVector(0.0, *dUhatdh,   dLambda);
      dDeltaUstepdh->addVector(1.0, dUhat, dlambda1dh);
      dDeltaLambdaStepdh =dlambda1dh; 
      theParam->activate(false);
    }
  }


  theModel->applyLoadDomain(currentLambda);    
  if (theModel->updateDomain() < 0)
    return -1;

  return 0;
}


int
ArcLength::update(const Vector &dX)
{

  AnalysisModel *theModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();    
  assert(theModel != nullptr && theLinSOE != nullptr);

  // have to do this as the SOE is going to change dX
  // when we call solve?
  deltaUbar = dX;

  // determine dUhat
  theLinSOE->solve(phat, dUhat);   

  // determine the coeeficients of our quadratic equation
  a = alpha2 + (dUhat^dUhat);
  b = alpha2*deltaLambdaStep 
    + (dUhat^deltaUbar)
    + (deltaUstep^dUhat);
  b *= 2.0;
  c = 2.0*(deltaUstep^deltaUbar) + (deltaUbar^deltaUbar);

  // check for a solution to quadratic
  b24ac = b*b - 4.0*a*c;

  if (b24ac < 0) {
    opserr << "ArcLength - imaginary roots due to multiple instability";
    opserr << " directions - initial load increment was too large\n";
    opserr << "a: " << a << " b: " << b << " c: " << c << " b24ac: " << b24ac << "\n";
    return -1;
  }

  double a2 = 2.0*a;
  if (a2 == 0.0) {
    opserr << "ArcLength - zero denominator";
    opserr << " alpha was set to 0.0 and zero reference load\n";
    return -2;
  }

  // determine the roots of the quadratic
  double dLambda;
  {
    double sqrtb24ac = sqrt(b24ac);
    double dlambda1 = (-b + sqrtb24ac)/a2;
    double dlambda2 = (-b - sqrtb24ac)/a2;

    double val    =  dUhat^deltaUstep;

    double theta1 = (deltaUstep^deltaUstep)
                  + (deltaUbar^deltaUstep)
                  + dlambda1*val;
#if 0
    // choose dLambda based on angle between incremental displacement before
    // and after this step -- want positive
    if (theta1 > 0)
      dLambda = dlambda1;
    else
      dLambda = dlambda2;
#else 
    double theta2 =
        (deltaUstep ^ deltaUstep)
      + (deltaUbar  ^ deltaUstep)
      + dlambda2 * val;

    if (true) {
      theta1 += alpha2*deltaLambdaStep*(deltaLambdaStep + dlambda1);
      theta2 += alpha2*deltaLambdaStep*(deltaLambdaStep + dlambda2);
    }
    if (theta1 > 0.0 && theta2 > 0.0) {
      const double dLambdaLinear = (b != 0.0) ? -c / b : 0.0;

      if (std::abs(dlambda1 - dLambdaLinear)
        <= std::abs(dlambda2 - dLambdaLinear))
        dLambda = dlambda1;
      else
        dLambda = dlambda2;
    }
    else if (theta1 > 0.0) {
      dLambda = dlambda1;
    }
    else if (theta2 > 0.0) {
      dLambda = dlambda2;
    }
    else {
      // No root satisfies the C&H angle test.
      return -1;
    }
#endif
  }
  dLAMBDA2 = dLambda;

  // determine delta U(i)
  deltaU = deltaUbar;    
  deltaU.addVector(1.0, dUhat,dLambda);
  
  // update dU and dlambda
  deltaUstep      += deltaU;
  deltaLambdaStep += dLambda;
#ifdef NO_ARC_LAMBDA
  double currentLambda = theModel->getCurrentDomainTime();
#endif
  currentLambda   += dLambda;

  // update the model
  theModel->incrDisp(deltaU);    
  theModel->applyLoadDomain(currentLambda);

  theModel->updateDomain();

  // set the X soln in linearSOE to be deltaU for convergence Test
  theLinSOE->setX(deltaU);

  numLastIter++;
  return 0;
}



int 
ArcLength::domainChanged()
{
  // we first create the Vectors needed
  AnalysisModel *theModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();    
  assert(theModel != nullptr && theLinSOE != nullptr);

  int size = theModel->getNumEqn(); // ask model in case N+1 space

  dUhat.resize(size);
  deltaUbar.resize(size);
  deltaU.resize(size);
  deltaUstep.resize(size);
  phat.resize(size);

  if (deltaUstep2 == nullptr || deltaUstep2->Size() != size) { 
    if (deltaUstep2 != nullptr)
      delete deltaUstep2;  
    deltaUstep2 = new Vector(size);
  }
  if (dDeltaUstepdh == nullptr || dDeltaUstepdh->Size() != size) { 
    if (dDeltaUstepdh != nullptr)
      delete dDeltaUstepdh;  
    dDeltaUstepdh = new Vector(size);
  }
  if (dphatdh == nullptr || dphatdh->Size() != size) { 
    if (dphatdh != nullptr)
      delete dphatdh;  
    dphatdh = new Vector(size);
  }
  if (dUhatdh == nullptr || dUhatdh->Size() != size) { 
    if (dUhatdh != nullptr)
      delete dUhatdh;  
    dUhatdh = new Vector(size);
  } 
  if (dUIJdh == nullptr || dUIJdh->Size() != size) { 
    if (dUIJdh != nullptr)
      delete dUIJdh;  
    dUIJdh = new Vector(size);
  }
  if (Residual == nullptr || Residual->Size() != size) { 
    if (Residual != nullptr)
      delete Residual;  
    Residual = new Vector(size);
  }
  if (sensU == nullptr || sensU->Size() != size) { 
    if (sensU != nullptr)
      delete sensU;  
    sensU = new Vector(size);
  }

  Domain *theDomain = theModel->getDomainPtr();
  int numGrads = theDomain->getNumParameters();

  if (dLAMBDAdh == 0 || dLAMBDAdh->Size() != (numGrads)) { 
    if (dLAMBDAdh != nullptr)
      delete dLAMBDAdh;
    dLAMBDAdh = new Vector(numGrads);
  }

  //
  // now we have to determine phat
  // do this by incrementing lambda by 1, applying load
  // and getting phat from unbalance.
  //
  // need to save currentLambda because calling
  // applyLoadDomain will change it (in the domain)
  if (false) {
    // Original
#ifdef NO_ARC_LAMBDA
    double 
#endif
    currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);
    this->formUnbalance(phat); // NOTE: this assumes unbalance at last was 0
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);
  } 
  else {
#ifdef NO_ARC_LAMBDA
    double 
#endif
    currentLambda = theModel->getCurrentDomainTime();
    theModel->applyLoadDomain(1.0);

    switch (reference_type) {
    case ReferencePattern::Point:
      theLinSOE->zeroB();
      this->formNodalUnbalance(phat);
      break;
    case ReferencePattern::Body:
      theLinSOE->zeroB();
      this->formElementResidual(phat);
      break;
    case ReferencePattern::Full:
      this->formUnbalance(phat);
      break;
    }

    theModel->setCurrentDomainTime(currentLambda);    
  }


  // check there is a reference load
  bool haveLoad = false;
  for (int i=0; i<size; i++)
    if ( phat(i) != 0.0 ) {
      haveLoad = true;
      break;
    }

  if (haveLoad == false) {
    opserr << "WARNING ArcLength::domainChanged() - zero reference load";
    return -1;
  }

  return 0;  
}




void
ArcLength::Print(OPS_Stream &s, int flag)
{
  s << "\t ArcLength\n";
}



// Added by Abbas
// obtain the derivative of the tangent displacement (dUhatdh)
//
//
void
ArcLength::formTangDispSensitivity(int gradNumber)
{
   AnalysisModel *theModel=this->getAnalysisModel();
   LinearSOE *theLinSOE = this->getLinearSOE(); 

   // To get the structural stiffness Matrix

   dphatdh->Zero();
   this->formTangent();

   if (theLinSOE->solve(*dphatdh, *dUhatdh)<0)
     opserr << "SOE failed to obtained dUhatdh ";

   // if the parameter is a load parameter.
   // add the dPext/dh contributions

   static Vector oneDimVectorWithOne(1);
   oneDimVectorWithOne(0) = 1.0;
   static ID oneDimID(1);

   LoadPattern *loadPatternPtr;
   Domain *theDomain = theModel->getDomainPtr();
   LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

   while ((loadPatternPtr = thePatterns()) != 0) {
      const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(gradNumber);
      int sizeRandomLoads = randomLoads.Size();
      if (sizeRandomLoads == 1) {
	   // No random loads in this load pattern
      }
      else {
         // Random loads: add contributions to the 'B' vector
         int numRandomLoads = (int)(sizeRandomLoads/2);
         for (int i=0; i<numRandomLoads*2; i=i+2) {
            int nodeNumber = (int)randomLoads(i);
            int dofNumber = (int)randomLoads(i+1);
            Node *aNode = theDomain->getNode(nodeNumber);
            DOF_Group *aDofGroup = aNode->getDOF_GroupPtr();
            const ID &anID = aDofGroup->getID();
            int relevantID = anID(dofNumber-1);
            oneDimID(0) = relevantID;
            theLinSOE->addB(oneDimVectorWithOne, oneDimID);
            (*dphatdh)=theLinSOE->getB();
         }
      }
   }

   if (theLinSOE->solve()<0) {
      opserr << "SOE failed to obtained dUhatdh ";
   }
}


// form dLambda for each time step dLambda
double 
ArcLength::formdLambdaDh(int gradNumber)
{
  
  double dUhatTdUhat  = (dUhat^dUhat);
  double dUhatTUhatdh = dUhat^(*dUhatdh);

  if (dLAMBDA==0.0 ) {
    dlambda1dh=0.0;

  } else {
  // double ALPHA2=alpha2*alpha2;
    double denomerator=pow((dUhatTdUhat+alpha2),2.0);
    dlambda1dh = signLastDeltaLambdaStep *1.0/(dLAMBDA)*(-arcLength*dUhatTUhatdh/(denomerator));

   //dlambda1dh=1.0/(2.0)*dLAMBDA*(-sqrt(arcLength)*dUhatTUhatdh/(dUhatTdUhat+alpha2));

   //dlambda1dh *=signLastDeltaLambdaStep;

  }

  if (dLAMBDAdh != 0) {
    (*dLAMBDAdh)(gradNumber) = (*dLAMBDAdh)(gradNumber) + dlambda1dh;
    return (*dLAMBDAdh)(gradNumber);

  } else {
    return 0.0;
  }
}


// dLambdadh of the subsequent iterations (J>1)
double 
ArcLength::getLambdaSensitivity(int gradNumber)
{

 // determine the coeeficients of our quadratic equation
  //  double a = alpha2 + (dUhat^dUhat);
 //   double b = alpha2*deltaLambdaStep 
  //    + (dUhat^deltaUbar)
 //     + (deltaUstep^dUhat);
 //   b *= 2.0;
 //   double c = 2*(deltaUstep^deltaUbar) + (deltaUbar^deltaUbar);
    // check for a solution to quadratic
  //  double b24ac = b*b - 4.0*a*c;
  if (b24ac < 0) {
    opserr << "ArcLength::update() - imaginary roots due to multiple instability";
    opserr << " directions - initial load increment was too large\n";
    opserr << "a: " << a << " b: " << b << " c: " << c << " b24ac: " << b24ac << endln;
    return -1;
  }

  double a2 = 2.0*a;
  if (a2 == 0.0) {
    return -2;
  }			       

  double dAdh = 2.0*(dUhat^(*dUhatdh));
  double dBdh = 2.0*(((*dUIJdh)^dUhat)+(deltaUbar^(*dUhatdh))+((*deltaUstep2)^(*dUhatdh))+((*dDeltaUstepdh)^dUhat)+alpha2*dDeltaLambdaStepdh);
  double dCdh = 2.0*(((*deltaUstep2)^(*dUIJdh))+((*dDeltaUstepdh)^deltaUbar)+(deltaUbar^(*dUIJdh)));//+2.0*(deltaUstep^(*dDeltaUstepdh));//+alpha2*dDeltaLambdaStepdh*((phat)^(phat)));


  double sqrtb24ac = sqrt(b24ac);
  double dSqrtb24acdh = (2.0*b*dBdh-4.0*(a*dCdh+dAdh*c))/(2.0*sqrtb24ac);
  double dlambda1     = (-b + sqrtb24ac)/a2;
  double dlambdaj1dh  = (a2*(-dBdh+dSqrtb24acdh)-((-b+sqrtb24ac)*2.0*dAdh))/(4.0*a*a);

  // double dlambda2 = (-b - sqrtb24ac)/a2;
  double dlambdaj2dh= (a2*(-dBdh-dSqrtb24acdh)-((-b-sqrtb24ac)*2.0*dAdh))/(4.0*a*a);
 
  double val = dUhat^(*deltaUstep2);
  // double theta1 = ((*deltaUstep2)^(*deltaUstep2)) + (deltaUbar^(*deltaUstep2));
  // theta1 += dlambda1*val;

  double dTheta1dh = 2.0*((*deltaUstep2)^(*dDeltaUstepdh))+(deltaUbar^(*dDeltaUstepdh))+((*dUIJdh)^(*deltaUstep2));
  //    double theta2 = theta1 + dlambda2*val;
  double dvaldh= (dUhat^(*dDeltaUstepdh))+((*dUhatdh)^(*deltaUstep2));
  dTheta1dh +=dlambdaj1dh*val+dlambda1*dvaldh;
  // choose dLambda based on angle between incremental displacement before
  // and after this step -- want positive
  
  if (dTheta1dh > 0)
    dlambdaJdh = dlambdaj1dh;
  else
    dlambdaJdh = dlambdaj2dh;

  //  double dResultdh=a*2*dlambdaJdh*dLAMBDA2+dAdh*(dLAMBDA2*dLAMBDA2)+b*dlambdaJdh+dBdh*dLAMBDA2+dCdh;

  // determine delta U(i)

  //  sensU->addVector(1.0,*dUhatdh,dLAMBDA2);
  deltaU = deltaUbar;    
  deltaU.addVector(1.0, dUhat,dLAMBDA2);

  // update dU and dlambda
  (*deltaUstep2) += deltaU;

  dDeltaUstepdh->addVector(1.0, *dUhatdh,dLAMBDA2);
  dDeltaUstepdh->addVector(1.0, dUhat,dlambdaJdh);
  (*dDeltaUstepdh) +=(*dUIJdh);


  dDeltaLambdaStepdh +=dlambdaJdh;

 // dDeltaLambdaStepdh +=dlambdaJdh;
 // Now update Lambda_ij
  if (dLAMBDAdh !=0) {
     (*dLAMBDAdh)(gradNumber) = (*dLAMBDAdh)(gradNumber)+ dlambdaJdh;   
     return (*dLAMBDAdh)(gradNumber);
  } else { 
    return 0.0;
  }
}

void
ArcLength::formResidualDispSensitivity( int gradNumber)
{
  // AnalysisModel *theModel=this->getAnalysisModel();
  // int size=theModel->getNumEqn();
  // static Matrix dKdh(size,size);
  // dKdh=this->getdKdh(gradNumber);
  // dUIJdh-> addMatrixVector(1.0,dKdh,*deltaUbar,-1.0);  
}


int
ArcLength::formIndependentSensitivityRHS()
{
  return 0;
}


int
ArcLength::formSensitivityRHS(int passedGradNumber)
{
  // Set a couple of data members
  sensitivityFlag = 1;

  gradNumber = passedGradNumber;

  // get model
  AnalysisModel* theAnalysisModel = this->getAnalysisModel();
  LinearSOE* theSOE = this->getLinearSOE();

  // Loop through elements
  FE_Element *elePtr;
  FE_EleIter &theEles = theAnalysisModel->getFEs(); 
  while((elePtr = theEles()) != nullptr) {
    theSOE->addB(elePtr->getResidual(this) ,elePtr->getID()  );
  }


  (*Residual) = theSOE->getB();

  double CallDlambda1dh=(*dLAMBDAdh)(gradNumber);
  Residual->addVector(1.0, phat, CallDlambda1dh); //needed to calculate dLambdadh
 // Residual->addVector(1.0,*dphatdh, currentLambda ); //needed to calculate dLambdadh

/////////////
//int size=theAnalysisModel->getNumEqn();
//static Matrix dKdh(size,size);
//dKdh=this->getdKdh(gradNumber);
//Residual-> addMatrixVector(1.0,dKdh,*deltaUbar,-1.0);

  theSOE->setB(*Residual);


  // Loop through the loadPatterns and add the dPext/dh contributions
  static Vector oneDimVectorWithOne(1);
  oneDimVectorWithOne(0) = 1.0;
  static ID oneDimID(1);

  Node *aNode;
  DOF_Group *aDofGroup;
  int nodeNumber, dofNumber, relevantID, i, sizeRandomLoads, numRandomLoads;
  LoadPattern *loadPatternPtr;
  Domain *theDomain = theAnalysisModel->getDomainPtr();
  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
  while ((loadPatternPtr = thePatterns()) != nullptr) {
    const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(gradNumber);
    sizeRandomLoads = randomLoads.Size();
    if (sizeRandomLoads == 1) {
      // No random loads in this load pattern
    }
    else {
      // Random loads: add contributions to the 'B' vector
      numRandomLoads = (int)(sizeRandomLoads/2);
      for (i=0; i<numRandomLoads*2; i=i+2) {
        nodeNumber = (int)randomLoads(i);
        dofNumber = (int)randomLoads(i+1);
        aNode = theDomain->getNode(nodeNumber);
        aDofGroup = aNode->getDOF_GroupPtr();
        const ID &anID = aDofGroup->getID();
        relevantID = anID(dofNumber-1);
        oneDimID(0) = relevantID;
        theSOE->addB(oneDimVectorWithOne, oneDimID);
      }
    }
  }

  theSOE->setB(*Residual);

  //reset sensitivity flag
  sensitivityFlag = 0;

  return 0;
}


int
ArcLength::saveSensitivity(const Vector &v, int gradNum, int numGrads)
{
  AnalysisModel* theAnalysisModel = this->getAnalysisModel();

  DOF_GrpIter &theDOFGrps = theAnalysisModel->getDOFs();
  DOF_Group 	*dofPtr;

  while ( (dofPtr = theDOFGrps() ) != nullptr)
    dofPtr->saveDispSensitivity(v,gradNum,numGrads);

  return 0;
}


int
ArcLength::saveLambdaSensitivity(double dlambdadh, int gradNum, int numGrads)
{
  AnalysisModel* theAnalysisModel = this->getAnalysisModel();
  Domain *theDomain = theAnalysisModel->getDomainPtr();

  LoadPattern *lpPtr;
  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

  while ( (lpPtr = thePatterns() ) != nullptr)
    lpPtr->saveLoadFactorSensitivity(dlambdadh, gradNum, numGrads);

  return 0;
}


int
ArcLength::commitSensitivity(int gradNum, int numGrads)
{
  AnalysisModel* theAnalysisModel = this->getAnalysisModel();

  // Loop through the FE_Elements and set unconditional sensitivities
  FE_Element *elePtr;
  FE_EleIter &theEles = theAnalysisModel->getFEs();    
  while((elePtr = theEles()) != nullptr) {
    elePtr->commitSensitivity(gradNum, numGrads);
  }
  return 0;
}



bool 
ArcLength::computeSensitivityAtEachIteration()
{
  return true;
}


int 
ArcLength::computeSensitivities()
{
  AnalysisModel *theModel = this->getAnalysisModel();   
  Domain *theDomain=theModel->getDomainPtr();//Abbas
  LinearSOE *theSOE = this->getLinearSOE();

  // Zero out the old right-hand side of the SOE
  theSOE->zeroB();

  // Form the part of the RHS which are indepent of parameter
  this->formIndependentSensitivityRHS();

  // De-activate all parameters
  ParameterIter &paramIter = theDomain->getParameters();
  Parameter *theParam;
  while ((theParam = paramIter()) != nullptr)
     theParam->activate(false);

  // Compute sensitivity wrt each parameter
  int  numGrads = theDomain->getNumParameters();
  paramIter = theDomain->getParameters();
  while ((theParam = paramIter()) != nullptr) {

    // Activate this parameter
    theParam->activate(true);

    // Zero the RHS vector
    theSOE->zeroB();

    // Get the grad index for this parameter
    int gradIndex = theParam->getGradIndex();
  
    // Form the RHS
    this->formTangDispSensitivity(gradIndex);
    this->formSensitivityRHS(gradIndex);

    this->formTangent();
    theSOE->solve();
    *dUIJdh = theSOE->getX();// sensitivity of the residual displacement


  // this->formTangDispSensitivity(gradIndex);
    double dlamdh = this->getLambdaSensitivity(gradIndex);

    // To obtain the response sensitivity 
    theSOE->setB(*Residual);
    theSOE->solve();
    (*sensU) = theSOE->getX();
  //  this->formResidualDispSensitivity(gradNumber);
  // (*sensU)=(*dUIJdh);

  //  (*sensU) +=(*dUIJdh);

    //  sensU->addVector(1.0,*dUhatdh,dLAMBDA2);
    // sensU->addVector(1.0,*dUhat,dlamdh);
    //(*sensU) =(*dUIJdh);

    //(*dDeltaUstepdh) +=(*dUIJdh);

    //  (*sensU) +=(*dDeltaUstepdh);

    // Save sensitivity to nodes
    this->saveSensitivity( (*sensU), gradIndex, numGrads );
    
    
    this->saveLambdaSensitivity(dlamdh, gradIndex, numGrads);

    // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
    this->commitSensitivity(gradIndex, numGrads);

    // De-activate this parameter for next sensitivity calc
    theParam->activate(false);

  }
  // end of if statment to be run only one time during the iteration process.

  return 0;
}


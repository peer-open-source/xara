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
// Written: fmk 
// Created: 07/99
// Revision: A
//
#include <MinUnbalDispNorm.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <ID.h>
#include <stdlib.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <EquiSolnAlgo.h>
#include <TaggedObjectStorage.h>
#include <Matrix.h>


MinUnbalDispNorm::MinUnbalDispNorm(double lambda1, int specNumIter,
		     double min, double max, int signFirstStep)
:StaticIntegrator(INTEGRATOR_TAGS_MinUnbalDispNorm),
 dLambda1LastStep(lambda1), 
 specNumIncrStep(specNumIter), numIncrLastStep(specNumIter),
 deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0), dUhatdh(0), 
 dLambdaj(0.0),
 phat(nullptr), deltaLambdaStep(0.0), currentLambda(0.0), 
 dLambda1min(min), dLambda1max(max), signLastDeterminant(1), 
 signFirstStepMethod(signFirstStep),
 dLambdaStepDh(0.0),dUIJdh(0),Dlambdadh(0.0),dphatdh(0),Residual2(0),
 signLastDeltaLambdaStep(1), sensitivityFlag(0),Residual(0), 
 dlambdadh(0.0), dLambda(0.0), sensU(0),d_deltaU_dh(0),gradNumber(0),dLAMBDAdh(0)
{
  // to avoid divide-by-zero error on first update() ensure numIncr != 0
  if (specNumIncrStep == 0) {
    opserr << "WARNING LoadControl::LoadControl() - numIncr set to 0, 1 assumed\n";
    specNumIncrStep = 1.0;
    numIncrLastStep = 1.0;
  }
}

MinUnbalDispNorm::~MinUnbalDispNorm()
{
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

  if (dUhatdh !=0)
    delete dUhatdh;
  if (dUIJdh !=0)
    delete dUIJdh; 
  if (Residual !=0)
    delete Residual;
  if (sensU !=0)
    delete sensU;
  if (Residual2 !=0)
    delete Residual2;
  if (dLAMBDAdh !=0) 
    delete dLAMBDAdh;
  if (dphatdh !=0)
    delete dphatdh;

  dLAMBDAdh=0;
  dUhatdh=0;
}

int
MinUnbalDispNorm::newStep()
{
  // get pointers to AnalysisModel and LinearSOE
  AnalysisModel *theModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();    
  if (theModel == 0 || theLinSOE == 0) {
    opserr << "WARNING MinUnbalDispNorm::newStep() ";
    opserr << "No AnalysisModel or LinearSOE has been set\n";
    return -1;
  }

  // get the current load factor
  currentLambda = theModel->getCurrentDomainTime();

  // determine dUhat
  this->formTangent();
  theLinSOE->setB(*phat);
  if (theLinSOE->solve() < 0) {
    opserr << "MinUnbalanceDispNorm::newStep(void) - failed in solver\n";
    return -1;
  }
  (*deltaUhat) = theLinSOE->getX();
  Vector &dUhat = *deltaUhat;

  // determine delta lambda(1) == dlambda
  double expon   = 1.0;
  double dLambda = dLambda1LastStep*pow(specNumIncrStep/numIncrLastStep, expon);

  // check aaint min and max values specified in constructor
  if (dLambda < dLambda1min)
    dLambda = dLambda1min;

  else if (dLambda > dLambda1max)
    dLambda = dLambda1max;

  dLambda1LastStep = dLambda;

  if (signFirstStepMethod == SIGN_LAST_STEP) {
    if (deltaLambdaStep < 0)
      signLastDeltaLambdaStep = -1;
    else
      signLastDeltaLambdaStep = +1;
    
    dLambda *= signLastDeltaLambdaStep; // base sign of load change
                                        // on what was happening last step
  } else {

    double det = theLinSOE->getDeterminant();
    int signDeterminant = 1;
    if (det < 0)
      signDeterminant = -1;
  
    dLambda *= signDeterminant * signLastDeterminant;
  
    signLastDeterminant = signDeterminant;
  }

  /*
  double work = (*phat)^(dUhat);
  int signCurrentWork = 1;
  if (work < 0) signCurrentWork = -1;

  if (signCurrentWork != signLastDeltaStep)
  */

  deltaLambdaStep = dLambda;
  currentLambda  += dLambda;
  numIncrLastStep = 0;

  // determine delta U(1) == dU
  (*deltaU) = dUhat;
  (*deltaU) *= dLambda;
  (*deltaUstep) = (*deltaU);


  if (this->activateSensitivity() == true) { 
    Domain *theDomain=theModel->getDomainPtr();
    ParameterIter &paramIter = theDomain->getParameters();
    Parameter *theParam;

    // De-activate all parameters
    while ((theParam = paramIter()) != 0)
      theParam->activate(false);
     
    // Now, compute sensitivity wrt each parameter
    // int numGrads = theDomain->getNumParameters();
    
    
    paramIter = theDomain->getParameters();
    while ((theParam = paramIter()) != 0) {
      // Activate this parameter
      theParam->activate(true);
      // Get the grad index for this parameter
      gradNumber = theParam->getGradIndex();

      this->formTangDispSensitivity(dUhatdh,gradNumber);
      this->formdLambdaDh(gradNumber);

      sensU->addVector(1.0, *dUhatdh ,dLambda);

      theParam->activate(false);
    } 
  }
  ///////////////Abbas/////////////////////////////

  // update model with delta lambda and delta U
  theModel->incrDisp(*deltaU);    
  theModel->applyLoadDomain(currentLambda);
  if (theModel->updateDomain() < 0) {
    opserr << "MinUnbalDispNorm::newStep - model failed to update for new dU\n";
    return -1;
  }

  return 0;
}

int
MinUnbalDispNorm::update(const Vector &dU)
{
   AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();    
    if (theModel == 0 || theLinSOE == 0) {
      opserr << "WARNING MinUnbalDispNorm::update() ";
      opserr << "No AnalysisModel or LinearSOE has been set\n";
      return -1;
    }

    (*deltaUbar) = dU; // have to do this as the SOE is gonna change

    // determine dUhat    
    theLinSOE->setB(*phat);
    theLinSOE->solve();
    (*deltaUhat) = theLinSOE->getX();    

    // determine delta lambda(i)
    double a = (*deltaUhat)^(*deltaUbar);
    double b = (*deltaUhat)^(*deltaUhat);
    if (b == 0.0) {
      opserr << "MinUnbalDispNorm::update() - zero denominator\n";
      return -1;
    }

    double dLambda = -a/b;
    dLambdaj = dLambda; //Abbas
    // determine delta U(i)
    (*deltaU) = (*deltaUbar);    
    deltaU->addVector(1.0, *deltaUhat,dLambda);
    
    // update dU and dlambda
    (*deltaUstep)   += *deltaU;
    deltaLambdaStep += dLambda;
    currentLambda   += dLambda;
    // update the model
    theModel->incrDisp(*deltaU);    
    theModel->applyLoadDomain(currentLambda);    

    if (theModel->updateDomain() < 0) {
      opserr << "MinUnbalDispNorm::update - model failed to update for new dU\n";
      return -1;
    }
    
    // set the X soln in linearSOE to be deltaU for convergence Test
    theLinSOE->setX(*deltaU);

    numIncrLastStep++;

    return 0;
}



int 
MinUnbalDispNorm::domainChanged()
{
   // we first create the Vectors needed
   AnalysisModel *theModel = this->getAnalysisModel();
   LinearSOE *theLinSOE = this->getLinearSOE();    
   if (theModel == 0 || theLinSOE == 0) {
       opserr << "WARNING MinUnbalDispNorm::update() ";
       opserr << "No AnalysisModel or LinearSOE has been set\n";
       return -1;
   }    
   int size = theModel->getNumEqn(); // ask model in case N+1 space

   if (deltaUhat == 0 || deltaUhat->Size() != size) {
      if (deltaUhat != 0)
          delete deltaUhat;   // delete the old
      deltaUhat = new Vector(size);
   }

   if (deltaUbar == 0 || deltaUbar->Size() != size) {
      if (deltaUbar != 0)
          delete deltaUbar;   // delete the old
      deltaUbar = new Vector(size);
   }

   
   if (deltaU == 0 || deltaU->Size() != size) {
      if (deltaU != 0)
          delete deltaU;   // delete the old
      deltaU = new Vector(size);
   }

   if (deltaUstep == 0 || deltaUstep->Size() != size) { 
      if (deltaUstep != 0)
          delete deltaUstep;  
      deltaUstep = new Vector(size);
   }

   if (phat == nullptr || phat->Size() != size) { 
      if (phat != nullptr)
          delete phat;  
      phat = new Vector(size);
   }    

   if (dphatdh == 0 || dphatdh->Size() != size) { 
      if (dphatdh != 0)
        delete dphatdh;  
      dphatdh = new Vector(size);
   }


  if (dUhatdh == 0 || dUhatdh->Size() != size) { 
    if (dUhatdh != 0)
      delete dUhatdh;  
    dUhatdh = new Vector(size);
  } 
  
  if (dUIJdh == 0 || dUIJdh->Size() != size) { 
     if (dUIJdh != 0)
        delete dUIJdh;  
     dUIJdh = new Vector(size);
  }

  if (Residual == 0 || Residual->Size() != size) { 
     if (Residual != 0)
        delete Residual;  
     Residual = new Vector(size);
  } 
  
  if (Residual2 == 0 || Residual2->Size() != size) { 
    if (Residual2 != 0)
      delete Residual2;  
    Residual2 = new Vector(size);
  } 

  if (sensU == 0 || sensU->Size() != size) { 
     if (sensU != 0)
        delete sensU;  
     sensU = new Vector(size);
  } 

  Domain *theDomain=theModel->getDomainPtr();//Abbas
  int numGrads = theDomain->getNumParameters();

  if (dLAMBDAdh == 0 || dLAMBDAdh->Size() != (numGrads)) { 
    if (dLAMBDAdh != 0 )  
      delete dLAMBDAdh;
    dLAMBDAdh = new Vector(numGrads);
  } 
 
   // now we have to determine phat
   // do this by incrementing lambda by 1, applying load
   // and getting phat from unbalance.
  if (false) {
    currentLambda = theModel->getCurrentDomainTime();
    currentLambda += 1.0;
    theModel->applyLoadDomain(currentLambda);
    this->formUnbalance(); // NOTE: this assumes unbalance at last was 0
    (*phat) = theLinSOE->getB();// - *phat;
    currentLambda -= 1.0;
    theModel->setCurrentDomainTime(currentLambda);

  } else {
    // need to save currentLambda because calling
    // applyLoadDomain will change it (in the domain)
    currentLambda = theModel->getCurrentDomainTime();
    theModel->applyLoadDomain(1.0);
    theLinSOE->zeroB();
    this->formNodalUnbalance();
    (*phat) = theLinSOE->getB();
    theModel->setCurrentDomainTime(currentLambda);
  }

  // check there is a reference load
  bool hasLoad = false;
  for (int i=0; i<size; i++)
    if ( (*phat)(i) != 0.0 ) {
      hasLoad = true;
      break; // i = size;
    }

  if (hasLoad == false) {
    opserr << "WARNING zero reference load";
    return -1;
  }

  return 0;
}

int
MinUnbalDispNorm::sendSelf(int cTag, Channel &theChannel)
{
  Vector data(8);
  data(0) = dLambda1LastStep;
  data(1) = specNumIncrStep;
  data(2) = numIncrLastStep;
  data(3) = deltaLambdaStep;
  data(4) = currentLambda;
  if (signLastDeltaLambdaStep == 1)
    data(5)  = 1.0;
  else
    data(5) = 0.0;
  data(6) = dLambda1min;
  data(7) = dLambda1max;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "MinUnbalDispNorm::sendSelf() - failed to send the data\n";
      return -1;
  }
  return 0;
}


int
MinUnbalDispNorm::recvSelf(int cTag,
		    Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  Vector data(8);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
      opserr << "MinUnbalDispNorm::sendSelf() - failed to send the data\n";
      return -1;
  }      

  // set the data

  dLambda1LastStep = data(0);
  specNumIncrStep = data(1);
  numIncrLastStep = data(2);
  deltaLambdaStep = data(3);
  currentLambda = data(4);
  if (data(5)== 1.0)
    signLastDeltaLambdaStep = 1;
  else
    signLastDeltaLambdaStep = -1;
  dLambda1min = data(6);
  dLambda1max = data(7);

  return 0;
}

void
MinUnbalDispNorm::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != nullptr) {
	double cLambda = theModel->getCurrentDomainTime();
	s << "\t MinUnbalDispNorm - currentLambda: " << cLambda;
    } else 
	s << "\t MinUnbalDispNorm - no associated AnalysisModel\n";
}


///////////////////////////Sensitivity Begin///////////////////////

// obtain the derivative of the tangent displacement (dUhatdh)
Vector *
MinUnbalDispNorm::formTangDispSensitivity(Vector *dUhatdh,int gradNumber)
{
   LinearSOE *theLinSOE = this->getLinearSOE(); 
   dUhatdh->Zero();
   dphatdh->Zero();

   //call the tangent (K)
   this->formTangent();
   theLinSOE->setB(*dphatdh);
   if(theLinSOE->solve()<0) {
      opserr<<"SOE failed to obtained dUhatdh ";
      exit(-1); // TODO: why exit?
   }
   (*dUhatdh)=theLinSOE->getX();


   
   // if the parameter is a load parameter.
   ////////////////////////////////////////////////////////
   // Loop through the loadPatterns and add the dPext/dh contributions

   static Vector oneDimVectorWithOne(1);
   oneDimVectorWithOne(0) = 1.0;
   static ID oneDimID(1);
   Node *aNode;
   DOF_Group *aDofGroup;

   int nodeNumber, dofNumber, relevantID, i, sizeRandomLoads, numRandomLoads;
   
   LoadPattern *loadPatternPtr;
   AnalysisModel *theModel = this->getAnalysisModel();   
   Domain *theDomain = theModel->getDomainPtr();
   LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
  
   while ((loadPatternPtr = thePatterns()) != nullptr) {

     const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(gradNumber);
      sizeRandomLoads = randomLoads.Size();
      if (sizeRandomLoads == 1) {
	      // No random loads in this load pattern
        continue;
      }
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
         theLinSOE->addB(oneDimVectorWithOne, oneDimID);
         (*dphatdh)=theLinSOE->getB();
         // dphatdh->addMatrixVector(1.0,dKdh,*deltaUhat,1.0);
      }
   }

   if(theLinSOE->solve()<0) {
     opserr << "SOE failed to obtained dUhatdh ";
     exit(-1);
   }

   (*dUhatdh)=theLinSOE->getX();

   return dUhatdh;
}

// form dLambda for each time step dLambda
double 
MinUnbalDispNorm::formdLambdaDh(int gradNumber)
{
  // Here dLambda1dh=0 because its constant.
  if (dLAMBDAdh != 0)
    return (*dLAMBDAdh)(gradNumber);
  return 0.0;
}


 // dLambdadh of the subsequent iterations (J>1)
double 
MinUnbalDispNorm::getLambdaSensitivity(int gradNumber)
{

 
 //  Vector &dufRdh=*dUIJdh;// component of the dUfrDh: derivative of the residual displacement

   double temp= (*deltaUhat)^(*deltaUhat);
   double denomerator= pow(temp,2.0);
   double a= (*deltaUhat)^(*dUIJdh);
   double b= (*dUhatdh)^(*deltaUbar);
   double c= (*deltaUhat)^(*deltaUbar);
   double d= (*deltaUhat)^(*dUhatdh);
   double Numerator= -(temp*(a+b)-(c*2.0*d));
   Dlambdadh = Numerator/denomerator;  //   

   // Now update Lambda_ij
   if(dLAMBDAdh !=0) {
     (*dLAMBDAdh)(gradNumber) = (*dLAMBDAdh)(gradNumber)+ Dlambdadh;
     return (*dLAMBDAdh)(gradNumber);
   } else {
     return 0.0;
   }
}


int
MinUnbalDispNorm::formIndependentSensitivityRHS()
{
   return 0;
}

   
int
MinUnbalDispNorm::formSensitivityRHS(int passedGradNumber)
{

  sensitivityFlag = 1;
  //this->Activatesensitivity();//Abbas
  // Set a couple of data members
  gradNumber = passedGradNumber;

  // get model
  AnalysisModel* theAnalysisModel = this->getAnalysisModel();
  LinearSOE* theSOE = this->getLinearSOE();
  //theSOE->zeroB();

  // Loop through elements
  FE_Element *elePtr;
  FE_EleIter &theEles = theAnalysisModel->getFEs(); 

  while ((elePtr = theEles()) != 0) {
    theSOE->addB(elePtr->getResidual(this) ,elePtr->getID()  );
  }


  (*Residual)=theSOE->getB();

 int size=theAnalysisModel->getNumEqn();
  Matrix dKdh(size,size);
  dKdh.Zero();
  //dKdh=this->getdKdh(gradNumber);
    

   //call the tangent (K)
 //  this->formTangent();


//   Residual->addMatrixVector(1.0,dKdh,*deltaUbar,-1.0);

  // double CallDlambda1dh=this->getLambdaSensitivity(gradNumber);
     
  double CallDlambda1dh = (*dLAMBDAdh)(gradNumber);
  Residual->addVector(1.0,*phat, CallDlambda1dh ); 
  Residual->addVector(1.0,*dphatdh,currentLambda);

  theSOE->setB(*Residual);
  
  //
  // Add the dPext/dh contributions
  //
  static Vector oneDimVectorWithOne(1);
  oneDimVectorWithOne(0) = 1.0;
  static ID oneDimID(1);
  

  Domain *theDomain = theAnalysisModel->getDomainPtr();
  LoadPattern *loadPatternPtr;
  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
  while((loadPatternPtr = thePatterns()) != nullptr) {
    const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(gradNumber);
    int sizeRandomLoads = randomLoads.Size();
    if (sizeRandomLoads == 1) {
      // No random loads in this load pattern
    }
    else {
      // Random loads: add contributions to the 'B' vector
      int numRandomLoads = (int)(sizeRandomLoads/2);
      for (int i=0; i<numRandomLoads*2; i=i+2) {
        Node *aNode;
        DOF_Group *aDofGroup;
        int nodeNumber = (int)randomLoads(i);
        int dofNumber = (int)randomLoads(i+1);
        aNode = theDomain->getNode(nodeNumber);
        aDofGroup = aNode->getDOF_GroupPtr();
        const ID &anID = aDofGroup->getID();
        oneDimID(0) = anID(dofNumber-1);
        theSOE->addB(oneDimVectorWithOne, oneDimID);
      }
    }
    //  (*Residual) =theSOE->getB();
  }

  theSOE->setB(*Residual);
 
  //reset sensitivity flag
  sensitivityFlag = 0;
  return 0;
}


int
MinUnbalDispNorm::saveSensitivity(const Vector &v, int gradNum, int numGrads)
{
   AnalysisModel* theAnalysisModel = this->getAnalysisModel();

   DOF_GrpIter &theDOFGrps = theAnalysisModel->getDOFs();
   DOF_Group 	*dofPtr;

   while ( (dofPtr = theDOFGrps() ) != 0)  {
      dofPtr->saveDispSensitivity(v,gradNum,numGrads);
   }

   return 0;
}

int
MinUnbalDispNorm::saveLambdaSensitivity(double dlambdadh, int gradNum, int numGrads)
{
   AnalysisModel* theAnalysisModel = this->getAnalysisModel();
   Domain *theDomain = theAnalysisModel->getDomainPtr();

   LoadPattern *lpPtr;
   LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
   while ( (lpPtr = thePatterns() ) != 0)
     lpPtr->saveLoadFactorSensitivity(dlambdadh, gradNum, numGrads);

   return 0;
}

   int 
MinUnbalDispNorm::commitSensitivity(int gradNum, int numGrads)
{

   AnalysisModel* theAnalysisModel = this->getAnalysisModel();

   // Loop through the FE_Elements and set unconditional sensitivities
   FE_Element *elePtr;
   FE_EleIter &theEles = theAnalysisModel->getFEs();    
   while((elePtr = theEles()) != 0) {
      elePtr->commitSensitivity(gradNum, numGrads);
   }
   return 0;
}



bool 
MinUnbalDispNorm::computeSensitivityAtEachIteration()
{
   return true;     
}


int 
MinUnbalDispNorm::computeSensitivities(void)
{

  LinearSOE *theSOE = this->getLinearSOE();
  
  // Zero out the old right-hand side of the SOE
  theSOE->zeroB();

  // Form the part of the RHS which are indepent of parameter
  this->formIndependentSensitivityRHS();

  AnalysisModel *theModel = this->getAnalysisModel();   
  Domain *theDomain=theModel->getDomainPtr();//Abbas
  ParameterIter &paramIter = theDomain->getParameters();
  Parameter *theParam;

  // De-activate all parameters
  while ((theParam = paramIter()) != 0)
    theParam->activate(false);

  // Now, compute sensitivity wrt each parameter
  int  numGrads = theDomain->getNumParameters();
  
  paramIter = theDomain->getParameters();

 
  while ((theParam = paramIter()) != 0) {
    // Activate this parameter
    theParam->activate(true);
    
    // Zero the RHS vector
    theSOE->zeroB();
    
    // Get the grad index for this parameter
    int gradIndex = theParam->getGradIndex();
    // this->formResidualDispSensitivity( dUIJdh, gradIndex);


    // Form the RHS


    this->formSensitivityRHS(gradIndex);

    this->formTangent();
	 //theSOE->setX(*dUIJdh);
    theSOE->solve();
    *dUIJdh=theSOE->getX();// sensitivity of the residual displacement

    this->formTangDispSensitivity(dUhatdh,gradIndex);
    double dlamdh = this->getLambdaSensitivity(gradIndex);
    //double dlamdh=(*dLAMBDAdh)(gradIndex);

    //   this->formTangDispSensitivity(dUhatdh,gradIndex);
   
    // To obtain the response sensitivity 
    theSOE->setB(*Residual);

    theSOE->solve();
    

  //  dUIJdh->addVector(1.0, *dUhatdh,dLambdaj);
//    dUIJdh->addVector(1.0, *deltaUhat, Dlambdadh);
//(*sensU)=(*dUIJdh); //theSOE->getX();
  //sensU->addVector(1.0,*dUIJdh,1.0);
    (*sensU) = theSOE->getX();
    //(*sensU) =(*dUIJdh);


    // Save sensitivity to nodes
    this->saveSensitivity( (*sensU), gradIndex, numGrads );
    this->saveLambdaSensitivity(dlamdh, gradIndex, numGrads);
    //dUIJdh->Zero();
    // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
    this->commitSensitivity(gradIndex, numGrads);
    // De-activate this parameter for next sensitivity calc
    theParam->activate(false);
    theSOE->zeroB();//reset the SOE to zero ;Abbas
 
  } 
 
  return 0;
}

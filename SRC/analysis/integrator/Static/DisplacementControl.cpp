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
// Description: This file contains the class definition for DisplacementControl.
// DisplacementControl is an algorithmic class for performing a static analysis
// using the arc length scheme, that is within a load step the following
// constraint is enforced: dU^TdU + alpha^2*dLambda^2 = DisplacementControl^2
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and DisplacementControl is a control parameter.
//
// File: ~/analysis/integrator/DisplacementControl.C
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
#include <DisplacementControl.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <ID.h>
#include <stdlib.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <EquiSolnAlgo.h>
#include <Matrix.h>
#include <TaggedObjectStorage.h>


DisplacementControl::DisplacementControl(int node, int dof, 
                                         double increment, 
                                         Domain *domain,
                                         int numIncr,
                                         double min, 
                                         double max, 
                                         int tang) 
:StaticIntegrator(INTEGRATOR_TAGS_DisplacementControl),
   theNode(node), theDof(dof), theIncrement(increment), theDomain(domain),
   theDofID(-1),
   deltaUhat(0), deltaUbar(0), deltaU(0), deltaUstep(0),dUhatdh(0),
   phat(0), deltaLambdaStep(0.0), currentLambda(0.0), dLambdaStepDh(0.0),dUIJdh(0),Dlambdadh(0.0),
   specNumIncrStep(numIncr), numIncrLastStep(numIncr),
   minIncrement(min), maxIncrement(max),
// sensitivityFlag(0),
// gradNumber(0),
   Residual(0),dlambdadh(0.0),
   dLambda(0.0),  sensU(0),d_deltaU_dh(0),Residual2(0),
   dLAMBDAdh(0),dphatdh(0)
{
  tangFlag = tang;

   // to avoid divide-by-zero error on first update() ensure numIncr != 0
   if (numIncr == 0) {
      opserr << "WARNING DisplacementControl::DisplacementControl() -";
      opserr << " numIncr set to 0, 1 assumed\n";
      specNumIncrStep = 1.0;
      numIncrLastStep = 1.0;
   }

}

DisplacementControl::~DisplacementControl()
{
    // delete any vector object created
   if (deltaUhat != 0)
      delete deltaUhat;
   if (deltaU != 0)
     delete deltaU;
   if (deltaUstep != 0)
     delete deltaUstep;
   if (deltaUbar != 0)
     delete deltaUbar;
   if (phat != 0)
     delete phat;
   if(dUhatdh !=0)
     delete dUhatdh;
   if(dUIJdh !=0)
     delete dUIJdh; 
   if(Residual !=0)
     delete Residual;
   if(sensU !=0)
     delete sensU;
   if(Residual2 !=0)
     delete Residual2;
   if(dLAMBDAdh !=0) 
     delete dLAMBDAdh;
   if(dphatdh !=0)
      delete dphatdh;

   dLAMBDAdh=0;
   dUhatdh=0;
}
 
int
DisplacementControl::newStep()
{

  if (theDofID == -1) {
    opserr << "DisplacementControl::newStep - dof is fixed or constrained (or domainChanged has not been called!)\n";
    return -1;
  }

  // get pointers to AnalysisModel and LinearSOE
  AnalysisModel *theModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();    
  if (theModel == 0 || theLinSOE == 0) {
     opserr << "WARNING DisplacementControl::newStep ";
     opserr << "No AnalysisModel or LinearSOE has been set\n";
     return -1;
  }

  // determine increment for this step
  double gamma = 1.0;
  double factor = pow(specNumIncrStep/numIncrLastStep, gamma);
  theIncrement *= factor;

  if (theIncrement < minIncrement)
     theIncrement = minIncrement;

  else if (theIncrement > maxIncrement)
     theIncrement = maxIncrement;


  // get the current load factor
  currentLambda = theModel->getCurrentDomainTime();

  // determine dUhat
  this->formTangent(tangFlag);
  theLinSOE->setB(*phat);
  if (theLinSOE->solve() < 0) {
     opserr << "DisplacementControl::newStep(void) - failed in solver\n";
     return -1;
  }

  (*deltaUhat) = theLinSOE->getX();
  Vector &dUhat = *deltaUhat;     // this is the Uft in the nonlinear lecture notes
  double dUahat = dUhat(theDofID);// this is the component of the Uft in our nonlinear lecture notes

  if (dUahat == 0.0) {
     opserr << "WARNING DisplacementControl::newStep() ";
     opserr << "dUahat is zero -- zero reference displacement at control node DOF\n";
     return -1;
  }

  // determine delta lambda(1) == dlambda    
  double dlambda = theIncrement/dUahat; // this is the dlambda of the 1st step

 // calldLambda1dh=theIncrement;
  deltaLambdaStep = dlambda;
  currentLambda += dlambda;

  // determine delta U(1) == dU
  (*deltaU) = dUhat;
  (*deltaU) *= dlambda;// this is eq(4) in the paper {dU}_1=dLAmbda1*Uft.
  (*deltaUstep) = (*deltaU);



 if (this->activateSensitivity()==true) { 
   Domain *theDomain=theModel->getDomainPtr();
   ParameterIter &paramIter = theDomain->getParameters();
   Parameter *theParam;

   // De-activate all parameters 
   while ((theParam = paramIter()) != nullptr)
     theParam->activate(false);
   
   // Now, compute sensitivity wrt each parameter
   paramIter = theDomain->getParameters();
   while ((theParam = paramIter()) != nullptr) {
     // Activate this parameter
     theParam->activate(true);
     // Get the grad index for this parameter
     int grad = theParam->getGradIndex();
     this->setGradIndex(grad);
     
     this->formTangDispSensitivity(dUhatdh,grad);

     this->formdLambdaDh(grad);
     theParam->activate(false);
   } 
 }
 ///////////////Abbas/////////////////////////////

  // update model with delta lambda and delta U
  theModel->incrDisp(*deltaU); 
  theModel->applyLoadDomain(currentLambda);
  if (theModel->updateDomain() < 0) {
     opserr << "DisplacementControl::newStep - model failed to update for new dU\n";
     return -1;
  }

  numIncrLastStep = 0;
  return 0;
}

// Update iteration
int DisplacementControl::update(const Vector &dU)
{

   if (theDofID == -1) {
      opserr << "DisplacementControl::update() - domainChanged has not been called\n";
      return -1;
   }
   AnalysisModel *theModel = this->getAnalysisModel();
   LinearSOE *theLinSOE = this->getLinearSOE();    
   if (theModel == 0 || theLinSOE == 0) {
      opserr << "WARNING DisplacementControl::update() ";
      opserr << "No AnalysisModel or LinearSOE has been set\n";
      return -1;
   }

   (*deltaUbar) = dU; // have to do this as the SOE is gonna change
   double dUabar = (*deltaUbar)(theDofID);//dUbar is the vector of residual displacement and dUabar is its component

   // determine dUhat    
   theLinSOE->setB(*phat);
   theLinSOE->solve();
   (*deltaUhat) = theLinSOE->getX();    

   double dUahat = (*deltaUhat)(theDofID);
   if (dUahat == 0.0) {
      opserr << "WARNING DisplacementControl::update() ";
      opserr << "dUahat is zero -- zero reference displacement at control node DOF\n";
      return -1;
   }

   // determine delta lambda(1) == dlambda    
   dLambda = -dUabar/dUahat;// this dLambda i,j

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
   if (theModel->updateDomain() < 0) {
      opserr << "DisplacementControl::update - model failed to update for new dU\n";
      return -1;
   }

   // set the X soln in linearSOE to be deltaU for convergence Test
   theLinSOE->setX(*deltaU);

   numIncrLastStep++;

   return 0;
}



int 
DisplacementControl::domainChanged()
{
   // we first create the Vectors needed
   AnalysisModel *theModel = this->getAnalysisModel();
   LinearSOE *theLinSOE = this->getLinearSOE(); 
    if (theModel == 0 || theLinSOE == 0) {
      opserr << "WARNING DisplacementControl::domainChanged ";
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

   if (phat == 0 || phat->Size() != size) { 
      if (phat != 0)
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

   if (Residual == nullptr || Residual->Size() != size) { 
      if (Residual != nullptr)
         delete Residual;  
      Residual = new Vector(size);
   } 
   
   if (Residual2 == nullptr || Residual2->Size() != size) { 
     if (Residual2 != nullptr)
       delete Residual2;  
     Residual2 = new Vector(size);
   } 

   if (sensU == nullptr || sensU->Size() != size) { 
      if (sensU != nullptr)
         delete sensU;  
      sensU = new Vector(size);
   } 


   Domain *theDomain = theModel->getDomainPtr();//Abbas
   int numGrads = theDomain->getNumParameters();

   if (dLAMBDAdh == 0 || dLAMBDAdh->Size() != (numGrads)) { 
     if (dLAMBDAdh != 0 )  
       delete dLAMBDAdh;
     dLAMBDAdh = new Vector(numGrads);
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

   // check there is a reference load
   int haveLoad = 0;
   for (int i=0; i<size; i++)
      if ( (*phat)(i) != 0.0 ) {
         haveLoad = 1;
         break;
      }

   if (haveLoad == 0) {
      opserr << "WARNING DisplacementControl::domainChanged - zero reference load";
      return -1;
   }

   // lastly we determine the id of the nodal dof
   // TODO: EXTRA CODE TO DO SOME ERROR CHECKING REQUIRED

   Node *theNodePtr = theDomain->getNode(theNode);
   if (theNodePtr == nullptr) {
      opserr << "DisplacementControl::domainChanged - no node\n";
      return -1;
   }

   DOF_Group *theGroup = theNodePtr->getDOF_GroupPtr();
   if (theGroup == 0) {
      return 0;
   }
   const ID &theID = theGroup->getID();
   theDofID = theID(theDof);
   return 0;
}

int
DisplacementControl::sendSelf(int cTag,
      Channel &theChannel)
{
   // TODO: sendSelf
   return 0;
}


int
DisplacementControl::recvSelf(int cTag,
      Channel &theChannel, FEM_ObjectBroker &theBroker)
{
   // TODO: recvSelf
   return 0;
}

void
DisplacementControl::Print(OPS_Stream &s, int flag)
{
   // TODO: Print
}


///////////////////////////Sensitivity Begin///////////////////////////////////
//Added by Abbas
//obtain the derivative of the tangent displacement (dUhatdh)
   Vector *
DisplacementControl::formTangDispSensitivity(Vector *dUhatdh,int gradNumber)
{
   LinearSOE *theLinSOE = this->getLinearSOE(); 
   dUhatdh->Zero();
   dphatdh->Zero();

   // To get the structural stiffness Matrix using the full general system of equations
   //...............................................................
 //  this->formTangent();
 //  Matrix K(size,size);
 //   K.Zero();
  // K=this->ActivateSensitivity();
  // // to print K to the screen
   //...............................................................


  // form dKdh*deltaUhat
  // dUhatdh->addMatrixVector(0.0,dKdh,*deltaUhat,-1.0);
   
   
   //call the tangent (K)
   this->formTangent(tangFlag);
   theLinSOE->setB(*dphatdh);
   if (theLinSOE->solve()<0) {
      opserr<<"SOE failed to obtained dUhatdh ";
      exit(-1);
   }
   (*dUhatdh)=theLinSOE->getX();
   
   // if the parameter is a load parameter.
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

      }
      else {
        // opserr<<"there is sensitivity load parameter"<<endln;//Abbas.............
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


         }
      }
   }

 if(theLinSOE->solve()<0) {
   opserr<<"SOE failed to obtained dUhatdh ";
   exit(-1);
 }

    (*dUhatdh)=theLinSOE->getX();


/////////////////////////////////////////////////////////
   return dUhatdh;
}

// form dLambda for each time step dLambda
double 
DisplacementControl::formdLambdaDh(int gradNumber)
{
  Vector &duHatdh=*dUhatdh;

  //the component of the duHatdh vector
  double duHatdh_Comp = duHatdh(theDofID);

  // call the deltaUhat component
  Vector &UFT=(*deltaUhat);

  double UFT_Comp = UFT(theDofID);// cll the component of the dUhat again
  if(UFT_Comp == 0.0)
    dlambdadh = 0.0;// to avoid dividing by zero
  else
    dlambdadh = -(theIncrement*duHatdh_Comp)/(UFT_Comp*UFT_Comp);

  if(dLAMBDAdh != 0) {
    (*dLAMBDAdh)(gradNumber) = (*dLAMBDAdh)(gradNumber) + dlambdadh;
    return (*dLAMBDAdh)(gradNumber);

  } else {
    return 0.0;
  }
}

 // dLambdadh of the subsequent iterations (J>1)
double 
DisplacementControl::getLambdaSensitivity(int gradNumber)
{

  //call the tangent displacement component (deltaUhat)
  // LinearSOE *theLinSOE = this->getLinearSOE(); 

  //  this->formTangent();
  //  theLinSOE->setB(*phat);
  //   if (theLinSOE->solve() < 0) {
  //      opserr << "DisplacementControl::newStep(void) - failed in solver\n";
  //      return -1;
  //   }

  //  (*deltaUhat) = theLinSOE->getX();
 
   Vector &UFT=*deltaUhat;
   double UFT_Comp=UFT(theDofID);

   //call the derivative of the tangent displacement component (dUftdh).
   Vector &duHatdh=*dUhatdh;
   double duHatdh_Comp=duHatdh(theDofID);

   //Call the residual displacement component
   Vector &dufR=*deltaUbar;
   double dufR_Comp=dufR(theDofID);

   Vector &dufRdh=*dUIJdh;// component of the dUfrDh: derivative of the residual displacement
   double dufRdh_Comp=dufRdh(theDofID);

   if(UFT_Comp==0.0 )
     Dlambdadh=0.0;
   else
     // dLambdadh_ij( the sensitivity of the load component to the parameter h)
     Dlambdadh=((-dufRdh_Comp*UFT_Comp  +(dufR_Comp*duHatdh_Comp)))/(UFT_Comp*UFT_Comp);  //   
 
   // Now update Lambda_ij
   if(dLAMBDAdh !=0) {
     (*dLAMBDAdh)(gradNumber) = (*dLAMBDAdh)(gradNumber)+ Dlambdadh;
     return (*dLAMBDAdh)(gradNumber);
   } else {
     return 0.0;
   }
}


int
DisplacementControl::formIndependentSensitivityRHS()
{
   return 0;
}


   
int
DisplacementControl::formSensitivityRHS(int grad)
{
  this->setResidualType(ResidualType::StaticSensitivity);
  this->setGradIndex(grad);

  // get model
  LinearSOE* theSOE = this->getLinearSOE();
  AnalysisModel* theAnalysisModel = this->getAnalysisModel();

  // Loop through elements
  FE_Element *elePtr;
  FE_EleIter &theEles = theAnalysisModel->getFEs(); 
  while ((elePtr = theEles()) != nullptr)
    theSOE->addB(elePtr->getResidual(this) , elePtr->getID()  );


  (*Residual) = theSOE->getB();

  // needed to calculate dLambdadh
  Residual->addVector(1.0,*phat, (*dLAMBDAdh)(grad) ); 

  Residual->addVector(1.0,*dphatdh, currentLambda ); // needed to calculate dLambdadh
  
  //  Matrix dKdh(size,size);
  //         dKdh.Zero();
  //         dKdh=this->getdKdh();
  //         Residual->addMatrixVector(1.0,dKdh,*deltaUbar,-1.0);
  
  
  Residual2->addVector(1.0,*phat, (*dLAMBDAdh)(grad));// needed to calculate dUdh
  //  Residual2->addVector(1.0,*dphatdh, currentLambda );
  theSOE->setB(*Residual);
  


  // Loop through the loadPatterns and add the dPext/dh contributions
  static Vector oneDimVectorWithOne(1);
  oneDimVectorWithOne(0) = 1.0;
  static ID oneDimID(1);
  
  LoadPattern *loadPatternPtr;
  Domain *theDomain = theAnalysisModel->getDomainPtr();
  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
  while((loadPatternPtr = thePatterns()) != 0) {
    const Vector &randomLoads = loadPatternPtr->getExternalForceSensitivity(grad);
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
        const ID &anID = theDomain->getNode(nodeNumber)
                                  ->getDOF_GroupPtr()
                                  ->getID();
        int relevantID = anID(dofNumber-1);
        oneDimID(0) = relevantID;
        theSOE->addB(oneDimVectorWithOne, oneDimID);
      }
    }
    //  (*Residual) =theSOE->getB();
  }
  
  
  theSOE->setB(*Residual);
  
  //  (*Residual2)=theSOE->getB();  
 
  // reset residual type
  this->setResidualType(ResidualType::StaticUnbalance);

  return 0;
}


int
DisplacementControl::saveSensitivity(const Vector &v, int gradNum, int numGrads)
{
   AnalysisModel* theAnalysisModel = this->getAnalysisModel();

   DOF_GrpIter &theDOFGrps = theAnalysisModel->getDOFs();
   DOF_Group         *dofPtr;

   while ( (dofPtr = theDOFGrps() ) != nullptr)
      dofPtr->saveDispSensitivity(v,gradNum,numGrads);

   return 0;
}

int
DisplacementControl::saveLambdaSensitivity(double dlambdadh, int gradNum, int numGrads)
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
DisplacementControl::commitSensitivity(int gradNum, int numGrads)
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
DisplacementControl::computeSensitivityAtEachIteration()
{
  return true;     
}


int 
DisplacementControl::computeSensitivities()
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
  while ((theParam = paramIter()) != nullptr)
    theParam->activate(false);
  
  // Now, compute sensitivity wrt each parameter
  int  numGrads = theDomain->getNumParameters();
  paramIter = theDomain->getParameters();

 
  while ((theParam = paramIter()) != nullptr) {
    // Activate this parameter
    theParam->activate(true);
    
    // Zero the RHS vector
    theSOE->zeroB();
    
    // Get the grad index for this parameter
    int gradIndex = theParam->getGradIndex();
                  
    // this->formResidualDispSensitivity( dUIJdh, gradIndex);//./.
 
    // Form the RHS
    
    //  this->formTangDispSensitivity(dUhatdh,gradIndex);
    this->formSensitivityRHS(gradIndex);
    this->formTangent(tangFlag);
    theSOE->solve();
    *dUIJdh = theSOE->getX();// sensitivity of the residual displacement
 
    this->formTangDispSensitivity(dUhatdh,gradIndex);
    double dlamdh = this->getLambdaSensitivity(gradIndex);

    // To obtain the response sensitivity 
    theSOE->setB(*Residual);
    theSOE->solve();
       
    //Vector *x=new Vector(size);
    //(*x)=theSOE->getX();
    //  x->addVector(1.0,*deltaUhat,Dlambdadh);
    //  x->addVector(1.0,*dUhatdh,dLambda);
    (*sensU) = theSOE->getX();
    //   sensU->addVector(1.0,*deltaUhat,Dlambdadh);
    //(*sensU) +=(*x);

    // Save sensitivity to nodes
    this->saveSensitivity( (*sensU), gradIndex, numGrads );
    this->saveLambdaSensitivity(dlamdh, gradIndex, numGrads);
    
    // Commit unconditional history variables
    this->commitSensitivity(gradIndex, numGrads);

    // De-activate this parameter for next sensitivity calc
    theParam->activate(false);

    theSOE->zeroB();
 
  } 
  // end of if statement to be run only one time during the iteration process.
  //  CallParam=0;

  return 0;
}

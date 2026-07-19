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
// Description: This file contains the implementation of IncrementalIntegrator.
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#include <IncrementalIntegrator.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Matrix.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <cmath>

IncrementalIntegrator::IncrementalIntegrator(int clasTag)
 : MovableObject(clasTag)
 , statusFlag(CURRENT_TANGENT)
 , eigenVectors(0), eigenValues(0), dampingForces(0)
 , isDiagonal(false),diagMass(0)
 , mV(0),tmpV1(0),tmpV2(0),
   theSOE(0), theAnalysisModel(0), theTest(0)
 // sensitivity 
 , analysisTypeTag(0)
 , SensitivityKey(false)
{
  
}

IncrementalIntegrator::~IncrementalIntegrator()
{
  if (eigenValues != 0)
    delete eigenValues;
  if (eigenVectors != 0)
    delete [] eigenVectors;
  if (dampingForces != 0)
    delete dampingForces;
  if (mV != 0)
    delete mV;
  if (tmpV1 != 0)
    delete tmpV1;
  if (tmpV2 != 0)
    delete tmpV2;
}

void
IncrementalIntegrator::setLinks(AnalysisModel &theModel, LinearSOE &theLinSOE, ConvergenceTest *theConvergenceTest)
{
  theAnalysisModel = &theModel;
  theSOE = &theLinSOE;
  theTest = theConvergenceTest;
}


int 
IncrementalIntegrator::formTangent(int statFlag, double iFact, double cFact)
{
  iFactor = iFact;
  cFactor = cFact;
  return this->formTangent(statFlag);
}

int 
IncrementalIntegrator::formTangent(int statFlag)
{
  int result = 0;
  statusFlag = statFlag;

  if (theAnalysisModel == nullptr || theSOE == nullptr) {
    return -1;
  }

  // zero the A matrix of the linearSOE
  theSOE->zeroA();


  // Nodes
  DOF_GrpIter &theDOFs = theAnalysisModel->getDOFs();
  DOF_Group *dofPtr;
  while ((dofPtr = theDOFs()) != nullptr) {
    if (theSOE->addA(dofPtr->getTangent(this),dofPtr->getID()) <0) {
      opserr << "TransientIntegrator::formTangent() - failed to addA:dof\n";
      result = -1;
    }
  }    

  // the loops to form and add the tangents are broken into two for 
  // efficiency when performing parallel computations - CHANGE

  // loop through the FE_Elements adding their contributions to the tangent
  FE_Element *elePtr;
  FE_EleIter &theEles2 = theAnalysisModel->getFEs();    
  while((elePtr = theEles2()) != nullptr)
    if (theSOE->addA(elePtr->getTangent(this), elePtr->getID()) < 0) {
      opserr << "WARNING IncrementalIntegrator::formTangent -";
      opserr << " failed in addA for ID " << elePtr->getID();            
      result = -3;
    }

  return result;
}


int
IncrementalIntegrator::formIndependentSensitivityLHS(int statFlag)
{
  return this->formTangent(statFlag);
}


int
IncrementalIntegrator::getLastResponse(Vector &result, const ID &id)
{  
  if (theSOE == 0) {
    opserr << "WARNING IncrementalIntegrator::getLastResponse() -";
    opserr << "no LineaerSOE object associated with this object\n";        
    return -1;
  }

  int res = 0; 
  int size = theSOE->getNumEqn() -1;
  const Vector &X = theSOE->getX();
  for (int i=0; i<id.Size(); i++) {
    int loc = id(i);
    if (loc < 0)
      result(i) = 0.0;
    else if (loc <= size) {
      result(i) = X(loc);        
    }
    else {
      opserr << "WARNING IncrementalIntegrator::getLastResponse() -";
      opserr << "location " << loc << "in ID outside bounds ";
      opserr << size << "\n";        
      res = -2;
    }
  }            
  return res;
}


int
IncrementalIntegrator::initialize()
{
  return 0;
}

int
IncrementalIntegrator::commit() 
{
  if (theAnalysisModel == 0) {
    opserr << "WARNING IncrementalIntegrator::commit() -";
    opserr << "no AnalysisModel object associated with this object\n";        
    return -1;
  }    

  return theAnalysisModel->commitDomain();
}


int
IncrementalIntegrator::revertToLastStep() 
{
  return 0;
}   


int
IncrementalIntegrator::revertToStart()
{
  opserr << "ERROR: revertToStart() method not yet implemented " << endln
         << " for the chosen type of integrator. " << endln;
  
  return 0;
}    

LinearSOE *
IncrementalIntegrator::getLinearSOE() const
{
  return theSOE;
}   

ConvergenceTest *
IncrementalIntegrator::getConvergenceTest() const
{
  return theTest;
}   

AnalysisModel *
IncrementalIntegrator::getAnalysisModel() const
{
  return theAnalysisModel;
}


int 
IncrementalIntegrator::formNodalUnbalance()
{
  // loop through the DOF_Groups and add the unbalance
  DOF_GrpIter &theDOFs = theAnalysisModel->getDOFs();
  DOF_Group *dofPtr;
  int res = 0;

  while ((dofPtr = theDOFs()) != nullptr) {
    if (theSOE->addB(dofPtr->getUnbalance(this),dofPtr->getID()) <0) {
        opserr << "WARNING IncrementalIntegrator::formNodalUnbalance -";
        opserr << " failed in addB for ID " << dofPtr->getID();
        res = -2;
    }
  }
      
  return res;
}


int
IncrementalIntegrator::formElementResidual()
{
  // loop through the FE_Elements and add the residual

  int res = 0;    

  FE_Element *elePtr;
  FE_EleIter &theEles2 = theAnalysisModel->getFEs();    
  while ((elePtr = theEles2()) != nullptr) {
    if (theSOE->addB(elePtr->getResidual(this),elePtr->getID()) < 0) {
      opserr << "WARNING IncrementalIntegrator::formElementResidual -";
      opserr << " failed in addB for ID " << elePtr->getID();
      res = -2;
    }
  }

  return res;            
}



int 
IncrementalIntegrator::setupModal(const Vector *modalDampingValues)
{
  int numModes = modalDampingValues->Size();

  const Vector &eigenvalues = theAnalysisModel->getEigenvalues();
  int numEigen = eigenvalues.Size();

  if (numEigen < numModes) 
    numModes = numEigen;

  int numDOF = theSOE->getNumEqn();

  if (eigenValues == 0 || *eigenValues != eigenvalues) {
    if (eigenValues != 0)
      delete eigenValues;
    if (eigenVectors != 0)
      delete [] eigenVectors;
    if (dampingForces != 0)
      delete dampingForces;
    if (mV != 0)
      delete mV;
    if (tmpV1 != 0)
      delete tmpV1;
    if (tmpV2 != 0)
      delete tmpV2;
    
    eigenValues   = new Vector(eigenvalues);
    dampingForces = new Vector(numDOF);
    eigenVectors  = new double[numDOF*numModes];
    mV = new Vector(numDOF);
    tmpV1 = new Vector(numDOF);
    tmpV2 = new Vector(numDOF);

    DOF_GrpIter &theDOFs2 = theAnalysisModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs2()) != 0) { 
      const Matrix &dofEigenvectors =dofPtr->getEigenvectors();
      const ID &dofID = dofPtr->getID();
      for (int j=0; j<numModes; j++) {
        for (int i=0; i<dofID.Size(); i++) {
          int id = dofID(i);
          if (id >= 0) 
            eigenVectors[j*numDOF + id] = dofEigenvectors(i,j);
        }
      }
    }

    double *eigenVectors2 = new double[numDOF*numModes];

    for (int i=0; i<numModes; i++) {
      double *eigenVectorI = &eigenVectors[numDOF*i];    
      double *mEigenVectorI = &eigenVectors2[numDOF*i];    
      Vector v1(eigenVectorI,numDOF);
      Vector v2(mEigenVectorI,numDOF);
      this->doMv(v1, v2);    
    }

    if (eigenVectors != nullptr)
      delete [] eigenVectors;

    eigenVectors = eigenVectors2;
  }

  return 0;
}


int 
IncrementalIntegrator::addModalDampingForce(const Vector *modalDampingValues)
{
  int res = 0;
  
  if (modalDampingValues == nullptr)
    return 0;

  int numModes = modalDampingValues->Size();

  const Vector &eigenvalues = theAnalysisModel->getEigenvalues();
  int numEigen = eigenvalues.Size();

  if (numEigen < numModes) {
    numModes = numEigen;
    opserr << "WARNING: Having to reset numModes to : " << numModes 
           << "as not enough eigenvalues. NOTE if you have done something to require new analysis or have not issued eigen command\n";
  }

  int numDOF = theSOE->getNumEqn();

  if (eigenValues == 0 || *eigenValues != eigenvalues) {
    this->setupModal(modalDampingValues);
  }

  const Vector &vel = this->getVel();

  dampingForces->Zero();

  for (int i=0; i<numModes; i++) {

    double eigenvalue = (*eigenValues)(i);
    double modalDampingValue = (*modalDampingValues)(i);
    if (eigenvalue > 0 && modalDampingValue != 0.0) {
      double wn = sqrt(eigenvalue);

      double *eigenVectorI = &eigenVectors[numDOF*i];
      double beta = 0.0;
      
      for (int j=0; j<numDOF; j++) {
        double eij = eigenVectorI[j];
        if (eij != 0) {
          beta += eij * vel(j);
        }
      }

      beta = -2.0 * modalDampingValue * wn * beta;

      // Fdamp[j] = e[i][j] 
      for (int j=0; j<numDOF; j++) {
        double eij = eigenVectorI[j];
        if (eij != 0)
          (*dampingForces)(j) += beta * eij;
      }
    }
  }

  theSOE->setB(*dampingForces);
  
  return res;
}


const Vector &
IncrementalIntegrator::getVel() {
  opserr << "IncrementalIntegrator::getVel() - not implemented for this integrator\n";
  return theSOE->getX();
}

int 
IncrementalIntegrator::doMv(const Vector &v, Vector &res) {

  int n = v.Size();
  if (isDiagonal == true) {
    for (int i=0; i<n; i++)
      res[i] = diagMass[i]*v[i];
    return 0;
  }

  res.Zero();

  // loop over the FE_Elements
  FE_Element *elePtr;
  FE_EleIter &theEles = theAnalysisModel->getFEs();    
  while((elePtr = theEles()) != nullptr) {
    const Vector &b = elePtr->getM_Force(v, 1.0);
    res.Assemble(b, elePtr->getID(), 1.0);
  }

  // loop over the DOF_Groups
  DOF_Group *dofPtr;
  DOF_GrpIter &theDofs = theAnalysisModel->getDOFs();
  while ((dofPtr = theDofs()) != 0) {
    const Vector &a = dofPtr->getM_Force(v, 1.0);      
    res.Assemble(a, dofPtr->getID(), 1.0);
  }
  return 0;
}

double IncrementalIntegrator::getCFactor()
{
  return 0;
}



int
IncrementalIntegrator::computeSensitivities()
{
  return -1;
}

bool
IncrementalIntegrator::computeSensitivityAtEachIteration()
{
  return false ;
}


bool 
IncrementalIntegrator::shouldComputeAtEachStep()
{
  return (analysisTypeTag == 1);
}


int
IncrementalIntegrator::setGradientType(int flag)
{
  analysisTypeTag = flag;
  return 0;
}


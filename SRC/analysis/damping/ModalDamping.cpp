//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2026, Gustavo A. Araújo R.
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 3-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-3-Clause
//
//===----------------------------------------------------------------------===//
//
// Written: Gustavo A. Araújo R.
//          Claudio M. Perez
//          Barbara Simpson
//          Stanford University
//
// Created: 06/2026
//

// Logging.h provides opserr, and ovoids the noise of OPS_Globals.h
#include <Logging.h>
#include <DOF_Group.h>
#include "ModalDamping.h"
#include <TransientIntegrator.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Matrix.h>
#include <Vector.h>
#include <cmath>
#include <assert.h>
#include <WoodburyUpdate.h>

ModalDamping::ModalDamping(AnalysisModel& theModel, 
                           const Vector& modalDampingValues,
                           int ndf)
 : dampingForces(new Vector(ndf))
 , theAnalysisModel(&theModel)
 , woodbury(nullptr)
 , numDOF(ndf)
 , numModes(modalDampingValues.Size())
 , eigenVectors(new double[ndf*modalDampingValues.Size()])
 , Q(eigenVectors, ndf, modalDampingValues.Size())
 , V(modalDampingValues.Size())
 , X(ndf)
{
  this->setupModal(modalDampingValues);

  woodbury = new WoodburyUpdate(V.Size(), ndf);
}


ModalDamping::~ModalDamping()
{
  if (eigenVectors != 0)
    delete [] eigenVectors;
  if (dampingForces != 0)
    delete dampingForces;
  if (woodbury != nullptr)
    delete woodbury;
}


int
ModalDamping::update(TransientIntegrator &integrator, LinearSOE &theLinSOE)
{
  assert(theLinSOE.getNumEqn() == numDOF);
  return woodbury->formWoodburyBasis(V, Q, integrator.getCFactor(), theLinSOE);
}

const Vector&
ModalDamping::updateX(const Vector& dX, LinearSOE& system)
{
  X = dX;
  this->applyTangent(X);
  system.setX(X);
  return X;
}


int 
ModalDamping::setupModal(const Vector &modalDampingValues)
{
  // int numModes = modalDampingValues->Size();

  const Vector &eigenvalues = theAnalysisModel->getEigenvalues();
  // int numEigen = eigenvalues.Size();

  // if (numEigen < numModes) 
  //   numModes = numEigen;


  if (true) {//currentStamp != eigenStamp) {

    // Store modal damping factors
    const int n_damp = numModes;

    V.resize(n_damp);
    for (int i=0; i<n_damp; i++)
      V(i) = std::sqrt(eigenvalues(i))*2.0*modalDampingValues(i);

    //
    // Store mass-weighted eigenvectors
    //

    double *eigenVectors2 = new double[numDOF*numModes];

    // Store eigenvectors
    DOF_GrpIter &theDOFs2 = theAnalysisModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs2()) != 0) { 
      const Matrix &dofEigenvectors = dofPtr->getEigenvectors();
      const ID &dofID = dofPtr->getID();
      for (int j=0; j<numModes; j++) {
        for (int i=0; i<dofID.Size(); i++) {
          int id = dofID(i);
          if (id >= 0) 
            eigenVectors2[j*numDOF + id] = dofEigenvectors(i,j);
        }
      }
    }


    // Now mass weight the eigenvectors

    for (int i=0; i<numModes; i++) {
      double *eigenVectorI = &eigenVectors2[numDOF*i];    
      double *mEigenVectorI = &eigenVectors[numDOF*i];    
      Vector v1(eigenVectorI,numDOF);
      Vector v2(mEigenVectorI,numDOF);
      theAnalysisModel->applyInertia(v1, v2);    
    }

    delete [] eigenVectors2;

    // eigenVectors = eigenVectors2;
  }

  return 0;
}

int 
ModalDamping::applyTangent(Vector& dX)
{
  return woodbury->applyWoodburyCorrection(Q, dX);
}


int 
ModalDamping::applyResidual(TransientIntegrator &integrator, LinearSOE &theSOE)
{
  int res = 0;

  const Vector *modalDampingValues = theAnalysisModel->getModalDampingFactors();
  
  if (modalDampingValues == nullptr)
    return 0;

  int numModes = V.Size();

  const Vector &eigenvalues = theAnalysisModel->getEigenvalues();
  int numEigen = eigenvalues.Size();

  if (numEigen < numModes) {
    numModes = numEigen;
    opserr << "WARNING: Having to reset numModes to : " << numModes 
           << "as not enough eigenvalues. NOTE if you have done something to require new analysis or have not issued eigen command\n";
  }

  // int numDOF = theSOE->getNumEqn();

  // int currentStamp = 0; // TODO
  // if (currentStamp != eigenStamp)
  //   this->setupModal(modalDampingValues);


  const Vector &vel = integrator.getVel();

  dampingForces->Zero();

  for (int i=0; i<numModes; i++) {

    if (V[i] > 0) {
      // double wn = sqrt(eigenvalue);
      
      // getCFactor?
      double *eigenVectorI = &eigenVectors[numDOF*i];
      double beta = 0.0;
      
      for (int j=0; j<numDOF; j++) {
        double eij = eigenVectorI[j];
        if (eij != 0) {
          beta += eij * vel(j);
        }
      }

      beta = -V[i] * beta;

      // Fdamp[j] = e[i][j] 
      for (int j=0; j<numDOF; j++) {
        double eij = eigenVectorI[j];
        if (eij != 0)
          (*dampingForces)(j) += beta * eij;
      }
    }
  }

  // NOTE: cannot addB, so this must be invoked before anything else 
  // populates theSOE::B.
  theSOE.setB(*dampingForces);
  
  return res;
}


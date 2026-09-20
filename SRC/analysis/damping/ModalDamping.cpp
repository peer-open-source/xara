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
 , pass_solve(false)
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
  pass_solve = true;
  theLinSOE.setForwardUpdate(nullptr);
  int status = woodbury->rebuild(V, Q, integrator.getCFactor(), theLinSOE);
  theLinSOE.setForwardUpdate(this);
  pass_solve = false;
  return status;
}



int 
ModalDamping::setupModal(const Vector &modalDampingValues)
{
  return 0;
}


int 
ModalDamping::solve(const Vector& b, Vector& dX)
{
  return 0;
}


int 
ModalDamping::applyTangent(Vector& dX)
{
  return 0;
}



int 
ModalDamping::apply(const Vector &vel, Vector &res)
{
  return 0;
}

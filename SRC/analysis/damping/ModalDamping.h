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
// Created: 06/2026
//
#pragma once

#include <Integrator.h>
#include <MovableObject.h>
#include <Vector.h>
#include <Matrix.h>


class LinearSOE;
class AnalysisModel;
class FE_Element;
class DOF_Group;
class Vector;
class OPS_Stream;
class TransientIntegrator;
class WoodburyUpdate;

#include <OPS_Stream.h>
#include <LinearSOE.h>


class ModalDamping
{
  public:
    typedef int TangentFlagType; 

    ModalDamping(AnalysisModel& theModel, 
                 const Vector& modalDampingValues,
                 int ndf);

    virtual ~ModalDamping();

    int update(TransientIntegrator &, LinearSOE&);
    int applyTangent(Vector& dX);
    int applyResidual(TransientIntegrator &, LinearSOE &);

    const Vector& updateX(const Vector& dX, LinearSOE& system);
    
    void Print(OPS_Stream& s, int flag) const {
      s << "ModalDamping(" << numDOF << ", " << numModes << ")\n";
    }

  private:

    // Setup Q and V
    int setupModal(const Vector &modalDampingValues);

    double   *const eigenVectors;
    Vector   *dampingForces;
    WoodburyUpdate *woodbury;

    Vector V; // 2*modalDampingValues*eigenvalues
    Matrix Q; // mass-weighted eigenvectors
    Vector X;


    AnalysisModel *theAnalysisModel;

    const int numDOF;
    const int numModes;
};


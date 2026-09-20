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
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Written: Gustavo A. Araújo R.
// Created: 06/2026
//
#pragma once

#include <Matrix.h>
#include <Vector.h>
class LinearSOE;
class OPS_Stream;

class WoodburyUpdate
{
  public:
    explicit WoodburyUpdate(int n_damp, int n_dof);
    ~WoodburyUpdate();

    int rebuild(const Vector& V,
                Matrix& Q,
                double cFactor,
                LinearSOE&);

    int applyWoodburyCorrection(const Matrix& Qmat, 
                                Vector& ddX);


  private:

    const int numDOF;
    const int numModes;

    Matrix Z;
    Matrix G;

    Vector *workV1;
    Vector *workV2;
};


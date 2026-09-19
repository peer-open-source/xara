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
//
//
// Written: Gustavo A. Araújo R.
// Created: 06/2026
//
#include <Logging.h>
#include <WoodburyUpdate.h>
#include <Domain.h>
#include <Matrix.h>
#include <Vector.h>
#include <LinearSOE.h>

#include <assert.h>
#include <cstring>

#include <OPS_Stream.h>

WoodburyUpdate::WoodburyUpdate(int n_damp, int n_dof)
  : numDOF(n_dof),
    numModes(n_damp),
    Z(n_dof, n_damp),
    G(n_damp, n_damp),
    workV1(new Vector(n_damp)),
    workV2(new Vector(n_damp))
{

}


WoodburyUpdate::~WoodburyUpdate()
{
  if (workV1 != nullptr)
    delete workV1;

  if (workV2 != nullptr)
    delete workV2;
}



int
WoodburyUpdate::rebuild(const Vector& V,
                        Matrix& Qmat,
                        double cFactor,
                        LinearSOE& theSOE)
{
  return 0;
}


int
WoodburyUpdate::applyWoodburyCorrection(const Matrix& Qmat,
                                        Vector& ddX)
{
  return 0;
}

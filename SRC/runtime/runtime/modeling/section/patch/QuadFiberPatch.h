//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Claudio M. Perez
//
// Adapted from: QuadFiberPatch.h
// Written by Remo M. de Souza
// December 1998
//
#pragma once

#include <MatrixND.h>
#include "FiberPatch.h"
#include <Vector.h>

class Matrix;

namespace OpenSees {

class Cell;

class QuadFiberPatch : public FiberPatch {
public:
  QuadFiberPatch();
  QuadFiberPatch(int materialID, int numSubdivIJ, int numSubdivJK, const MatrixND<4,2>& vertexCoords);

  ~QuadFiberPatch();

  int getMaterialID() const;
  int getNumCells() const;
  FiberCell** getCells() const;

protected:
private:
  int matID;
  int nDivIJ, nDivJK;
  MatrixND<4,2> vertCoord;
};

} // namespace OpenSees

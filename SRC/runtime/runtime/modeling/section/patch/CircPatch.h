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
// Adapted from CircPatch.h
// Written by Remo M. de Souza
// December 1998
//
#ifndef CircPatch_h
#define CircPatch_h

#include "FiberPatch.h"
#include <Vector.h>
#include <VectorND.h>

class Matrix;

namespace OpenSees {
class FiberCell;
class CircPatch : public FiberPatch {
public:
  CircPatch(int material,
            int numSubdivCircunf, int numSubdivRadial,
            const VectorND<2>& centerPosition, 
            double internRadius, double externRadius,
            double initialAngle, double finalAngle);

  ~CircPatch();

  int getMaterialID() const;
  int getNumCells() const;
  FiberCell** getCells() const;

private:
  int matID;
  int nDivCirc, nDivRad;
  const VectorND<2> centerPosit;
  double intRad, extRad;
  double initAng, finalAng;
};
} // namespace OpenSees
#endif

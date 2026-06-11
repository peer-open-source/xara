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
// File: StraightFiberLayer.h
// Written by Remo M. de Souza
// December 1998
//
#pragma once

#include "FiberLayer.h"
#include <VectorND.h>
#include <FiberCell.h>

namespace OpenSees {

class StraightFiberLayer : public FiberLayer {
public:
  StraightFiberLayer(int materialID, int numReinfBars, 
                     double reinfBarArea,
                     const VectorND<2>& initialPosition, 
                     const VectorND<2>& finalPosition);

  virtual ~StraightFiberLayer() {};

  int getNumReinfBars() const;
  std::vector<FiberCell> getReinfBars() const;

private:
  int nReinfBars;
  VectorND<2> initPosit;
  VectorND<2> finalPosit;
};

} // namespace OpenSees

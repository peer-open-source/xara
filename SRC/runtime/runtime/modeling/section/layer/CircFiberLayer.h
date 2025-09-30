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
// File: CircFiberLayer.h
// Written by Remo M. de Souza
// December 1998
//
#pragma once

#include <FiberLayer.h>
#include <VectorND.h>
#include <FiberCell.h>

namespace OpenSees {

class CircFiberLayer : public FiberLayer {
public:
  // Constructor for an arc
  CircFiberLayer(int material, int n, double area,
                 const VectorND<2>& center, double radius, 
                 double initialAngle,
                 double finalAngle);
  // Constructor for full circle
  CircFiberLayer(int material, int n, double area,
                 const VectorND<2>& center, double radius);

  virtual ~CircFiberLayer() {};

  int getNumReinfBars() const;
  std::vector<FiberCell> getReinfBars() const;

private:
  int nReinfBars;
  VectorND<2> centerPosit;
  double arcRad;
  double initAng;
  double finalAng;
};
} // namespace OpenSees

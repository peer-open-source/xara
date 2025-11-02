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
// File: FiberLayer.h
// Written by Remo M. de Souza
// December 1998
//
#pragma once

#include <vector>
#include <FiberCell.h>

namespace OpenSees {
class FiberLayer {
public:
  FiberLayer(int material, double area) : material(material), area(area) {}
  virtual ~FiberLayer() {};

  int getMaterialID() const {
    return material;
  }

  virtual int getNumReinfBars() const = 0;
  virtual std::vector<FiberCell> getReinfBars() const = 0;

protected:
  double getCellArea() const {
    return area;
  }
private:
  int material;
  double area;
};
} // namespace OpenSees

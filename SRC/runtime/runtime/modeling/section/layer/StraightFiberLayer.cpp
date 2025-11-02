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
// File: StraightFiberLayer.C
// Written by Remo M. de Souza
// December 1998
//
#include <VectorND.h>

#include <FiberCell.h>
#include "StraightFiberLayer.h"
using namespace OpenSees;

StraightFiberLayer::StraightFiberLayer(int material, int n, double area,
                                       const VectorND<2>& xi, 
                                       const VectorND<2>& xj)
 : FiberLayer(material, area),
   nReinfBars(n),
   initPosit(xi),
   finalPosit(xj)
{
}

int
StraightFiberLayer::getNumReinfBars() const
{
  return nReinfBars;
}

std::vector<FiberCell> 
StraightFiberLayer::getReinfBars() const
{
  VectorND<2> barPosit;
  std::vector<FiberCell> bars(nReinfBars);

  if (nReinfBars == 1) {
    VectorND<2> location {
      barPosit(0) = (initPosit(0) + finalPosit(0))/2.0,
      barPosit(1) = (initPosit(1) + finalPosit(1))/2.0
    };
    bars[0] = FiberCell(this->getMaterialID(), this->getCellArea(), location);
  }

  else if (nReinfBars > 1) {
    double dy = (finalPosit(0) - initPosit(0)) / (nReinfBars - 1);
    double dz = (finalPosit(1) - initPosit(1)) / (nReinfBars - 1);

    // reinfBars = new ReinfBar[nReinfBars];

    for (int i = 0; i < nReinfBars; i++) {
      VectorND<2> location {
         initPosit(0) + dy * i,
         initPosit(1) + dz * i
      };

      bars[i] = FiberCell(this->getMaterialID(), this->getCellArea(), location);
    }
  }

  return bars;
}

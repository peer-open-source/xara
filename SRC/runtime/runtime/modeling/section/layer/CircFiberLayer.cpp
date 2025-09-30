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
// File: CircFiberLayer.C
// Written by Remo M. de Souza
// December 1998
//
#include <VectorND.h>
#include <cmath>
#include <string>

#include "CircFiberLayer.h"
#include <FiberCell.h>
using namespace OpenSees;


CircFiberLayer::CircFiberLayer(int material, 
                               int numReinfBars, 
                               double area,
                               const VectorND<2>& centerPosition, 
                               double arcRadius, double initialAngle,
                               double finalAngle)
 : FiberLayer(material, area),
   nReinfBars(numReinfBars),
   centerPosit(centerPosition),
   arcRad(arcRadius),
   initAng(initialAngle),
   finalAng(finalAngle)
{
}

CircFiberLayer::CircFiberLayer(int material, int numReinfBars, double area,
                               const VectorND<2>& centerPosition, double radius)
 : FiberLayer(material, area),
   nReinfBars(numReinfBars),
   centerPosit(centerPosition),
   arcRad(radius),
   initAng(0.0),
   finalAng(0.0)
{
  // Figure out final angle so that complete circle does not put
  // two bars at the same location
  if (nReinfBars > 0)
    finalAng = 360.0 - 360.0 / nReinfBars;
}

int
CircFiberLayer::getNumReinfBars() const
{
  return nReinfBars;
}

std::vector<FiberCell>
CircFiberLayer::getReinfBars() const
{
  std::vector<FiberCell> bars(nReinfBars);

  double pi = std::acos(-1.0);

  if (nReinfBars > 0) {
    double initAngRad  = pi * initAng / 180.0;
    double finalAngRad = pi * finalAng / 180.0;

    double dtheta;
    if (nReinfBars > 1)
      dtheta = (finalAngRad - initAngRad) / (nReinfBars - 1);
    else
      dtheta = 0.0; // Doesn't really matter what this is


    for (int i = 0; i < nReinfBars; i++) {
      double theta = initAngRad + dtheta * i;
      VectorND<2> position {
           centerPosit(0) + arcRad * cos(theta),
           centerPosit(1) + arcRad * sin(theta)
      };
      bars[i] = FiberCell(this->getMaterialID(), this->getCellArea(), position);
    }
  }

  return bars;
}

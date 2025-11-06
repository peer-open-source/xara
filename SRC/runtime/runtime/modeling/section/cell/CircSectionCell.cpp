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
// Written: fmk
// March 2014

#include <cmath>
#include <CircSectionCell.h>

using namespace OpenSees;

CircSectionCell::CircSectionCell(double r1, double r2, 
                                 double alpha, double theta, 
                                 double offsetX,
                                 double offsetY)
{

  double a = alpha / 2.0;
  double At =   a * r2*r2;
  double ct = 2.0 * r2*std::sin(a)/(3.0 * a);
  double A1 =   a * r1*r1;
  double c1 = 2.0 * r1*std::sin(a)/(3.0 * a);

  area     = At - A1;
  double c = (At * ct - A1 * c1) / area;

  location[0] = std::cos(theta) * c + offsetX;
  location[1] = std::sin(theta) * c + offsetY;
}
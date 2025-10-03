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
//
// Adapted from Remo M. de Souza's QuadCell.C
// December 1998
//
#include <MatrixND.h>
#include "QuadFiberCell.h"
using namespace OpenSees;


QuadFiberCell::QuadFiberCell(const MatrixND<4,2>& vertCoord)
{

  // double x0 = vertCoord(0, 0);
  // double y0 = vertCoord(0, 1);
  // double x1 = vertCoord(1, 0);
  // double y1 = vertCoord(1, 1);
  // double x2 = vertCoord(2, 0);
  // double y2 = vertCoord(2, 1);
  // double x3 = vertCoord(3, 0);
  // double y3 = vertCoord(3, 1);

  // double area = ((x2 - x1) * (y0 - y1) - (x0 - x1) * (y2 - y1) + (x0 - x3) * (y2 - y3) -
  //                (x2 - x3) * (y0 - y3)) /
  //               2.0;

  area = 0;

  for (int i = 0; i < 4; i++) {
    int i1 = (i + 1) % 4;
    double yi     = vertCoord(i, 0);
    double zi     = vertCoord(i, 1);
    double yi1    = vertCoord(i1, 0);
    double zi1    = vertCoord(i1, 1);

    area += (zi1 - zi) * (yi1 + yi);
  }
  area /= 2.0;

  //
  // centroid
  //

  double CGy = 0.0, CGz = 0.0;

  for (int i = 0; i < 4; i++) {
    int i1 = (i + 1) % 4;

    double yi  = vertCoord(i, 0);
    double zi  = vertCoord(i, 1);
    double yi1 = vertCoord(i1, 0);
    double zi1 = vertCoord(i1, 1);

    double dyi = yi1 - yi;
    double dzi = zi1 - zi;

    double integ = yi * zi + (yi * dzi + zi * dyi) / 2.0 + dyi * dzi / 3.0;

    CGy -= dyi * integ;
    CGz += dzi * integ;
  }

  CGy /= area;
  CGz /= area;

  location[0] = CGy;
  location[1] = CGz;

}

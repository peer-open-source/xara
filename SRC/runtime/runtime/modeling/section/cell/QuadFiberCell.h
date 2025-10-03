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
// File: QuadCell.h
//
// Written by Remo M. de Souza
// December 1998
#pragma once

#include "FiberCell.h"
#include <MatrixND.h>

namespace OpenSees {

class QuadFiberCell : public FiberCell {
public:
  QuadFiberCell(const MatrixND<4,2>& vertexCoords);
};

} // namespace OpenSees

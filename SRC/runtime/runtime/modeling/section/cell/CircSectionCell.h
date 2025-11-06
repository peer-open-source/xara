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
#pragma once

#include <FiberCell.h>

namespace OpenSees {

class CircSectionCell : public FiberCell {
public:
  CircSectionCell(double r2, double r1, double alpha, double theta, double centerX, double centerY);
};

} // namespace OpenSees

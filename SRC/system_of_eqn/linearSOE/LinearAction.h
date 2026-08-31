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
#include <Vector.h>
#include <LinearSOE.h>


class LinearAction
{
  public:
    virtual ~LinearAction() = default;

    virtual int apply(const Vector& x, Vector& b) = 0;
    virtual int solve(const Vector& b, Vector& x) = 0;
};

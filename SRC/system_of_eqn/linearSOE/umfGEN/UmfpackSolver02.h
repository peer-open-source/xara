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
// Description: This file contains the class definition for
// UmfpackSolver02.
// It solves an SparseGenCSC by calling UMFPACK.
//
#pragma once

#include <LinearSOESolver.h>
#include "../../../../OTHER/UMFPACK/umfpack.h"

#include <array>

class Channel;
class FEM_ObjectBroker;
class SparseGenCSC;
class Vector;

class UmfpackSolver02 : public LinearSOESolver
{
public:
    explicit UmfpackSolver02(bool doDet = false);
    ~UmfpackSolver02() override;

    UmfpackSolver02(const UmfpackSolver02 &) = delete;
    UmfpackSolver02 &operator=(const UmfpackSolver02 &) = delete;
    UmfpackSolver02(UmfpackSolver02 &&) = delete;
    UmfpackSolver02 &operator=(UmfpackSolver02 &&) = delete;

    int solve() override;
    int solve(const Vector &, Vector &) override;
    int setSize() override;

    int setLinearSOE(SparseGenCSC &);

    double getDeterminant() override;
    bool requireDeterminant() override { doDeterminant = true; return true; }

private:
    void clearSymbolic();

    void *symbolic;
    std::array<double, UMFPACK_CONTROL> control;
    std::array<double, UMFPACK_INFO> info;
    SparseGenCSC *theSOE;
    double determinant;
    bool doDeterminant;
};

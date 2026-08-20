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

    int sendSelf(int commitTag, Channel &) override;
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;

    double getDeterminant() override;

private:
    void clearSymbolic();

    void *symbolic;
    std::array<double, UMFPACK_CONTROL> control;
    std::array<double, UMFPACK_INFO> info;
    SparseGenCSC *theSOE;
    double determinant;
    bool doDeterminant;
};

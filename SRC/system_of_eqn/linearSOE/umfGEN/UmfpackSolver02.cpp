//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2026, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Created: 2026-07
//
#include <UmfpackSolver02.h>
#include <SparseGenCSC.h>

#include <Channel.h>
#include <Constants.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>

UmfpackSolver02::UmfpackSolver02(bool doDet)
  : LinearSOESolver(-1), symbolic(nullptr),
    control(), info(), theSOE(nullptr), determinant(0.0), doDeterminant(doDet)
{
}

UmfpackSolver02::~UmfpackSolver02()
{
  this->clearSymbolic();
}

int
UmfpackSolver02::solve()
{
  if (theSOE == nullptr)
    return -1;

  return this->solve(theSOE->getB(), theSOE->X);
}

int
UmfpackSolver02::solve(const Vector &B, Vector &X)
{
  if (theSOE == nullptr)
    return -1;

  const int size = X.Size();
  const int nnz = static_cast<int>(theSOE->getRowIndices().size());
  if (size == 0 || nnz == 0)
    return 0;

  if (symbolic == nullptr)
    return -1;

  if (B.Size() != size || static_cast<int>(theSOE->getPointers().size()) != size + 1)
    return -1;

  void *numeric = nullptr;
  int status = umfpack_di_numeric(theSOE->getPointers().data(), 
                                  theSOE->getRowIndices().data(),
                                  theSOE->getValues().data(), 
                                  symbolic, 
                                  &numeric,
                                  control.data(), info.data());
  if (status != UMFPACK_OK) {
    if (numeric != nullptr)
      umfpack_di_free_numeric(&numeric);
    return -1;
  }

  status = umfpack_di_solve(UMFPACK_A, 
                            theSOE->getPointers().data(), 
                            theSOE->getRowIndices().data(),
                            theSOE->getValues().data(), 
                            &X(0),
                            &B(0),
                            numeric,
                            control.data(), info.data());

  if (doDeterminant)
    umfpack_di_get_determinant(&determinant, nullptr, numeric, info.data());

  umfpack_di_free_numeric(&numeric);
  return status == UMFPACK_OK ? 0 : -1;
}

int
UmfpackSolver02::setSize()
{
  this->clearSymbolic();
  umfpack_di_defaults(control.data());
  control[UMFPACK_PIVOT_TOLERANCE] = 1.0;
  control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;

  if (theSOE == nullptr) {
    return -1;
  }

  const int size = theSOE->getX().Size();
  const int nnz = static_cast<int>(theSOE->Ai.size());
  if (size == 0 || nnz == 0)
    return 0;

  if (static_cast<int>(theSOE->getPointers().size()) != size + 1)
    return -1;

  const int status = umfpack_di_symbolic(size, size,
                                         theSOE->getPointers().data(),
                                         theSOE->getRowIndices().data(),
                                         theSOE->getValues().data(),
                                         &symbolic, 
                                         control.data(), 
                                         info.data());
  if (status != UMFPACK_OK) {
    this->clearSymbolic();
    return -1;
  }

  return 0;
}

int
UmfpackSolver02::setLinearSOE(SparseGenCSC &soe)
{
  theSOE = &soe;
  return 0;
}


double
UmfpackSolver02::getDeterminant()
{
  return doDeterminant ? determinant : OpenSees::Constants::nan;
}

void
UmfpackSolver02::clearSymbolic()
{
  if (symbolic != nullptr) {
    umfpack_di_free_symbolic(&symbolic);
  }
}

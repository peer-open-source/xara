//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2026, Gustavo A. Araújo R.
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Written: Gustavo A. Araújo R.
// Created: 06/2026
//
#ifndef WoodburyUpdate_h
#define WoodburyUpdate_h


class LinearSOE;
class Matrix;
class Vector;
class OPS_Stream;


// Wraps a structural LinearSOE for modalDampingW: 
//   A_s in inner matrix;
//   A_modal = Q*diag(dbar)*Q^T
// applied in solve (Woodbury) and formAp. See WoodburyUpdate.cpp for formulas.
class WoodburyUpdate
{
  public:
    explicit WoodburyUpdate(int n_damp, int n_dof);
    ~WoodburyUpdate();


    // Build Z = A_s^{-1}*Q and G = Q^T*Z + diag(1/dbar) for current A_s.
    int formWoodburyBasis(const Vector& V,
                          Matrix& Q,
                          double cFactor,
                          LinearSOE&);

    // x = x_s - Z * G^{-1} * Q^T * x_s  (inner must hold x_s = A_s^{-1} b)
    int applyWoodburyCorrection(const Matrix& Qmat, 
                                // const Vector& dX, 
                                Vector& ddX);

    void Print(OPS_Stream& s, int flag) const {
      s << "Woodbury()\n";
    }


  private:

    const int numDOF;
    const int numModes;

    Matrix *woodburyZ;     // Z = A_s^{-1} * Q
    Matrix *woodburyG;     // G = Q^T * Z + diag(1/dbar)

    Vector *workV1;        // scratch (modal-sized)
    Vector *workV2;        // y = G^{-1} * Q^T * x
};


#endif

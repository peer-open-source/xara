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
#include <Logging.h>
#include <WoodburyUpdate.h>
#include <Domain.h>
#include <Matrix.h>
#include <Vector.h>
#include <LinearSOE.h>

// #include <math.h>
#include <assert.h>
#include <cstring>

#include <OPS_Stream.h>
#define LinSOE_TAGS_WoodburyUpdate 0

// modalDampingW: A x = b with  A = A_s + Q D Q^T  (A_s assembled in inner SOE only).
//
//   Q = M Phi   (mass-weighted mode shapes; columns with xi_i = 0 omitted)
//   D = diag(dbar),   dbar_i = c_2 * 2*xi_i*omega_i,   omega_i = sqrt(lambda_i)
//
//   x_s = A_s^{-1} b
//   Z = A_s^{-1} Q
//   G = D^{-1} + Q^T Z  
//   x = x_s - Z * G^{-1} * Q^T * x_s
//
// formAp: A*p = A_s*p + Q*D*Q^T*p.  Z,G rebuilt when A_s changes (zeroA/addA).

WoodburyUpdate::WoodburyUpdate(int n_damp, int n_dof)
    : numDOF(n_dof),
      numModes(n_damp),
      woodburyZ(new Matrix(n_dof, n_damp)),
      woodburyG(new Matrix(n_damp, n_damp)),
      workV1(new Vector(n_damp)),
      workV2(new Vector(n_damp))
{

}


WoodburyUpdate::~WoodburyUpdate()
{

    if (woodburyZ != nullptr) {
        delete woodburyZ;
        woodburyZ = nullptr;
    }
    if (woodburyG != nullptr) {
        delete woodburyG;
        woodburyG = nullptr;
    }
    if (workV1 != nullptr) {
        delete workV1;
        workV1 = nullptr;
    }
    if (workV2 != nullptr) {
        delete workV2;
        workV2 = nullptr;
    }
}



// formWoodburyBasis
// Form Z and G for the current A_s (call after tangent assembly; cFactor from integrator).
int
WoodburyUpdate::formWoodburyBasis(const Vector& V,
                                   Matrix& Qmat,
                                   double cFactor,
                                   LinearSOE& innerSOE)
{


    assert(numDOF == innerSOE.getNumEqn());
    assert(Qmat.noRows() == numDOF);
    assert(Qmat.noCols() == numModes);

    const int n_damp = V.Size();

    if (numDOF <= 0 || n_damp <= 0)
        return -1;

    // No modal tangent contribution this assembly.
    if (cFactor == 0.0)
        return 0;

    // Q stored column-wise in eigenVectors; Z in woodburyZ (numDOF x numModes).
    /*const*/ double *eigenVectors = &Qmat(0,0);

    const Vector bSave(innerSOE.getB());
    int activeCol = 0;
    for (int i = 0; i < n_damp; ++i) {
        if (V[i] <= 0.0)
            continue;

        const Vector qcol(&eigenVectors[activeCol * numDOF], numDOF);
        // Z(:,activeCol) = A_s^{-1} * q_col
        if (innerSOE.setB(qcol) < 0) {
            innerSOE.setB(bSave);
            return -2;
        }
        if (innerSOE.solve() < 0) {
            innerSOE.setB(bSave);
            return -3;
        }
        Vector &xcol = const_cast<Vector &>(innerSOE.getX());
        std::memcpy(&(*woodburyZ)(0, activeCol), &xcol(0),
                    static_cast<size_t>(numDOF) * sizeof(double));

        ++activeCol;
    }
    innerSOE.setB(bSave);

    if (activeCol <= 0)
        return 0;


    assert(numModes == activeCol);

    // G = Q^T * Z + diag(1/dbar)
    woodburyG->Zero();
    if (woodburyG->addMatrixTransposeProduct(0.0, Qmat, *woodburyZ, 1.0) < 0)
        return -4;

    for (int i = 0; i < numModes; ++i)
        (*woodburyG)(i, i) += 1.0 / (V[i]*cFactor);


    return 0;
}


// Given inner x_s = A_s^{-1} b in innerSOE, set vectX = A_eff^{-1} b:
//   r = Q^T * x_s,   y = G^{-1} * r,   x = x_s - Z * y
int
WoodburyUpdate::applyWoodburyCorrection(const Matrix& Qmat, 
                                        // const Vector& dX, 
                                        Vector& ddX)
{
    assert(numDOF == ddX.Size());
    assert(Qmat.noRows() == numDOF);
    assert(Qmat.noCols() == numModes);

    // const Vector &innerX = innerSOE->getX();
    // if (innerX.Size() != vectX->Size()) {
    //     opserr << "WARNING WoodburyUpdate::applyWoodburyCorrection() - inner/wrapper X size mismatch\n";
    //     return -1;
    // }

//
//
//

    // No positive eigenmodes: solve A_s x = b only
    if (numModes <= 0)
        return 0;

    // workV1 = Q^T * x_s
    if (workV1->addMatrixTransposeVector(0.0, Qmat, ddX, 1.0) < 0)
        return -1;

    // workV2 = G^{-1} * workV1 = G^{-1} * Q^T * x_s
    if (woodburyG->Solve(*workV1, *workV2) < 0)
        return -1;

    // x = x_s - Z * workV2 = x_s - Z * G^{-1} * Q^T * x_s
    // ddX = dX;
    if (ddX.addMatrixVector(1.0, *woodburyZ, *workV2, -1.0) < 0)
        return -1;

    return 0;
}




// int
// WoodburyUpdate::saveSparseA(OPS_Stream &output, int baseIndex)
// {
//     if (innerSOE == nullptr)
//         return -1;
//     return innerSOE->saveSparseA(output, baseIndex);
// }

// int
// WoodburyUpdate::getSparseA(ID &rowIndices, ID &colIndices, Vector &values,
//                            int baseIndex)
// {
//     if (innerSOE == nullptr)
//         return -1;
//     return innerSOE->getSparseA(rowIndices, colIndices, values, baseIndex);
// }


#if 0

// Ap = A_eff * p = A_s*p + Q*diag(dbar)*Q^T*p
int
WoodburyUpdate::formAp(const Vector &p, Vector &Ap)
{
    int res = innerSOE->formAp(p, Ap);
    if (res < 0)
        return res;
    return applyModalMatvec(p, Ap);
}

// A*p += Q * diag(dbar) * Q^T * p  (caller must have set A*p = A_s*p via inner formAp first)
int
WoodburyUpdate::applyModalMatvec(const Vector &p, Vector &Ap)
{
    if (!basisValid || numModes <= 0)
        return 0;

    if (workV1 == nullptr)
        return -1;

    Matrix Qmat(eigenVectors, numDOF, numModes);
    // workV1 = Q^T * p
    if (workV1->addMatrixTransposeVector(0.0, Qmat, p, 1.0) < 0)
        return -1;

    // workV1 = diag(dbar) * Q^T * p
    for (int i = 0; i < numModes; ++i)
        (*workV1)(i) *= dBar[i];

    // A*p += Q * workV1 = A*p + Q * diag(dbar) * Q^T * p
    if (Ap.addMatrixVector(1.0, Qmat, *workV1, 1.0) < 0)
        return -1;

    return 0;
}

#endif


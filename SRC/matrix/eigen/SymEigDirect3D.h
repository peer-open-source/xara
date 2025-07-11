// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2025
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 8.0.2025.05.10
//
#pragma once

// The document
// https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
// describes algorithms for solving the eigensystem associated with a 3x3
// symmetric real-valued matrix.
//
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

#include <Vector3D.h>

namespace OpenSees {

template <typename T> 
class SortEigenstuff {
public:
  void
  operator()(int32_t sortType, bool isRotation, 
             Vector3D& eval,
             std::array<Vector3D, 3>& evec)
  {
    if (sortType != 0) {
      // Sort the eigenvalues to eval[0] <= eval[1] <= eval[2].
      std::array<size_t, 3> index{};
      if (eval[0] < eval[1]) {
        if (eval[2] < eval[0]) {
          // even permutation
          index[0] = 2;
          index[1] = 0;
          index[2] = 1;
        } else if (eval[2] < eval[1]) {
          // odd permutation
          index[0]   = 0;
          index[1]   = 2;
          index[2]   = 1;
          isRotation = !isRotation;
        } else {
          // even permutation
          index[0] = 0;
          index[1] = 1;
          index[2] = 2;
        }
      } else {
        if (eval[2] < eval[1]) {
          // odd permutation
          index[0]   = 2;
          index[1]   = 1;
          index[2]   = 0;
          isRotation = !isRotation;
        } else if (eval[2] < eval[0]) {
          // even permutation
          index[0] = 1;
          index[1] = 2;
          index[2] = 0;
        } else {
          // odd permutation
          index[0]   = 1;
          index[1]   = 0;
          index[2]   = 2;
          isRotation = !isRotation;
        }
      }

      if (sortType == -1) {
        // The request is for eval[0] >= eval[1] >= eval[2]. This
        // requires an odd permutation, (i0,i1,i2) -> (i2,i1,i0).
        std::swap(index[0], index[2]);
        isRotation = !isRotation;
      }

      Vector3D unorderedEVal                = eval;
      std::array<Vector3D, 3> unorderedEVec = evec;
      for (size_t j = 0; j < 3; ++j) {
        size_t i = index[j];
        eval[j]  = unorderedEVal[i];
        evec[j]  = unorderedEVec[i];
      }
    }

    // Ensure the ordered eigenvectors form a right-handed basis.
    if (!isRotation) {
      for (size_t j = 0; j < 3; ++j) {
        evec[2][j] = -evec[2][j];
      }
    }
  }
};


template <typename T, int32_t sortType> 
class SymEigDirect3D {
public:
  // The input matrix must be symmetric, so only the unique elements
  // must be specified: a00, a01, a02, a11, a12, and a22.  The
  // eigenvalues are sorted in ascending order: eval0 <= eval1 <= eval2.

  void
  operator()(T a00, T a01, T a02, T a11, T a12, T a22, 
             Vector3D& eval,
             std::array<Vector3D, 3>& evec) const
  {
    // Precondition the matrix by factoring out the maximum absolute
    // value of the components.  This guards against floating-point
    // overflow when computing the eigenvalues.
    T max0          = std::max(std::fabs(a00), std::fabs(a01));
    T max1          = std::max(std::fabs(a02), std::fabs(a11));
    T max2          = std::max(std::fabs(a12), std::fabs(a22));
    T maxAbsElement = std::max(std::max(max0, max1), max2);
    if (maxAbsElement == (T)0) {
      // A is the zero matrix.
      eval[0] = (T)0;
      eval[1] = (T)0;
      eval[2] = (T)0;
      evec[0] = {(T)1, (T)0, (T)0};
      evec[1] = {(T)0, (T)1, (T)0};
      evec[2] = {(T)0, (T)0, (T)1};
      return;
    }

    T invMaxAbsElement = (T)1 / maxAbsElement;
    a00 *= invMaxAbsElement;
    a01 *= invMaxAbsElement;
    a02 *= invMaxAbsElement;
    a11 *= invMaxAbsElement;
    a12 *= invMaxAbsElement;
    a22 *= invMaxAbsElement;

    T norm = a01 * a01 + a02 * a02 + a12 * a12;
    if (norm > (T)0) {
      // Compute the eigenvalues of A.

      // In the PDF mentioned previously, B = (A - q*I)/p, where
      // q = tr(A)/3 with tr(A) the trace of A (sum of the diagonal
      // entries of A) and where p = sqrt(tr((A - q*I)^2)/6).
      T q = (a00 + a11 + a22) / (T)3;

      // The matrix A - q*I is represented by the following, where
      // b00, b11 and b22 are computed after these comments,
      //   +-           -+
      //   | b00 a01 a02 |
      //   | a01 b11 a12 |
      //   | a02 a12 b22 |
      //   +-           -+
      T b00 = a00 - q;
      T b11 = a11 - q;
      T b22 = a22 - q;

      // The is the variable p mentioned in the PDF.
      T p = std::sqrt((b00 * b00 + b11 * b11 + b22 * b22 + norm * (T)2) / (T)6);

      // We need det(B) = det((A - q*I)/p) = det(A - q*I)/p^3.  The
      // value det(A - q*I) is computed using a cofactor expansion
      // by the first row of A - q*I.  The cofactors are c00, c01
      // and c02 and the determinant is b00*c00 - a01*c01 + a02*c02.
      // The det(B) is then computed finally by the division
      // with p^3.
      T c00 = b11 * b22 - a12 * a12;
      T c01 = a01 * b22 - a12 * a02;
      T c02 = a01 * a12 - b11 * a02;
      T det = (b00 * c00 - a01 * c01 + a02 * c02) / (p * p * p);

      // The halfDet value is cos(3*theta) mentioned in the PDF. The
      // acos(z) function requires |z| <= 1, but will fail silently
      // and return NaN if the input is larger than 1 in magnitude.
      // To avoid this problem due to rounding errors, the halfDet
      // value is clamped to [-1,1].
      T halfDet = det * (T)0.5;
      halfDet   = std::min(std::max(halfDet, (T)-1), (T)1);

      // The eigenvalues of B are ordered as
      // beta0 <= beta1 <= beta2.  The number of digits in
      // twoThirdsPi is chosen so that, whether float or double,
      // the floating-point number is the closest to theoretical
      // 2*pi/3.
      T angle             = std::acos(halfDet) / (T)3;
      T const twoThirdsPi = (T)2.09439510239319549;
      T beta2             = std::cos(angle) * (T)2;
      T beta0             = std::cos(angle + twoThirdsPi) * (T)2;
      T beta1             = -(beta0 + beta2);

      // The eigenvalues of A are ordered as
      // alpha0 <= alpha1 <= alpha2.
      eval[0] = q + p * beta0;
      eval[1] = q + p * beta1;
      eval[2] = q + p * beta2;

      // Compute the eigenvectors so that the set
      // {evec[0], evec[1], evec[2]} is right handed and
      // orthonormal.
      if (halfDet >= (T)0) {
        ComputeEigenvector0(a00, a01, a02, a11, a12, a22, eval[2], evec[2]);
        ComputeEigenvector1(a00, a01, a02, a11, a12, a22, evec[2], eval[1], evec[1]);
        evec[0] = Cross(evec[1], evec[2]);
      } else {
        ComputeEigenvector0(a00, a01, a02, a11, a12, a22, eval[0], evec[0]);
        ComputeEigenvector1(a00, a01, a02, a11, a12, a22, evec[0], eval[1], evec[1]);
        evec[2] = Cross(evec[0], evec[1]);
      }
    } else {
      // The matrix is diagonal.
      eval[0] = a00;
      eval[1] = a11;
      eval[2] = a22;
      evec[0] = {(T)1, (T)0, (T)0};
      evec[1] = {(T)0, (T)1, (T)0};
      evec[2] = {(T)0, (T)0, (T)1};
    }

    // The preconditioning scaled the matrix A, which scales the
    // eigenvalues.  Revert the scaling.
    eval[0] *= maxAbsElement;
    eval[1] *= maxAbsElement;
    eval[2] *= maxAbsElement;

    SortEigenstuff<T>()(sortType, true, eval, evec);
  }

private:
  static Vector3D
  Multiply(T s, Vector3D const& U)
  {
    Vector3D product = {s * U[0], s * U[1], s * U[2]};
    return product;
  }

  static Vector3D
  Subtract(Vector3D const& U, Vector3D const& V)
  {
    Vector3D difference = {U[0] - V[0], U[1] - V[1], U[2] - V[2]};
    return difference;
  }

  static Vector3D
  Divide(Vector3D const& U, T s)
  {
    T invS                    = (T)1 / s;
    Vector3D division = {U[0] * invS, U[1] * invS, U[2] * invS};
    return division;
  }

  static T
  Dot(Vector3D const& U, Vector3D const& V)
  {
    T dot = U[0] * V[0] + U[1] * V[1] + U[2] * V[2];
    return dot;
  }

  static Vector3D
  Cross(Vector3D const& U, Vector3D const& V)
  {
    Vector3D cross = {U[1] * V[2] - U[2] * V[1], U[2] * V[0] - U[0] * V[2],
                              U[0] * V[1] - U[1] * V[0]};
    return cross;
  }

  void
  ComputeOrthogonalComplement(Vector3D const& W, Vector3D& U,
                              Vector3D& V) const
  {
    // Robustly compute a right-handed orthonormal set { U, V, W }.
    // The vector W is guaranteed to be unit-length, in which case
    // there is no need to worry about a division by zero when
    // computing invLength.
    T invLength;
    if (std::fabs(W[0]) > std::fabs(W[1])) {
      // The component of maximum absolute value is either W[0]
      // or W[2].
      invLength = (T)1 / std::sqrt(W[0] * W[0] + W[2] * W[2]);
      U         = {-W[2] * invLength, (T)0, +W[0] * invLength};
    } else {
      // The component of maximum absolute value is either W[1]
      // or W[2].
      invLength = (T)1 / std::sqrt(W[1] * W[1] + W[2] * W[2]);
      U         = {(T)0, +W[2] * invLength, -W[1] * invLength};
    }
    V = Cross(W, U);
  }

  void
  ComputeEigenvector0(T a00, T a01, T a02, T a11, T a12, T a22, T eval0,
                      Vector3D& evec0) const
  {
    // Compute a unit-length eigenvector for eigenvalue[i0].  The
    // matrix is rank 2, so two of the rows are linearly independent.
    // For a robust computation of the eigenvector, select the two
    // rows whose cross product has largest length of all pairs of
    // rows.
    Vector3D row0 {{a00 - eval0, a01, a02}};
    Vector3D row1 {{a01, a11 - eval0, a12}};
    Vector3D row2 {{a02, a12, a22 - eval0}};
    Vector3D r0xr1 = row0.cross(row1);
    Vector3D r0xr2 = row0.cross(row2);
    Vector3D r1xr2 = row1.cross(row2);
    double d0      = r0xr1.dot(r0xr1);
    double d1      = r0xr2.dot(r0xr2);
    double d2      = r1xr2.dot(r1xr2);

    T dmax       = d0;
    int32_t imax = 0;
    if (d1 > dmax) {
      dmax = d1;
      imax = 1;
    }
    if (d2 > dmax) {
      imax = 2;
    }

    if (imax == 0) {
      evec0 = Divide(r0xr1, std::sqrt(d0));
    } else if (imax == 1) {
      evec0 = Divide(r0xr2, std::sqrt(d1));
    } else {
      evec0 = Divide(r1xr2, std::sqrt(d2));
    }
  }

  void
  ComputeEigenvector1(T a00, T a01, T a02, T a11, T a12, T a22, Vector3D const& evec0,
                      T eval1, Vector3D& evec1) const
  {
    // Robustly compute a right-handed orthonormal set
    // { U, V, evec0 }.
    Vector3D U, V;
    ComputeOrthogonalComplement(evec0, U, V);

    // Let e be eval1 and let E be a corresponding eigenvector which
    // is a solution to the linear system (A - e*I)*E = 0.  The matrix
    // (A - e*I) is 3x3, not invertible (so infinitely many
    // solutions), and has rank 2 when eval1 and eval are different.
    // It has rank 1 when eval1 and eval2 are equal.  Numerically, it
    // is difficult to compute robustly the rank of a matrix.  Instead,
    // the 3x3 linear system is reduced to a 2x2 system as follows.
    // Define the 3x2 matrix J = [U V] whose columns are the U and V
    // computed previously.  Define the 2x1 vector X = J*E.  The 2x2
    // system is 0 = M * X = (J^T * (A - e*I) * J) * X where J^T is
    // the transpose of J and M = J^T * (A - e*I) * J is a 2x2 matrix.
    // The system may be written as
    //     +-                        -++-  -+       +-  -+
    //     | U^T*A*U - e  U^T*A*V     || x0 | = e * | x0 |
    //     | V^T*A*U      V^T*A*V - e || x1 |       | x1 |
    //     +-                        -++   -+       +-  -+
    // where X has row entries x0 and x1.

    Vector3D AU = {a00 * U[0] + a01 * U[1] + a02 * U[2],
                           a01 * U[0] + a11 * U[1] + a12 * U[2],
                           a02 * U[0] + a12 * U[1] + a22 * U[2]};

    Vector3D AV = {a00 * V[0] + a01 * V[1] + a02 * V[2],
                           a01 * V[0] + a11 * V[1] + a12 * V[2],
                           a02 * V[0] + a12 * V[1] + a22 * V[2]};

    T m00 = U[0] * AU[0] + U[1] * AU[1] + U[2] * AU[2] - eval1;
    T m01 = U[0] * AV[0] + U[1] * AV[1] + U[2] * AV[2];
    T m11 = V[0] * AV[0] + V[1] * AV[1] + V[2] * AV[2] - eval1;

    // For robustness, choose the largest-length row of M to compute
    // the eigenvector.  The 2-tuple of coefficients of U and V in the
    // assignments to eigenvector[1] lies on a circle, and U and V are
    // unit length and perpendicular, so eigenvector[1] is unit length
    // (within numerical tolerance).
    T absM00 = std::fabs(m00);
    T absM01 = std::fabs(m01);
    T absM11 = std::fabs(m11);
    T maxAbsComp;
    if (absM00 >= absM11) {
      maxAbsComp = std::max(absM00, absM01);
      if (maxAbsComp > (T)0) {
        if (absM00 >= absM01) {
          m01 /= m00;
          m00 = (T)1 / std::sqrt((T)1 + m01 * m01);
          m01 *= m00;
        } else {
          m00 /= m01;
          m01 = (T)1 / std::sqrt((T)1 + m00 * m00);
          m00 *= m01;
        }
        evec1 = Subtract(Multiply(m01, U), Multiply(m00, V));
      } else {
        evec1 = U;
      }
    } else {
      maxAbsComp = std::max(absM11, absM01);
      if (maxAbsComp > (T)0) {
        if (absM11 >= absM01) {
          m01 /= m11;
          m11 = (T)1 / std::sqrt((T)1 + m01 * m01);
          m01 *= m11;
        } else {
          m11 /= m01;
          m01 = (T)1 / std::sqrt((T)1 + m11 * m11);
          m11 *= m01;
        }
        evec1 = Subtract(Multiply(m11, U), Multiply(m01, V));
      } else {
        evec1 = U;
      }
    }
  }
};
} // namespace OpenSees
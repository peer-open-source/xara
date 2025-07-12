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
// symmetric real-valued matrix. The iterative algorithm is implemented
// by class SymmmetricEigensolver3x3. The noniterative algorithm is
// implemented by class NISymmetricEigensolver3x3.

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

template <typename T, bool aggressive, int32_t sortType> 
class SymEigIterate3D {
public:
  // The input matrix must be symmetric, so only the unique elements
  // must be specified: a00, a01, a02, a11, a12, and a22.
  //
  // If 'aggressive' is 'true', the iterations occur until a
  // superdiagonal entry is exactly zero.  If 'aggressive' is 'false',
  // the iterations occur until a superdiagonal entry is effectively
  // zero compared to the/ sum of magnitudes of its diagonal neighbors.
  // Generally, the nonaggressive convergence is acceptable.
  //
  // The order of the eigenvalues is specified by sortType:
  // -1 (decreasing), 0 (no sorting) or +1 (increasing).  When sorted,
  // the eigenvectors are ordered accordingly, and
  // {evec[0], evec[1], evec[2]} is guaranteed to/ be a right-handed
  // orthonormal set.  The return value is the number of iterations
  // used by the algorithm.

  int32_t
  operator()(T const& a00, T const& a01, T const& a02, T const& a11, T const& a12, T const& a22,
             Vector3D& eval, std::array<Vector3D, 3>& evec) const
  {
    // Compute the Householder reflection H0 and B = H0*A*H0, where
    // b02 = 0. H0 = {{c,s,0},{s,-c,0},{0,0,1}} with each inner
    // triple a row of H0.
    T const zero    = static_cast<T>(0);
    T const one     = static_cast<T>(1);
    T const half    = static_cast<T>(0.5);
    bool isRotation = false;
    T c = zero, s = zero;
    GetCosSin(a12, -a02, c, s);
    T term0 = c * a00 + s * a01;
    T term1 = c * a01 + s * a11;
    T term2 = s * a00 - c * a01;
    T term3 = s * a01 - c * a11;
    T b00   = c * term0 + s * term1;
    T b01   = s * term0 - c * term1;
    //T b02 = c * a02 + s * a12;  // 0
    T b11 = s * term2 - c * term3;
    T b12 = s * a02 - c * a12;
    T b22 = a22;

    // Maintain Q as the product of the reflections. Initially,
    // Q = H0. Updates by Givens reflections G are Q <- Q * G. The
    // columns of the final Q are the estimates for the eigenvectors.
    std::array<Vector3D, 3> Q = {{{c,  s, zero},
                                  {s, -c, zero}, 
                                  {zero, zero, one}}};

    // The smallest subnormal number is 2^{-alpha}. The value alpha is
    // 149 for 'float' and 1074 for 'double'.
    int32_t constexpr alpha = std::numeric_limits<T>::digits - std::numeric_limits<T>::min_exponent;
    int32_t i = 0, imax = 0, power = 0;
    T c2 = zero, 
      s2 = zero;

    if (std::fabs(b12) <= std::fabs(b01)) {
      // It is known that |currentB12| < 2^{-i/2} * |initialB12|.
      // Compute imax so that 0 is the closest floating-point number
      // to 2^{-imax/2} * |initialB12|.
      (void)std::frexp(b12, &power);
      imax = 2 * (power + alpha + 1);

      for (i = 0; i < imax; ++i) {
        // Compute the Givens reflection
        // G = {{c,0,-s},{s,0,c},{0,1,0}} where each inner triple
        // is a row of G.
        GetCosSin(half * (b00 - b11), b01, c2, s2);
        s = std::sqrt(half * (one - c2));
        c = half * s2 / s;

        // Update Q <- Q * G.
        for (size_t r = 0; r < 3; ++r) {
          term0   = c * Q[r][0] + s * Q[r][1];
          term1   = Q[r][2];
          term2   = c * Q[r][1] - s * Q[r][0];
          Q[r][0] = term0;
          Q[r][1] = term1;
          Q[r][2] = term2;
        }
        isRotation = !isRotation;

        // Update B <- Q^T * B * Q, ensuring that b02 is zero and
        // |b12| has strictly decreased.
        term0 = c * b00 + s * b01;
        term1 = c * b01 + s * b11;
        term2 = s * b00 - c * b01;
        term3 = s * b01 - c * b11;
        //b02 = s * c * (b11 - b00) + (c * c - s * s) * b01; // 0
        b00 = c * term0 + s * term1;
        b01 = s * b12;
        b11 = b22;
        b12 = c * b12;
        b22 = s * term2 - c * term3;

        if (Converged(aggressive, b00, b11, b01)) {
          // Compute the Householder reflection
          // H1 = {{c,s,0},{s,-c,0},{0,0,1}} where each inner
          // triple is a row of H1.
          GetCosSin(half * (b00 - b11), b01, c2, s2);
          s = std::sqrt(half * (one - c2));
          c = half * s2 / s;

          // Update Q <- Q * H1.
          for (size_t r = 0; r < 3; ++r) {
            term0   = c * Q[r][0] + s * Q[r][1];
            term1   = s * Q[r][0] - c * Q[r][1];
            Q[r][0] = term0;
            Q[r][1] = term1;
          }
          isRotation = !isRotation;

          // Compute the diagonal estimate D = Q^T * B * Q.
          term0 = c * b00 + s * b01;
          term1 = c * b01 + s * b11;
          term2 = s * b00 - c * b01;
          term3 = s * b01 - c * b11;
          b00   = c * term0 + s * term1;
          b11   = s * term2 - c * term3;
          break;
        }
      }
    }
    else {
      // It is known that |currentB01| < 2^{-i/2} * |initialB01|.
      // Compute imax so that 0 is the closest floating-point number
      // to 2^{-imax/2} * |initialB01|.
      (void)std::frexp(b01, &power);
      imax = 2 * (power + alpha + 1);

      for (i = 0; i < imax; ++i) {
        // Compute the Givens reflection
        // G = {{0,1,0},{c,0,-s},{s,0,c}} where each inner triple
        // is a row of G.
        GetCosSin(half * (b11 - b22), b12, c2, s2);
        s = std::sqrt(half * (one - c2));
        c = half * s2 / s;

        // Update Q <- Q * G.
        for (size_t r = 0; r < 3; ++r) {
          term0   = c * Q[r][1] + s * Q[r][2];
          term1   = Q[r][0];
          term2   = c * Q[r][2] - s * Q[r][1];
          Q[r][0] = term0;
          Q[r][1] = term1;
          Q[r][2] = term2;
        }
        isRotation = !isRotation;

        // Update B <- Q^T * B * Q, ensuring that b02 is zero and
        // |b01| has strictly decreased.
        term0 = c * b11 + s * b12;
        term1 = c * b12 + s * b22;
        term2 = s * b11 - c * b12;
        term3 = s * b12 - c * b22;
        //b02 = s * c * (b22 - b11) + (c * c - s * s) * b12;  // 0
        b22 = s * term2 - c * term3;
        b12 = -s * b01;
        b11 = b00;
        b01 = c * b01;
        b00 = c * term0 + s * term1;

        if (Converged(aggressive, b11, b22, b12)) {
          // Compute the Householder reflection
          // H1 = {{1,0,0},{0,c,s},{0,s,-c}} where each inner
          // triple is a row of H1.
          GetCosSin(half * (b11 - b22), b12, c2, s2);
          s = std::sqrt(half * (one - c2));
          c = half * s2 / s;

          // Update Q <- Q * H1.
          for (size_t r = 0; r < 3; ++r) {
            term0   = c * Q[r][1] + s * Q[r][2];
            term1   = s * Q[r][1] - c * Q[r][2];
            Q[r][1] = term0;
            Q[r][2] = term1;
          }
          isRotation = !isRotation;

          // Compute the diagonal estimate D = Q^T * B * Q.
          term0 = c * b11 + s * b12;
          term1 = c * b12 + s * b22;
          term2 = s * b11 - c * b12;
          term3 = s * b12 - c * b22;
          b11   = c * term0 + s * term1;
          b22   = s * term2 - c * term3;
          break;
        }
      }
    }

    eval = {b00, b11, b22};
    for (size_t row = 0; row < 3; ++row) {
      for (size_t col = 0; col < 3; ++col) {
        evec[row][col] = Q[col][row];
      }
    }

    SortEigenstuff<T>()(sortType, isRotation, eval, evec);
    return i;
  }

private:
  // Normalize (u,v) to (c,s) with c <= 0 when (u,v) is not (0,0).
  // If (u,v) = (0,0), the function returns (c,s) = (-1,0). When used
  // to generate a Householder reflection, it does not matter whether
  // (c,s) or (-c,-s) is returned. When generating a Givens reflection,
  // c = cos(2*theta) and s = sin(2*theta). Having a negative cosine
  // for the double-angle term ensures that the single-angle terms
  // c = cos(theta) and s = sin(theta) satisfy |c| < 1/sqrt(2) < |s|.
  static void
  GetCosSin(T const& u, T const& v, T& c, T& s)
  {
    T const zero = static_cast<T>(0);
    T length     = std::sqrt(u * u + v * v);
    if (length > zero) {
      c = u / length;
      s = v / length;
      if (c > zero) {
        c = -c;
        s = -s;
      }
    } else {
      c = static_cast<T>(-1);
      s = zero;
    }
  }

  static bool
  Converged(bool aggressive, T const& diagonal0, T const& diagonal1, T const& superdiagonal)
  {
    if (aggressive) {
      // Test whether the superdiagonal term is zero.
      return superdiagonal == static_cast<T>(0);
    } else {
      // Test whether the superdiagonal term is effectively zero
      // compared to its diagonal neighbors.
      T sum = std::fabs(diagonal0) + std::fabs(diagonal1);
      return sum + std::fabs(superdiagonal) == sum;
    }
  }
};
} // namespace OpenSees
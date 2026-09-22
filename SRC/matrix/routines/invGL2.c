//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2026, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// adapted from https://caps.gsfc.nasa.gov/simpson/software/m22inv_f90.txt
//
// David Simpson
// Claudio Perez
//
#include <math.h>

int
cmx_inv2(double *a, double * ainv, int *ok_flag__)
{
  // ****************************************************************************************
  //  m22inv  -  compute the inverse of a 2x2 matrix.
  //
  //  a      : (input)  2x2 matrix to be inverted 
  //  ainv   : (output) 2x2 inverse of matrix a 
  //
  //  ok_flag: (output) 0 if the input matrix could be inverted, 
  //           and -1 if the input matrix is singular. 
  //
  // ****************************************************************************************/

  // Parameter adjustments
  ainv -= 3;
  a    -= 3;

  // Function Body
  const double eps = 1e-10;
  const double det = a[3] * a[6] - a[5] * a[4];
  if (fabs(det) <= eps) {
    *ok_flag__ = -1;
  }

  double cofactor[4];

  cofactor[0] =  a[6];
  cofactor[2] = -a[4];
  cofactor[1] = -a[5];
  cofactor[3] =  a[3];

  for (int i__ = 1; i__ <= 2; ++i__)
    for (int j = 1; j <= 2; ++j)
      ainv[j + (i__ << 1)] = cofactor[i__ + (j << 1) - 3] / det;

  *ok_flag__ = 0;
  return 0;
}

// Solve a 2x2 system A x = b for x, with A stored column-major.
//
//   A = | a[0]  a[2] |     b = | b[0] |     x = | x[0] |
//       | a[1]  a[3] |         | b[1] |         | x[1] |
//
// *det is set to det(A) and is written even when A is singular, so the
// caller can assess conditioning. 
// x is computed only when A is nonsingular.
// x is permitted to alias a or b. 
// det may be NULL.
//
// Returns 0 on success, 1 if A is singular (det == 0).
//
int
cmx_solve2(const double *a, const double *b, double *x, double *det)
{
  const double d = a[0] * a[3] - a[2] * a[1];

  if (det)
    *det = d;

  if (d == 0.0)
    return 1;

  const double invd = 1.0 / d;

  // Numerators are formed before any store so x may alias a or b.
  const double x0 = b[0] * a[3] - a[2] * b[1];
  const double x1 = a[0] * b[1] - b[0] * a[1];

  x[0] = x0 * invd;
  x[1] = x1 * invd;

  return 0;
}

// Solve a 2x2 system A x = b for x, with A stored column-major.
//
//   A = | a[0]  a[2] |     b = | b[0] |     x = | x[0] |
//       | a[1]  a[3] |         | b[1] |         | x[1] |
//
// a, b, x, det must point to four DISTINCT, NON-OVERLAPPING objects.
// Unlike cmx_solve2, this variant permits no aliasing and det must NOT
// be NULL. Violating either is undefined behavior.
//
// *det receives det(A) unconditionally so the caller can judge
// conditioning. x is computed only when A is nonsingular.
//
// Returns 0 on success, 1 if A is singular (det == 0).
//
int
cmx_rsolve2(const double * restrict a,
            const double * restrict b,
            double       * restrict x)
{
  const double d = a[0] * a[3] - a[2] * a[1];

  if (d == 0.0)
    return 1;

  const double invd = 1.0 / d;

  x[0] = (b[0] * a[3] - a[2] * b[1]) * invd;
  x[1] = (a[0] * b[1] - b[0] * a[1]) * invd;

  return 0;
}
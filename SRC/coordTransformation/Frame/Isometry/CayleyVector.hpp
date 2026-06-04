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
// Written: Claudio M. Perez
// Created: 06/2026
//
#pragma once

#include <Vector3D.h>
#include <Matrix3D.h>
#include <GroupSO3.h>

namespace OpenSees {

class CayleyVector
{
public:
  explicit CayleyVector(const Vector3D& v)
    : v(v)
  {}

  const Vector3D&
  vector() const noexcept
  {
    return v;
  }

  Matrix3D
  map() const
  {
    const double c = scale();

    Matrix3D R{};
    R.addDiagonal(1.0);
    R.addSpin(v, c);
    R.addSpinSquare(v, 0.5*c);

    return R;
  }

  Matrix3D
  tan() const
  {
    //
    // J_C(v), defined by
    //
    //   d Cay(v) Cay(v)^T = [J_C(v) dv]_x
    //
    const double c = scale();

    Matrix3D J{};
    J.addDiagonal(c);
    J.addSpin(v, 0.5*c);

    return J;
  }

  Matrix3D
  taninv() const
  {
    //
    // K_C(v) = J_C(v)^{-1}
    //
    Matrix3D K{};
    K.addDiagonal(1.0);
    K.addSpin(v, -0.5);
    K.addTensorProduct(v, v, 0.25);

    return K;
  }

  Matrix3D
  dtan(const Vector3D& y) const
  {
    //
    // D J_C(v)[y]
    //
    // J_C(v) = c(v) * (I + 1/2 [v]_x)
    // c(v)   = 1 / (1 + v.v/4)
    //
    // dc[y] = -c^2 * (v.y)/2
    //
    const double c  = scale();
    const double dc = -0.5*c*c*v.dot(y);

    Matrix3D dJ{};
    dJ.addDiagonal(dc);
    dJ.addSpin(v, 0.5*dc);
    dJ.addSpin(y, 0.5*c);

    return dJ;
  }

  Matrix3D
  dtaninv(const Vector3D& y) const
  {
    //
    // D K_C(v)[y]
    //
    // K_C(v) = I - 1/2 [v]_x + 1/4 v v^T
    //
    Matrix3D dK{};
    dK.addSpin(y, -0.5);
    dK.addTensorProduct(y, v, 0.25);
    dK.addTensorProduct(v, y, 0.25);

    return dK;
  }

  Matrix3D
  dmap(const Vector3D& y) const
  {
    //
    // D Cay(v)[y] = [J_C(v)y]_x Cay(v)
    //
    return Hat(tan()*y)*map();
  }

  static CayleyVector
  fromRotation(const Matrix3D& R)
  {
    return CayleyVector{CayleyFromVersor(Versor::from_matrix(R))};
  }

private:
  double
  scale() const noexcept
  {
    return 1.0/(1.0 + 0.25*v.dot(v));
  }

private:
  Vector3D v;
};

} // namespace OpenSees
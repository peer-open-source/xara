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
/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Description: J2 plasticity material where stress
// components 22, 33, 13, and 23 are constrained to 0.
//
// Written: CMP
// Created: Feb 2026
//
#include <J2BeamThread3d.h>
#include <Channel.h>
#include <string.h>
#include <cmath>
#include <float.h>
#include <MatrixND.h>
#include <VectorND.h>
#include <string.h>

using namespace OpenSees;

#if 1
namespace {
struct J2Parameters {
  double Hiso, Hkin, E, nu, sigmaY;
};

static int
ReturnMap(const J2Parameters& params, 
          Vector3D& xsi,
          const Vector3D& Tepsilon,
          const double epsPn[3],
          double epsPn1[3], 
          const double& alphan,
          double& alphan1,
          double& dg
        )

{
  const double &E = params.E,
               &nu = params.nu,
               &sigmaY = params.sigmaY,
               &Hiso = params.Hiso,
               &Hkin = params.Hkin;

  double twoG = E / (1.0 + nu);
  double G    = 0.5 * twoG;
  Vector3D sig;
  sig[0] = E * (Tepsilon(0) - epsPn[0]);
  sig[1] = G * (Tepsilon(1) - epsPn[1]);
  sig[2] = G * (Tepsilon(2) - epsPn[2]);

  static constexpr double one3 = 1.0/3.0;
  static constexpr double two3 = 2.0 * one3;
  static const double root23   = std::sqrt(2.0/3.0);

  double two3Hkin = two3 * Hkin;

  // relative stress
  xsi[0] = sig[0] - Hkin * epsPn[0];
  xsi[1] = sig[1] - Hkin * epsPn[1]/3.0;
  xsi[2] = sig[2] - Hkin * epsPn[2]/3.0;

  double q = std::sqrt( 2.0/3.0 * xsi[0]*xsi[0]
                      + 2.0*xsi[1] * xsi[1]
                      + 2.0*xsi[2] * xsi[2]);

  double F = q - root23*(sigmaY + Hiso*alphan);

  if (F < -100 * DBL_EPSILON) {
    epsPn1[0] = epsPn[0];
    epsPn1[1] = epsPn[1];
    epsPn1[2] = epsPn[2];
    alphan1   = alphan;
    dg        = 0.0;
    return 0;
  }
  else {
    // Solve for dg
    dg = 0.0;

    VectorND<4> R{0.0, 0.0, 0.0, F};
    VectorND<4> x{xsi[0], xsi[1], xsi[2], dg};

    MatrixND<4, 4> J{};
    VectorND<4> dx{};

    int iter                    = 0;
    int maxIter                 = 25;
    constexpr static double tol = 1.0e-8;

    do {
      J(0, 0) = 1.0 + dg*two3*(E + Hkin);
      J(0, 1) = 0.0;
      J(0, 2) = 0.0;
      J(1, 0) = 0.0;
      J(1, 1) = 1.0 + dg * 2 * (G + Hkin/3.0);
      J(1, 2) = 0.0;
      J(2, 0) = 0.0;
      J(2, 1) = 0.0;
      J(2, 2) = 1.0 + dg * 2 * (G + Hkin/3.0);

      J(0, 3) = two3 * (E + Hkin) * x(0);
      J(1, 3) = 2.0 * (G + Hkin/3.0) * x(1);
      J(2, 3) = 2.0 * (G + Hkin/3.0) * x(2);

      //J(2,0) = x(0)*two3/q; J(2,1) = x(1)*2.0/q;
      J(3, 0) = (1.0 - 2.0/3.0*Hiso * dg) * x[0] * two3 / q;
      J(3, 1) = (1.0 - 2.0/3.0*Hiso * dg) * x[1] * 2.0 / q;
      J(3, 2) = (1.0 - 2.0/3.0*Hiso * dg) * x[2] * 2.0 / q;

      //J(2,2) = -root23*Hiso;
      J(3, 3) = -two3 * Hiso * q;

      J.solve(R, dx);
      x.addVector(1.0, dx, -1.0);

      dg = x(3);

      q  = std::sqrt(two3 * x(0)*x(0) + 2.0*x(1)*x(1) + 2.0*x(2)*x(2));

      R(0) = x(0) - xsi[0] + dg * two3 * (E + Hkin) * x(0);
      R(1) = x(1) - xsi[1] + dg * (twoG + two3Hkin) * x(1);
      R(2) = x(2) - xsi[2] + dg * (twoG + two3Hkin) * x(2);
      R(3) = q - root23 * (sigmaY + Hiso * (alphan + dg * root23 * q));

    } while ((R.norm() > sigmaY * tol) && (iter++ < maxIter));

    if (R.norm() > sigmaY * tol) {
      opserr << "J2BeamThread3d::getTangent -- maxIter reached, "
             << R.norm() << " > "
             << sigmaY*tol << "\n";
      return -1;
    }

    //
    alphan1 = alphan + dg * root23 * q;
    epsPn1[0] = epsPn[0] + dg * two3 * x(0);
    epsPn1[1] = epsPn[1] + dg * 2.0  * x(1);
    epsPn1[2] = epsPn[2] + dg * 2.0  * x(2);
    return 1;
  }

  return 0;
} 
} // end namespace
#endif

J2BeamThread3d::J2BeamThread3d(int tag, 
          double e, double g, double sy, double hi, double hk)
 : NDMaterial(tag, ND_TAG_J2BeamThread3d),
   E(e),
   nu(g),
   sigmaY(sy),
   Hiso(hi),
   Hkin(hk),
   parameterID(0),
   SHVs(0),
   Tepsilon{},
   sigma{},
   D_wrapper(D),
   s_wrapper(sigma),
   e_wrapper(Tepsilon),
   alphan(0.0),
   alphan1(0.0),
   dg_n1(0.0)
{
  epsPn[0] = 0.0;
  epsPn[1] = 0.0;
  epsPn[2] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
  epsPn1[2] = 0.0;
}

J2BeamThread3d::J2BeamThread3d()
 : NDMaterial(0, ND_TAG_J2BeamThread3d),
   E(0.0),
   nu(0.0),
   sigmaY(0.0),
   Hkin(0.0),
   parameterID(0),
   SHVs(0),
   Tepsilon{},
   sigma{},
   s_wrapper(sigma),
   e_wrapper(Tepsilon),
   alphan(0.0),
   alphan1(0.0),
   dg_n1(0.0)
{
  epsPn[0] = 0.0;
  epsPn[1] = 0.0;
  epsPn[2] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
  epsPn1[2] = 0.0;
}

J2BeamThread3d::~J2BeamThread3d()
{
  if (SHVs != nullptr)
    delete SHVs;
}

int
J2BeamThread3d::setTrialStrain(const Vector& strain)
{
  Tepsilon = strain;
  return 0;
}


int
J2BeamThread3d::setTrialStrainIncr(const Vector& strain)
{
  assert(false);
  return 0;
}


const Matrix&
J2BeamThread3d::getTangent()
{
  double twoG = E / (1.0 + nu);
  double G    = 0.5 * twoG;

  static constexpr double one3 = 1.0 / 3.0;
  static constexpr double two3 = 2.0 * one3;
  static const double root23   = std::sqrt(2.0/3.0);

  double two3Hkin = two3 * Hkin;
#if 1
  Vector3D sig;
  sig[0] = E * (Tepsilon(0) - epsPn[0]);
  sig[1] = G * (Tepsilon(1) - epsPn[1]);
  sig[2] = G * (Tepsilon(2) - epsPn[2]);
  Vector3D xsi;
  //xsi[0] = sig[0] - two3*Hkin*1.5*epsPn[0];
  //xsi[1] = sig[1] - two3*Hkin*0.5*epsPn[1];
  xsi[0] = sig[0] - Hkin * epsPn[0];
  xsi[1] = sig[1] - one3 * Hkin * epsPn[1];
  xsi[2] = sig[2] - one3 * Hkin * epsPn[2];

  double q = std::sqrt( 2.0/3.0 * xsi[0]*xsi[0]
                      + 2.0*xsi[1]*xsi[1]
                      + 2.0*xsi[2]*xsi[2]);
  double F = q - root23 * (sigmaY + Hiso*alphan);

  if (F < -100 * DBL_EPSILON) {
    D(0, 0) = E;
    D(1, 1) = G;
    D(2, 2) = G;
    D(0, 1) = D(1, 0) = 0.0;
    D(0, 2) = D(2, 0) = 0.0;
    D(1, 2) = D(2, 1) = 0.0;

    epsPn1[0] = epsPn[0];
    epsPn1[1] = epsPn[1];
    epsPn1[2] = epsPn[2];
    alphan1   = alphan;
  }
  else {
    // Solve for dg
    double dg = 0.0;

    VectorND<4> R{0.0, 0.0, 0.0, F};
    VectorND<4> x{xsi[0], xsi[1], xsi[2], dg};

    MatrixND<4, 4> J{};
    VectorND<4> dx{};

    int iter                    = 0;
    int maxIter                 = 25;
    constexpr static double tol = 1.0e-8;

    do {
      J(0, 0) = 1.0 + dg * two3 * (E + Hkin);
      J(0, 1) = 0.0;
      J(0, 2) = 0.0;
      J(1, 0) = 0.0;
      J(1, 1) = 1.0 + dg * 2.0 * (G + Hkin/3.0);
      J(1, 2) = 0.0;
      J(2, 0) = 0.0;
      J(2, 1) = 0.0;
      J(2, 2) = 1.0 + dg * 2.0 * (G + Hkin/3.0);

      J(0, 3) = two3 * (E + Hkin) * x(0);
      J(1, 3) = 2.0 * (G + Hkin/3.0) * x(1);
      J(2, 3) = 2.0 * (G + Hkin/3.0) * x(2);

      //J(2,0) = x(0)*two3/q; J(2,1) = x(1)*2.0/q;
      J(3, 0) = (1.0 - 2.0 / 3.0 * Hiso * dg) * x[0] * two3 / q;
      J(3, 1) = (1.0 - 2.0 / 3.0 * Hiso * dg) * x[1] * 2.0 / q;
      J(3, 2) = (1.0 - 2.0 / 3.0 * Hiso * dg) * x[2] * 2.0 / q;

      //J(2,2) = -root23*Hiso;
      J(3, 3) = -two3 * Hiso * q;

      J.solve(R, dx);
      x.addVector(1.0, dx, -1.0);

      dg    = x(3);
      dg_n1 = dg;

      q = std::sqrt(two3*x(0)*x(0) + 2.0*x(1)*x(1) + 2.0*x(2)*x(2));

      R(0) = x(0) - xsi[0] + dg * two3 * (E + Hkin) * x(0);
      R(1) = x(1) - xsi[1] + dg * (twoG + two3Hkin) * x(1);
      R(2) = x(2) - xsi[2] + dg * (twoG + two3Hkin) * x(2);
      R(3) = q - root23 * (sigmaY + Hiso * (alphan + dg * root23 * q));
    } while ((R.norm() > sigmaY * tol) && (iter++ < maxIter));

    if (R.norm() > sigmaY * tol) {
      opserr << "J2BeamThread3d::getTangent -- maxIter reached, "
             << R.norm() << " > "
             << sigmaY * tol << "\n";
    }

    alphan1   = alphan   + dg*root23*q;

    epsPn1[0] = epsPn[0] + dg * two3 * x(0);
    epsPn1[1] = epsPn[1] + dg * 2.0 * x(1);
    epsPn1[2] = epsPn[2] + dg * 2.0 * x(2);

    //J(2,0) = (1.0-two3*Hiso*dg)*x(0)*two3/q; J(2,1) = (1.0-two3*Hiso*dg)*x(1)*2.0/q;
    //J(2,2) = -two3*Hiso*q;
    //static Matrix invJ(3,3);
    //J.Invert(invJ);

    J(0, 0) = 1.0 + dg * two3 * E / (1.0 + dg * two3Hkin);
    J(0, 1) = 0.0;
    J(0, 2) = 0.0;
    J(1, 0) = 0.0;
    J(1, 1) = 1.0 + dg * twoG / (1.0 + dg * two3Hkin);
    J(1, 2) = 0.0;
    J(2, 0) = 0.0;
    J(2, 1) = 0.0;
    J(2, 2) = 1.0 + dg * twoG / (1.0 + dg * two3Hkin);

    J(0, 3) = (two3 * E - dg * two3 * E / (1.0 + dg * two3Hkin) * two3Hkin) * x(0);
    J(1, 3) = (twoG - dg * twoG / (1.0 + dg*two3Hkin) * two3Hkin) * x(1);
    J(2, 3) = (twoG - dg * twoG / (1.0 + dg*two3Hkin) * two3Hkin) * x(2);

    //J(2,0) = x(0)/q*two3/(1.0+dg*two3Hkin);
    //J(2,1) = x(1)/q* 2.0/(1.0+dg*two3Hkin);
    J(3, 0) = (1.0 - two3 * Hiso * dg) * x(0) / q * two3 / (1.0 + dg * two3Hkin);
    J(3, 1) = (1.0 - two3 * Hiso * dg) * x(1) / q * 2.0 / (1.0 + dg * two3Hkin);
    J(3, 2) = (1.0 - two3 * Hiso * dg) * x(2) / q * 2.0 / (1.0 + dg * two3Hkin);

    //J(2,2) = -(x(0)/q*two3/(1.0+dg*two3Hkin)*two3Hkin*x(0))
    //         -(x(1)/q* 2.0/(1.0+dg*two3Hkin)*two3Hkin*x(1));
    //J(2,2) = -q*two3Hkin/(1.0+dg*two3Hkin) - root23*Hiso;
    J(3, 3) = -q * two3Hkin / (1.0 + dg * two3Hkin) - two3 * Hiso * q;

    MatrixND<4, 4> invJ;
    J.invert(invJ);

    D(0, 0) = invJ(0, 0) * E;
    D(1, 0) = invJ(1, 0) * E;
    D(2, 0) = invJ(2, 0) * E;
    D(0, 1) = invJ(0, 1) * G;
    D(1, 1) = invJ(1, 1) * G;
    D(2, 1) = invJ(2, 1) * G;
    D(0, 2) = invJ(0, 2) * G;
    D(1, 2) = invJ(1, 2) * G;
    D(2, 2) = invJ(2, 2) * G;
  }

#else
    MatrixND<4, 4> J{};

    J(0, 0) = 1.0 + dg * two3 * E / (1.0 + dg * two3Hkin);
    J(0, 1) = 0.0;
    J(0, 2) = 0.0;
    J(1, 0) = 0.0;
    J(1, 1) = 1.0 + dg * twoG / (1.0 + dg * two3Hkin);
    J(1, 2) = 0.0;
    J(2, 0) = 0.0;
    J(2, 1) = 0.0;
    J(2, 2) = 1.0 + dg * twoG / (1.0 + dg * two3Hkin);

    J(0, 3) = (two3 * E - dg * two3 * E / (1.0 + dg * two3Hkin) * two3Hkin) * x(0);
    J(1, 3) = (twoG - dg * twoG / (1.0 + dg*two3Hkin) * two3Hkin) * x(1);
    J(2, 3) = (twoG - dg * twoG / (1.0 + dg*two3Hkin) * two3Hkin) * x(2);

    //J(2,0) = x(0)/q*two3/(1.0+dg*two3Hkin);
    //J(2,1) = x(1)/q* 2.0/(1.0+dg*two3Hkin);
    J(3, 0) = (1.0 - two3 * Hiso * dg) * x(0) / q * two3 / (1.0 + dg * two3Hkin);
    J(3, 1) = (1.0 - two3 * Hiso * dg) * x(1) / q * 2.0 / (1.0 + dg * two3Hkin);
    J(3, 2) = (1.0 - two3 * Hiso * dg) * x(2) / q * 2.0 / (1.0 + dg * two3Hkin);

    //J(2,2) = -(x(0)/q*two3/(1.0+dg*two3Hkin)*two3Hkin*x(0))
    //         -(x(1)/q* 2.0/(1.0+dg*two3Hkin)*two3Hkin*x(1));
    //J(2,2) = -q*two3Hkin/(1.0+dg*two3Hkin) - root23*Hiso;
    J(3, 3) = -q * two3Hkin / (1.0 + dg * two3Hkin) - two3 * Hiso * q;

    MatrixND<4, 4> invJ;
    J.invert(invJ);

    D(0, 0) = invJ(0, 0) * E;
    D(1, 0) = invJ(1, 0) * E;
    D(2, 0) = invJ(2, 0) * E;
    D(0, 1) = invJ(0, 1) * G;
    D(1, 1) = invJ(1, 1) * G;
    D(2, 1) = invJ(2, 1) * G;
    D(0, 2) = invJ(0, 2) * G;
    D(1, 2) = invJ(1, 2) * G;
    D(2, 2) = invJ(2, 2) * G;
#endif
  return D_wrapper;
}

const Matrix&
J2BeamThread3d::getInitialTangent()
{
  double G = 0.5*E/(1.0 + nu);

  D(0, 0) = E;
  D(1, 1) = G;
  D(2, 2) = G;
  D(0, 1) = D(1, 0) = 0.0;
  D(0, 2) = D(2, 0) = 0.0;
  D(1, 2) = D(2, 1) = 0.0;

  return D_wrapper;
}

const Vector&
J2BeamThread3d::getStress()
{

#if 1
  double G = 0.5 * E/(1.0 + nu);
  sigma(0) = E * (Tepsilon(0) - epsPn[0]);
  sigma(1) = G * (Tepsilon(1) - epsPn[1]);
  sigma(2) = G * (Tepsilon(2) - epsPn[2]);

  static constexpr double one3   = 1.0/3.0;
  static constexpr double two3   = 2.0/3.0;
  static const double root23 = std::sqrt(two3);

  double xsi[3];
  //xsi[0] = sigma(0) - two3*Hkin*1.5*epsPn[0];
  //xsi[1] = sigma(1) - two3*Hkin*0.5*epsPn[1];
  xsi[0] = sigma(0) - Hkin * epsPn[0];
  xsi[1] = sigma(1) - one3 * Hkin * epsPn[1];
  xsi[2] = sigma(2) - one3 * Hkin * epsPn[2];

  double q = std::sqrt(two3 * xsi[0] * xsi[0] + 2.0*xsi[1]*xsi[1] + 2.0*xsi[2]*xsi[2]);
  double F = q - root23 * (sigmaY + Hiso*alphan);

  if (F < -100 * DBL_EPSILON) {
    epsPn1[0] = epsPn[0];
    epsPn1[1] = epsPn[1];
    epsPn1[2] = epsPn[2];
    alphan1   = alphan;
  }
  else {
    // Solve for dg
    double dg = 0.0;

    VectorND<4> R;
    R(0) = 0.0;
    R(1) = 0.0;
    R(2) = 0.0;
    R(3) = F;
    VectorND<4> x;
    x(0) = xsi[0];
    x(1) = xsi[1];
    x(2) = xsi[2];
    x(3) = dg;

    static MatrixND<4, 4> J;
    static VectorND<4> dx;

    int iter    = 0;
    int maxIter = 25;
    while ((iter < maxIter) && (R.norm() > sigmaY*1.0e-14)) {
      iter++;

      J(0, 0) = 1.0 + dg * two3 * (E + Hkin);
      J(0, 1) = 0.0;
      J(0, 2) = 0.0;
      J(1, 0) = 0.0;
      J(1, 1) = 1.0 + dg * (2.0 * G + two3 * Hkin);
      J(1, 2) = 0.0;
      J(2, 0) = 0.0;
      J(2, 1) = 0.0;
      J(2, 2) = 1.0 + dg * (2.0 * G + two3 * Hkin);

      J(0, 3) = two3 * (E + Hkin) * x(0);
      J(1, 3) = (2.0 * G + two3 * Hkin) * x(1);
      J(2, 3) = (2.0 * G + two3 * Hkin) * x(2);

      //J(2,0) = x(0)*two3/q; J(2,1) = x(1)*2.0/q;
      J(3, 0) = (1.0 - two3 * Hiso * dg) * x(0) * two3 / q;
      J(3, 1) = (1.0 - two3 * Hiso * dg) * x(1) * 2.0 / q;
      J(3, 2) = (1.0 - two3 * Hiso * dg) * x(2) * 2.0 / q;

      //J(2,2) = -root23*Hiso;
      J(3, 3) = -two3 * Hiso * q;

      J.solve(R, dx);
      x = x - dx;

      dg    = x(3);
      dg_n1 = dg;

      q = std::sqrt(two3*x(0)*x(0) + 2.0*x(1)*x(1) + 2.0*x(2)*x(2));

      R(0) = x(0) - xsi[0] + dg * two3 * (E + Hkin) * x(0);
      R(1) = x(1) - xsi[1] + dg * (2.0 * G + two3 * Hkin) * x(1);
      R(2) = x(2) - xsi[2] + dg * (2.0 * G + two3 * Hkin) * x(2);
      R(3) = q - root23 * (sigmaY + Hiso * (alphan + dg * root23 * q));
    }

    if (iter == maxIter) {
      opserr << "J2BeamThread3d::getStress -- maxIter reached " << R.norm() << "\n";
    }

    alphan1 = alphan + dg * root23 * q;

    epsPn1[0] = epsPn[0] + dg * two3 * x(0);
    epsPn1[1] = epsPn[1] + dg * 2.0 * x(1);
    epsPn1[2] = epsPn[2] + dg * 2.0 * x(2);

    //sigma(0) = x(0) + two3*Hkin*1.5*epsPn1[0];
    //sigma(1) = x(1) + two3*Hkin*0.5*epsPn1[1];
    sigma(0) = x(0) +        Hkin * epsPn1[0];
    sigma(1) = x(1) + one3 * Hkin * epsPn1[1];
    sigma(2) = x(2) + one3 * Hkin * epsPn1[2];
  }
#else
  struct J2Parameters params{Hiso, Hkin, E, nu, sigmaY};
  Vector3D sig;
  int ret = ReturnMap(params, sig, Tepsilon, epsPn, epsPn1, alphan, alphan1, dg_n1);
  if (ret < 0) {
    opserr << "J2BeamThread3d::getStress -- ReturnMap failed\n";
  }
  sigma(0) = sig[0];
  sigma(1) = sig[1];
  sigma(2) = sig[2];
#endif
  return s_wrapper;
}

const Vector&
J2BeamThread3d::getStrain()
{
  return e_wrapper;
}

int
J2BeamThread3d::commitState()
{
  epsPn[0] = epsPn1[0];
  epsPn[1] = epsPn1[1];
  epsPn[2] = epsPn1[2];
  alphan   = alphan1;

  return 0;
}


int
J2BeamThread3d::revertToLastCommit()
{
  epsPn1[0] = epsPn[0];
  epsPn1[1] = epsPn[1];
  epsPn1[2] = epsPn[2];
  alphan1   = alphan;
  return 0;
}


int
J2BeamThread3d::revertToStart()
{
  Tepsilon.zero();

  epsPn[0] = 0.0;
  epsPn[1] = 0.0;
  epsPn[2] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
  epsPn1[2] = 0.0;

  alphan  = 0.0;
  alphan1 = 0.0;

  dg_n1 = 0.0;

  if (SHVs != 0)
    SHVs->Zero();

  return 0;
}

NDMaterial*
J2BeamThread3d::getCopy()
{
  return new J2BeamThread3d(this->getTag(), E, nu, sigmaY, Hiso, Hkin);
}

NDMaterial*
J2BeamThread3d::getCopy(const char* type)
{
  if (strcmp(type, this->getType()) == 0)
    return this->getCopy();

  return nullptr;
}

const char*
J2BeamThread3d::getType() const
{
  return "BeamFiber";
}

int
J2BeamThread3d::getOrder() const
{
  return 3;
}

const Vector&
J2BeamThread3d::getStressSensitivity(int gradIndex, bool conditional)
{
  static Vector sigma(3);

  sigma(0) = 0.0;
  sigma(1) = 0.0;
  sigma(2) = 0.0;

  double dEdh      = 0.0;
  double dsigmaYdh = 0.0;
  double dHkindh   = 0.0;
  double dHisodh   = 0.0;
  double dGdh      = 0.0;

  if (parameterID == 1) { // E
    dEdh = 1.0;
    dGdh = 0.5 / (1.0 + nu);
  }
  if (parameterID == 2) { // nu
    dGdh = -0.5 * E / (1.0 + 2.0 * nu + nu * nu);
  }
  if (parameterID == 5) {
    dsigmaYdh = 1.0;
  }
  if (parameterID == 6) {
    dHkindh = 1.0;
  }
  if (parameterID == 7) {
    dHisodh = 1.0;
  }

  double G = 0.5 * E / (1.0 + nu);

  double depsPdh[3];
  depsPdh[0]      = 0.0;
  depsPdh[1]      = 0.0;
  depsPdh[2]      = 0.0;
  double dalphadh = 0.0;
  if (SHVs != 0) {
    depsPdh[0] = (*SHVs)(0, gradIndex);
    depsPdh[1] = (*SHVs)(1, gradIndex);
    depsPdh[2] = (*SHVs)(2, gradIndex);
    dalphadh   = (*SHVs)(3, gradIndex);
  }

  static const double one3   = 1.0/3.0;
  static const double two3   = 2.0 * one3;
  static const double root23 = std::sqrt(two3);

  double xsi[3];
  xsi[0] = E * (Tepsilon(0) - epsPn1[0]) - Hkin * epsPn1[0];
  xsi[1] = G * (Tepsilon(1) - epsPn1[1]) - one3 * Hkin * epsPn1[1];
  xsi[2] = G * (Tepsilon(2) - epsPn1[2]) - one3 * Hkin * epsPn1[2];

  double q = std::sqrt(two3 * xsi[0] * xsi[0] + 2.0 * xsi[1] * xsi[1] + 2.0 * xsi[2] * xsi[2]);
  double F = q - root23 * (sigmaY + Hiso * alphan1);

  if (F <= -100 * DBL_EPSILON) {
    sigma(0) = dEdh * (Tepsilon(0) - epsPn1[0]) - E * depsPdh[0];
    sigma(1) = dGdh * (Tepsilon(1) - epsPn1[1]) - G * depsPdh[1];
    sigma(2) = dGdh * (Tepsilon(2) - epsPn1[2]) - G * depsPdh[2];
  } else {
    static Matrix J(4, 4);
    static Vector b(4);
    static Vector dx(4);

    double dg = dg_n1;

    J(0, 0) = 1.0 + dg * two3 * (E + Hkin);
    J(0, 1) = 0.0;
    J(0, 2) = 0.0;
    J(1, 0) = 0.0;
    J(1, 1) = 1.0 + dg * (2.0 * G + two3 * Hkin);
    J(1, 2) = 0.0;
    J(2, 0) = 0.0;
    J(2, 1) = 0.0;
    J(2, 2) = 1.0 + dg * (2.0 * G + two3 * Hkin);

    J(0, 3) = two3 * (E + Hkin) * xsi[0];
    J(1, 3) = (2.0 * G + two3 * Hkin) * xsi[1];
    J(2, 3) = (2.0 * G + two3 * Hkin) * xsi[2];

    //J(2,0) = xsi[0]*two3/q; J(2,1) = xsi[1]*2.0/q;
    J(3, 0) = (1.0 - two3 * Hiso * dg) * xsi[0] * two3 / q;
    J(3, 1) = (1.0 - two3 * Hiso * dg) * xsi[1] * 2.0 / q;
    J(3, 2) = (1.0 - two3 * Hiso * dg) * xsi[2] * 2.0 / q;

    //J(2,2) = -root23*Hiso;
    J(3, 3) = -two3 * Hiso * q;

    b(0) = dEdh * Tepsilon(0) - (E + Hkin) * depsPdh[0] - (dEdh + dHkindh) * epsPn1[0];
    b(1) =
        dGdh * Tepsilon(1) - (G + one3 * Hkin) * depsPdh[1] - (dGdh + one3 * dHkindh) * epsPn1[1];
    b(2) =
        dGdh * Tepsilon(2) - (G + one3 * Hkin) * depsPdh[2] - (dGdh + one3 * dHkindh) * epsPn1[2];
    b(3) = root23 * (dsigmaYdh + dHisodh * alphan1 + Hiso * dalphadh);

    J.Solve(b, dx);

    depsPdh[0] += dx(3) * two3 * xsi[0] + dg * two3 * dx(0);
    depsPdh[1] += dx(3) * 2.0 * xsi[1] + dg * 2.0 * dx(1);
    depsPdh[2] += dx(3) * 2.0 * xsi[2] + dg * 2.0 * dx(2);

    sigma(0) = dx(0) + Hkin * depsPdh[0] + dHkindh * epsPn1[0];
    sigma(1) = dx(1) + one3 * Hkin * depsPdh[1] + one3 * dHkindh * epsPn1[1];
    sigma(2) = dx(2) + one3 * Hkin * depsPdh[2] + one3 * dHkindh * epsPn1[2];
  }

  return sigma;
}

int
J2BeamThread3d::commitSensitivity(const Vector& depsdh, int gradIndex, int numGrads)
{
  if (SHVs == 0) {
    SHVs = new Matrix(4, numGrads);
  }

  if (gradIndex >= SHVs->noCols()) {
    return 0;
  }

  //return 0;

  double dEdh      = 0.0;
  double dsigmaYdh = 0.0;
  double dHkindh   = 0.0;
  double dHisodh   = 0.0;
  double dGdh      = 0.0;

  if (parameterID == 1) { // E
    dEdh = 1.0;
    dGdh = 0.5 / (1.0 + nu);
  }
  if (parameterID == 2) { // nu
    dGdh = -0.5 * E / (1.0 + 2.0 * nu + nu * nu);
  }
  if (parameterID == 5) {
    dsigmaYdh = 1.0;
  }
  if (parameterID == 6) {
    dHkindh = 1.0;
  }
  if (parameterID == 7) {
    dHisodh = 1.0;
  }

  double G = 0.5 * E / (1.0 + nu);

  double depsPdh[3];
  depsPdh[0]      = 0.0;
  depsPdh[1]      = 0.0;
  depsPdh[2]      = 0.0;
  double dalphadh = 0.0;
  if (SHVs != 0) {
    depsPdh[0] = (*SHVs)(0, gradIndex);
    depsPdh[1] = (*SHVs)(1, gradIndex);
    depsPdh[2] = (*SHVs)(2, gradIndex);
    dalphadh   = (*SHVs)(3, gradIndex);
  }

  static constexpr double one3   = 1.0 / 3.0;
  static constexpr double two3   = 2.0 * one3;
  static const double root23 = std::sqrt(two3);

  double xsi[3];
  xsi[0] = E * (Tepsilon(0) - epsPn1[0]) - Hkin * epsPn1[0];
  xsi[1] = G * (Tepsilon(1) - epsPn1[1]) - one3 * Hkin * epsPn1[1];
  xsi[2] = G * (Tepsilon(2) - epsPn1[2]) - one3 * Hkin * epsPn1[2];

  double q = std::sqrt(two3 * xsi[0] * xsi[0] + 2.0 * xsi[1] * xsi[1] + 2.0 * xsi[2] * xsi[2]);
  double F = q - root23 * (sigmaY + Hiso * alphan1);

  if (F <= -100 * DBL_EPSILON) {
    // Do nothing
  } else {
    static Matrix J(4, 4);
    static Vector b(4);
    static Vector dx(4);

    double dg = dg_n1;

    J(0, 0) = 1.0 + dg * two3 * (E + Hkin);
    J(0, 1) = 0.0;
    J(0, 2) = 0.0;
    J(1, 0) = 0.0;
    J(1, 1) = 1.0 + dg * (2.0 * G + two3 * Hkin);
    J(1, 2) = 0.0;
    J(2, 0) = 0.0;
    J(2, 1) = 0.0;
    J(2, 2) = 1.0 + dg * (2.0 * G + two3 * Hkin);

    J(0, 3) = two3 * (E + Hkin) * xsi[0];
    J(1, 3) = (2.0 * G + two3 * Hkin) * xsi[1];
    J(2, 3) = (2.0 * G + two3 * Hkin) * xsi[2];

    //J(2,0) = xsi[0]*two3/q; J(2,1) = xsi[1]*2.0/q;
    J(3, 0) = (1.0 - two3 * Hiso * dg) * xsi[0] * two3 / q;
    J(3, 1) = (1.0 - two3 * Hiso * dg) * xsi[1] * 2.0 / q;
    J(3, 2) = (1.0 - two3 * Hiso * dg) * xsi[2] * 2.0 / q;

    //J(2,2) = -root23*Hiso;
    J(3, 3) = -two3 * Hiso * q;

    b(0) =
        E * depsdh(0) + dEdh * Tepsilon(0) - (E + Hkin) * depsPdh[0] - (dEdh + dHkindh) * epsPn1[0];
    b(1) = G * depsdh(1) + dGdh * Tepsilon(1) - (G + one3 * Hkin) * depsPdh[1] -
           (dGdh + one3 * dHkindh) * epsPn1[1];
    b(2) = G * depsdh(2) + dGdh * Tepsilon(2) - (G + one3 * Hkin) * depsPdh[2] -
           (dGdh + one3 * dHkindh) * epsPn1[2];
    b(3) = root23 * (dsigmaYdh + dHisodh * alphan1 + Hiso * dalphadh);

    J.Solve(b, dx);

    dalphadh +=
        dx(3) * root23 * q +
        dg * root23 * (xsi[0] * two3 * dx(0) + xsi[1] * 2.0 * dx(1) + xsi[2]*2.0*dx(2)) / q;
    depsPdh[0] += dx(3) * two3 * xsi[0] + dg * two3 * dx(0);
    depsPdh[1] += dx(3) * 2.0 * xsi[1] + dg * 2.0 * dx(1);
    depsPdh[2] += dx(3) * 2.0 * xsi[2] + dg * 2.0 * dx(2);

    (*SHVs)(0, gradIndex) = depsPdh[0];
    (*SHVs)(1, gradIndex) = depsPdh[1];
    (*SHVs)(2, gradIndex) = depsPdh[2];
    (*SHVs)(3, gradIndex) = dalphadh;
  }

  return 0;
}


int
J2BeamThread3d::sendSelf(int commitTag, Channel& theChannel)
{
  int res = 0;

  static Vector data(6);

  data(0) = this->getTag();
  data(1) = E;
  data(2) = nu;
  data(3) = sigmaY;
  data(4) = Hiso;
  data(5) = Hkin;

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "J2BeamThread3d::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

int
J2BeamThread3d::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  int res = 0;

  static Vector data(6);

  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "J2BeamThread3d::recvSelf -- could not recv Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  E      = data(1);
  nu     = data(2);
  sigmaY = data(3);
  Hiso   = data(4);
  Hkin   = data(5);

  return res;
}

void
J2BeamThread3d::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "J2 Beam Fiber Material Model" << "\n";
    s << "\tE:  " << E << "\n";
    s << "\tnu:  " << nu << "\n";
    s << "\tsigmaY:  " << sigmaY << "\n";
    s << "\tHiso:  " << Hiso << "\n";
    s << "\tHkin:  " << Hkin << "\n";
  }
  else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"J2BeamFiber\", ";
    s << "\"E\": " << E << ", ";
    s << "\"nu\": " << nu << ", ";
    s << "\"sigmaY\": " << sigmaY << ", ";
    s << "\"Hiso\": " << Hiso << ", ";
    s << "\"Hkin\": " << Hkin;
    s << "}";
  }
  return;
}

int
J2BeamThread3d::setParameter(const char** argv, int argc, Parameter& param)
{
  if (strcmp(argv[0], "E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  } else if (strcmp(argv[0], "nu") == 0) {
    param.setValue(nu);
    return param.addObject(2, this);
  } else if (strcmp(argv[0], "sigmaY") == 0 || strcmp(argv[0], "fy") == 0 ||
             strcmp(argv[0], "Fy") == 0) {
    param.setValue(sigmaY);
    return param.addObject(5, this);
  } else if (strcmp(argv[0], "Hkin") == 0) {
    param.setValue(Hkin);
    return param.addObject(6, this);
  } else if (strcmp(argv[0], "Hiso") == 0) {
    param.setValue(Hiso);
    return param.addObject(7, this);
  }

  return -1;
}

int
J2BeamThread3d::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
  case 1:  E = info.theDouble; return 0;
  case 2:  nu = info.theDouble; return 0;
  case 5:  sigmaY = info.theDouble; return 0;
  case 6:  Hkin = info.theDouble; return 0;
  case 7:  Hiso = info.theDouble; return 0;
  default: return -1;
  }
}

int
J2BeamThread3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
// Please cite the following resource in any derivative works:
//
// Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//

//
// Written: Claudio M. Perez, 
//          Filip C. Filippou
//          University of California, Berkeley
//
// Developed with FEDEASLab [2].
//
// References:
// 
// [2] Filippou, F.C. (1998)
//     "FEDEASLab: Finite Elements for Design Evaluation and Analysis of Structures"
//
// [3] Battini, J.-M. and Pacoste, C. (2002) 
//     "Co-rotational beam elements with warping effects in instability problems", 
//     Computer Methods in Applied Mechanics and Engineering, 191(17–18), 
//     pp. 1755–1789. Available at: https://doi.org/10.1016/S0045-7825(01)00352-8.
//
#pragma once
#include <array>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <MatrixND.h>
#include <GroupSO3.h>
#include "EuclidIsometry.h"

class Node;

namespace OpenSees {

template <int nn>
class BattiniIsometry : public AlignedIsometry<nn>
{
public:
  BattiniIsometry(const Vector3D& vecxz)
  : AlignedIsometry<nn>(vecxz), 
    n(0)
  {
  }

  using AlignedIsometry<nn>::init;
  using AlignedIsometry<nn>::pres;

  Matrix3D
  update_basis(const Matrix3D& RI, const Matrix3D& RJ, const Vector3D& dx) 
  noexcept final
  {
    Matrix3D R;
    {
      Vector3D e1 = dx;
      e1 /= e1.norm();

      constexpr static Vector3D D2 {0,1,0};
      const Vector3D E2 = this->AlignedIsometry<nn>::R[init]*D2;
      q = RI*E2; //*R[init];
      q.addVector(0.5, RJ*E2, 0.5);


      Vector3D e3 = e1.cross(q);
      e3 /= e3.norm();

      const Vector3D e2 = e3.cross(e1);

      for (int i = 0; i < 3; i++) {
        R(i,0) = e1[i];
        R(i,1) = e2[i];
        R(i,2) = e3[i];
      }
    
      const Vector3D Q = R^q;
      const double Q2 = Q[1];

      n = Q[0]/Q2;
      // n = Q[0]/Q[1];

      Vector3D QI = R^(RI*E2);
      Vector3D QJ = R^(RJ*E2);
      // n11 = QI[0]/Q[1];
      // n12 = QI[1]/Q[1];
      // n21 = QJ[0]/Q[1];
      // n22 = QJ[1]/Q[1];
      n11 = QI[0]/Q2;
      n12 = QI[1]/Q2;
      n13 = QI[2]/Q2;

      n21 = QJ[0]/Q2;
      n22 = QJ[1]/Q2;
      n23 = QJ[2]/Q2;
    }
    return R;
  }


  MatrixND<3,6> 
  getRotationGradient(int node) final {
    MatrixND<3,6> Gb{};

    constexpr Vector3D axis{1, 0, 0};
    constexpr Matrix3D ix = Hat(axis);

    double Ln = this->getLength();

    if (node == 0) {
      Gb.template insert<0,0>( ix, -1.0/Ln);
      Gb(0,2) =    n/Ln;
      Gb(0,3) =  n12/2.0; // - n;
      Gb(0,4) = -n11/2.0;
    }
    else if (node == nn-1) {
      Gb.template insert<0,0>( ix,   1.0/Ln);
      Gb(0,2) = -  n/Ln;
      Gb(0,3) =  n22/2.0;
#if 0 // TODO
// original impl
      Gb(0,4) = -n12/2.0;
#else
      Gb(0,4) = -n21/2.0;
#endif
    }
  
    return Gb;
  }


  MatrixND<6*nn,6*nn>
  getHessian(const VectorND<6>& pw) final
  {
    MatrixND<6*nn,6*nn> H{};

    static_assert(nn >= 2, "Requires at least two nodes.");

    constexpr int I = 0;
    constexpr int J = nn - 1;

    const double L = this->getLength();

    const Vector3D m{pw[3], pw[4], pw[5]};
    const double m1 = m[0];

    constexpr static Vector3D 
                   e1{1.0, 0.0, 0.0},
                   e2{0.0, 1.0, 0.0},
                   e3{0.0, 0.0, 1.0};

    using Row6 = std::array<double,6>;

    auto sigma = [](int node) -> double {
      if (node == I)
        return -1.0;
      if (node == J)
        return  1.0;
      return 0.0;
    };

    auto addGRow = [](Row6& row,
                      const MatrixND<3,6>& G,
                      int r,
                      double scale) {
      for (int c = 0; c < 6; ++c)
        row[c] += scale * G(r,c);
    };

    auto addScaledRow = [](Row6& row,
                          const Row6& src,
                          double scale) {
      for (int c = 0; c < 6; ++c)
        row[c] += scale * src[c];
    };

    auto splitX = [](const Row6& row) -> Vector3D {
      return Vector3D{row[0], row[1], row[2]};
    };

    auto splitT = [](const Row6& row) -> Vector3D {
      return Vector3D{row[3], row[4], row[5]};
    };

    auto lambdaRow = [&](int b, const MatrixND<3,6>& Gw) -> Row6 {
      Row6 lam{};

      /*
      * lambda_b = dQ2/Q2
      *
      * =
      * 1/2 delta_bI (nu_I1 theta_3 - rho_I theta_1)
      * +
      * 1/2 delta_bJ (nu_J1 theta_3 - rho_J theta_1)
      * -
      * n domega_3.
      */
      addGRow(lam, Gw, 2, -n);

      if (b == I) {
        lam[3] += -0.5 * n13;
        lam[5] +=  0.5 * n11;
      }
      else if (b == J) {
        lam[3] += -0.5 * n23;
        lam[5] +=  0.5 * n21;
      }

      return lam;
    };

    auto directorRows = [&](int K,
                            int b,
                            const MatrixND<3,6>& Gw,
                            const Row6& lam,
                            Row6& dnu1,
                            Row6& dnu2) {
      dnu1 = {};
      dnu2 = {};

      double nu1 = 0.0;
      double nu2 = 0.0;
      double rho = 0.0;

      if (K == I) {
        nu1 = n11;
        nu2 = n12;
        rho = n13;
      }
      else {
        nu1 = n21;
        nu2 = n22;
        rho = n23;
      }

      /*
      * d nu_1 =
      * delta_bK(theta_2 rho - theta_3 nu_2)
      * - (omega_2 rho - omega_3 nu_2)
      * - nu_1 lambda.
      */
      if (b == K) {
        dnu1[4] +=  rho;
        dnu1[5] += -nu2;
      }

      addGRow(dnu1, Gw, 1, -rho);
      addGRow(dnu1, Gw, 2,  nu2);
      addScaledRow(dnu1, lam, -nu1);

      /*
      * d nu_2 =
      * delta_bK(theta_3 nu_1 - theta_1 rho)
      * - (omega_3 nu_1 - omega_1 rho)
      * - nu_2 lambda.
      */
      if (b == K) {
        dnu2[5] +=  nu1;
        dnu2[3] += -rho;
      }

      addGRow(dnu2, Gw, 2, -nu1);
      addGRow(dnu2, Gw, 0,  rho);
      addScaledRow(dnu2, lam, -nu2);
    };

    /*
    * B = ([e1]_x - n e1 e3^T) / L.
    *
    * B^T m = [ 0,
    *           m3/L,
    *          -(m2 + n m1)/L ].
    */
    const Vector3D Btm{
      0.0,
      m[2] / L,
      -(m[1] + n*m[0]) / L
    };

    const int end_nodes[2] = {I, J};

    for (int ib = 0; ib < 2; ++ib) {
      const int b = end_nodes[ib];
      const double sb = sigma(b);

      const MatrixND<3,6> Gw = this->getRotationGradient(b);
      const Row6 lam = lambdaRow(b, Gw);

      Row6 dI1{}, dI2{};
      Row6 dJ1{}, dJ2{};

      directorRows(I, b, Gw, lam, dI1, dI2);
      directorRows(J, b, Gw, lam, dJ1, dJ2);

      //
      // dn = 1/2 (d nu_I1 + d nu_J1).
      //
      Row6 dn{};
      for (int c = 0; c < 6; ++c)
        dn[c] = 0.5 * (dI1[c] + dJ1[c]);

      const Vector3D dnx = splitX(dn);
      const Vector3D dnt = splitT(dn);

      for (int ia = 0; ia < 2; ++ia) {
        const int a = end_nodes[ia];
        const double sa = sigma(a);

        Matrix3D Hxx{};
        Matrix3D Hxt{};
        Matrix3D Htx{};
        Matrix3D Htt{};

        /*
        * Common aligned-family translational-row blocks:
        *
        * H_ab^{xx}
        * =
        * -sigma_a [
        *   sigma_b/L (B^T m) e1^T
        *   + m1/L e3 dn_x^T
        * ].
        *
        * H_ab^{x theta} = -sigma_a m1/L e3 dn_theta^T.
        */
        Hxx.addTensorProduct(Btm, e1, -sa * sb / L);
        Hxx.addTensorProduct(e3, dnx, -sa * m1 / L);

        Hxt.addTensorProduct(e3, dnt, -sa * m1 / L);

        //
        // Battini rotational-row blocks:
        //
        // H_Kb^{theta x} = 1/2 m1 (e1 dnu_K2_x^T - e2 dnu_K1_x^T)
        //
        // H_Kb^{theta theta}  = 1/2 m1 (e1 dnu_K2_theta^T - e2 dnu_K1_theta^T).
        //
        const Row6& dK1 = (a == I) ? dI1 : dJ1;
        const Row6& dK2 = (a == I) ? dI2 : dJ2;


        Htx.addTensorProduct(e1, splitX(dK2),  0.5 * m1);
        Htx.addTensorProduct(e2, splitX(dK1), -0.5 * m1);

        Htt.addTensorProduct(e1, splitT(dK2),  0.5 * m1);
        Htt.addTensorProduct(e2, splitT(dK1), -0.5 * m1);

        H.assemble(Hxx, 6*a + 0, 6*b + 0, 1.0);
        H.assemble(Hxt, 6*a + 0, 6*b + 3, 1.0);
        H.assemble(Htx, 6*a + 3, 6*b + 0, 1.0);
        H.assemble(Htt, 6*a + 3, 6*b + 3, 1.0);
      }
    }

    return H;
  }

  MatrixND<nn*6,nn*6>
  getRotationJacobian(const VectorND<nn*6>&pwx) final 
  {
    if constexpr (nn != 2) {
      return MatrixND<nn*6,nn*6> {};
    }
    else {
      // TODO: Copied from RankinIsometry.h, need to adapt to Battini
      MatrixND<3,12> NWL{};
      const double Ln = this->getLength();

      constexpr static Matrix3D ex = Hat(Vector3D {1,0,0});

      for (int i=0; i<nn; i++)
        NWL.assemble(Hat(&pwx[i*6]), 0, i*6,  -1.0);


      MatrixND<12,3> Gamma{};
      Gamma.template insert<0,0>(ex,  1.0);
      Gamma(3,0) = -1.0;
      Gamma.template insert<6,0>(ex, -1.0);

      MatrixND<12,3> Psi{};
      Psi.template insert<3,0>(Eye3, 1.0);
      Psi.template insert<6,0>(ex,   -Ln);
      Psi.template insert<9,0>(Eye3, 1.0);

      const Matrix3D B = Gamma^Psi;
      Matrix3D A;
      B.invert(A);
      return Gamma*A.transpose()*NWL;
    }
  }

private:
  Vector3D q {0,1,0};
  // double n   = 0,
  //        n11 = 0,
  //        n12 = 1,
  //        n21 = 0,
  //        n22 = 1;

  double n   = 0,
         n11 = 0,
         n12 = 1,
         n13 = 0,
         n21 = 0,
         n22 = 1,
         n23 = 0;
};
} // namespace OpenSees
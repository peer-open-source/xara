//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//----------------------------------------------------------------------------//
//
// Please cite the following resource in any derivative works:
//
// [1] Perez, C.M., and Filippou F.C. (2024)
//     "On Nonlinear Geometric Transformations of Finite Elements" 
//     Int. J. Numer. Meth. Engrg.; https://doi.org/10.1002/nme.7506
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
// [1] Perez, C.M., and Filippou F.C. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024;
//     https://doi.org/10.1002/nme.7506
//
// [2] Filippou, F.C. (1998)
//     "FEDEASLab: Finite Elements for Design Evaluation and Analysis of Structures"
//
// [3] Nour-Omid, B. and Rankin, C.C. (1991) "Finite rotation analysis and 
//     consistent linearization using projectors", 
//     Computer Methods in Applied Mechanics and Engineering, 93(3), pp. 353–384. 
//     Available at: https://doi.org/10.1016/0045-7825(91)90248-5.
//
#pragma once
#include <Vector3D.h>
#include <Matrix3D.h>
#include <MatrixND.h>
#include "EuclidIsometry.h"

class Node;

namespace OpenSees {

template <int nn>
class RankinIsometry : public AlignedIsometry<nn>
{
public:
  RankinIsometry(const Vector3D& vecxz)
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
      const Vector3D e1 = dx/dx.norm();

      constexpr static Vector3D D2 {0,1,0};
      const Vector3D E2 = this->AlignedIsometry<nn>::R[init]*D2;
      q = RI*E2;
      Vector3D e3 = e1.cross(q);
      e3 /= e3.norm();

      const Vector3D e2 = e3.cross(e1);

      for (int i = 0; i < 3; i++) {
        R(i,0) = e1[i];
        R(i,1) = e2[i];
        R(i,2) = e3[i];
      }

      const Vector3D Q = R^q;
      n = Q[0]/Q[1];
    }
    return R;
  }

  MatrixND<6*nn,6*nn>
  getHessian(const VectorND<6>& pw) final
  {
    MatrixND<6*nn,6*nn> H{};

    static_assert(nn >= 2, "Requires at least two nodes.");

    constexpr int I = 0;
    constexpr int J = nn - 1;

    const double L = this->getLength();
    if (L == 0.0)
      return H;

    const Vector3D m{pw[3], pw[4], pw[5]};
    const double m1 = m[0];

    const Vector3D e1{1.0, 0.0, 0.0};
    const Vector3D e2{0.0, 1.0, 0.0};
    const Vector3D e3{0.0, 0.0, 1.0};

    auto sigma = [](int node) -> double {
      if (node == I)
        return -1.0;
      if (node == J)
        return  1.0;
      return 0.0;
    };

    //
    // B = ([e1]_x - n e1 e3^T) / L.
    //
    // Therefore:
    //
    //   B^T m = [ 0,
    //             m3/L,
    //            -(m2 + n m1)/L ].
    //
    const Vector3D Btm{
      0.0,
      m[2] / L,
      -(m[1] + n*m[0]) / L
    };

    const double nfac = 1.0 + n*n;

    const int end_nodes[2] = {I, J};

    for (int ia = 0; ia < 2; ++ia) {
      const int a = end_nodes[ia];
      const double sa = sigma(a);

      for (int ib = 0; ib < 2; ++ib) {
        const int b = end_nodes[ib];
        const double sb = sigma(b);

        /*
        * dn_b = dn_x^T v_b + dn_theta^T theta_b
        *
        * dn_x     = (1+n^2) sigma_b/L e2
        * dn_theta = -(1+n^2) delta_{bI} e3
        */
        Vector3D dnx{};
        Vector3D dnt{};

        dnx[1] = nfac * sb / L;

        if (b == I)
          dnt[2] = -nfac;

        Matrix3D Hxx{},Hxt{},Htx{},Htt{};

        /*
        * H_ab^{xx}
        * =
        * -sigma_a [
        *   sigma_b/L (B^T m) e1^T
        *   + m1/L e3 dn_x^T
        * ].
        */
        Hxx.addTensorProduct(Btm, e1, -sa * sb / L);
        Hxx.addTensorProduct(e3, dnx, -sa * m1 / L);

        /*
        * H_ab^{x theta}
        * =
        * -sigma_a m1/L e3 dn_theta^T.
        */
        Hxt.addTensorProduct(e3, dnt, -sa * m1 / L);

        /*
        * Rankin C_a is nonzero only for a = I:
        *
        * H_ab^{theta x} = -delta_{aI} m1 e2 dn_x^T
        *
        * H_ab^{theta theta}
        * =
        * -delta_{aI} m1 e2 dn_theta^T.
        */
        if (a == I) {
          Htx.addTensorProduct(e2, dnx, -m1);
          Htt.addTensorProduct(e2, dnt, -m1);
        }

        H.assemble(Hxx, 6*a + 0, 6*b + 0, 1.0);
        H.assemble(Hxt, 6*a + 0, 6*b + 3, 1.0);
        H.assemble(Htx, 6*a + 3, 6*b + 0, 1.0);
        H.assemble(Htt, 6*a + 3, 6*b + 3, 1.0);
      }
    }

    return H;
  }

  Matrix3D 
  getRotationSensitivity(
    std::array<Node*,nn> nodes
  ) final {
    Matrix3D dR{};

    std::array<int,nn> ix{};
    for (int i=0; i<nn; i++)
      ix[i] = nodes[i]->getCrdsSensitivity();
    const Matrix3D RI = MatrixFromVersor(nodes[0]->getTrialRotation());
    //
    // Coordinate sensitivity of the end-to-end reference chord.
    // The current chord has the same direct coordinate sensitivity,
    // with displacements, rotations, and offsets held fixed.
    //
    Vector3D dXdh{};

    auto addCoordinateSensitivity = [](Vector3D& v, int dir, double scale) {
      if (dir >= 1 && dir <= 3)
        v[dir - 1] += scale;
    };

    addCoordinateSensitivity(dXdh, ix[nn-1],  1.0);
    addCoordinateSensitivity(dXdh, ix[0],    -1.0);

    auto dot = [](const Vector3D& a, const Vector3D& b) {
      double s = 0.0;
      for (int i = 0; i < 3; ++i)
        s += a[i]*b[i];
      return s;
    };

    auto projectNormal = [&](const Vector3D& v, const Vector3D& a) -> Vector3D {
      Vector3D r = v;
      const double av = dot(a, v);
      for (int i = 0; i < 3; ++i)
        r[i] -= a[i]*av;
      return r;
    };

    if (dot(dXdh, dXdh) == 0.0)
      return dR;

    constexpr static Vector3D D1 {1,0,0};
    constexpr static Vector3D D2 {0,1,0};
    constexpr static Vector3D D3 {0,0,1};

    const Matrix3D& R0 = this->R[init];
    const Matrix3D& Rt = this->R[pres];

    const Vector3D E1 = R0*D1;
    const Vector3D E2 = R0*D2;

    const Vector3D e1 = Rt*D1;
    const Vector3D e3 = Rt*D3;

    const double L0 = this->L;
    const double Lt = this->getLength();

    const Vector3D y0 = this->vz.cross(E1);
    const double y0n = y0.norm();

    const Vector3D s = e1.cross(q);
    const double sn = s.norm();

    if (L0 == 0.0 || Lt == 0.0 || y0n == 0.0 || sn == 0.0)
      return dR;

    //
    // Initial basis sensitivity.
    //
    Vector3D dE1 = projectNormal(dXdh, E1);
    for (int i = 0; i < 3; ++i)
      dE1[i] /= L0;

    Vector3D dE2 = projectNormal(this->vz.cross(dE1), E2);
    for (int i = 0; i < 3; ++i)
      dE2[i] /= y0n;

    //
    // Current basis sensitivity.
    //
    Vector3D de1 = projectNormal(dXdh, e1);
    for (int i = 0; i < 3; ++i)
      de1[i] /= Lt;

    const Vector3D dq = RI*dE2;

    const Vector3D ds_a = de1.cross(q);
    const Vector3D ds_b = e1.cross(dq);

    Vector3D ds{};
    for (int i = 0; i < 3; ++i)
      ds[i] = ds_a[i] + ds_b[i];

    Vector3D de3 = projectNormal(ds, e3);
    for (int i = 0; i < 3; ++i)
      de3[i] /= sn;

    const Vector3D de2_a = de3.cross(e1);
    const Vector3D de2_b = e3.cross(de1);

    Vector3D de2{};
    for (int i = 0; i < 3; ++i)
      de2[i] = de2_a[i] + de2_b[i];

    for (int i = 0; i < 3; ++i) {
      dR(i,0) = de1[i];
      dR(i,1) = de2[i];
      dR(i,2) = de3[i];
    }

    return dR;
  }

  MatrixND<6*nn,6*nn>
  getRotationJacobian(const VectorND<6*nn>&pwx) final 
  {
    if constexpr (nn != 2) {
      // TODO: Implement for nn != 2
      return MatrixND<6*nn,6*nn> {};
    }
    else {
      
      MatrixND<3,12> NWL{};
      const double Ln = this->getLength();

      constexpr static Matrix3D ex = Hat(Vector3D {1,0,0});

      for (int i=0; i<nn; i++)
        NWL.assemble(Hat(&pwx[i*6]), 0, i*6,  -1.0);

  #if __cplusplus >= 202000L
      static constinit MatrixND<12,3> Gamma = MakeGamma();
      static constinit MatrixND<12,3> Psi0  = MakePsi();
      MatrixND<12,3> Psi = Psi0;
      Psi.template insert<6,0>(ex,  -Ln);
  #else
      MatrixND<12,3> Gamma{};
      Gamma.template insert<0,0>(ex,  1.0);
      Gamma(3,0) = -1.0;
      Gamma.template insert<6,0>(ex, -1.0);

      MatrixND<12,3> Psi{};
      Psi.template insert<3,0>(Eye3, 1.0);
      Psi.template insert<6,0>(ex,   -Ln);
      Psi.template insert<9,0>(Eye3, 1.0);
  #endif
      const Matrix3D B = Gamma^Psi;
      Matrix3D A;
      B.invert(A);
      return Gamma*A.transpose()*NWL;
    }
  }

  MatrixND<3,6> 
  getRotationGradient(int node) final {
    MatrixND<3,6> Gb{};

    constexpr Vector3D axis{1, 0, 0};
    constexpr Matrix3D ix = Hat(axis);

    const double Ln = this->getLength();

    if (node == 0) {
      Gb.template insert<0,0>( ix, -1.0/Ln);
      Gb(0,2) =  n/Ln;
      Gb(0,3) =   1.0;
      Gb(0,4) =    -n;
    }
    else if (node == nn-1) {
      Gb.template insert<0,0>( ix,  1.0/Ln);
      Gb(0,2) = -n/Ln;
      Gb(0,3) =  0.0;
    }
    return Gb;
  }

private:
  MatrixND<3,6>
  getBasisVariation(int ie, int node)
  {
    MatrixND<3,6> dei{};
    if (ie == 1) {
      Matrix3D A{};
      A(1,1) = A(2,2) = 1.0/this->getLength();
      if (node == 0) {
        dei.template insert<0,0>(A, -1.0);
      }
      else if (node == nn-1) {
        dei.template insert<0,0>(A,  1.0);
      }
    }

    else if (ie == 3) {
      Matrix3D A{};
      Vector3D q = this->getRotation()^this->q;
      Vector3D v = q.cross(Vector3D{1,0,0});

      A(0,0) = A(1,1) = 1.0/v.norm();
      Matrix3D Q = Hat(q);
      dei = A*(
        Hat(Vector3D{1,0,0})*Q
       +Q*this->getBasisVariation(1, node)
      );
    }
    return dei;
  }

private:
  Vector3D q;
  double n   = 0;


#if __cplusplus >= 202000L
  static inline consteval MatrixND<12,3> 
  MakePsi()
  {
    MatrixND<12,3> Psi{};
    Psi.template insert<3,0>(Eye3, 1.0);
    Psi.template insert<9,0>(Eye3, 1.0);
    return Psi;
  }

  static inline consteval MatrixND<12,3> 
  MakeGamma()
  {
    MatrixND<12,3> Gamma{};
    constexpr Matrix3D ex = Hat(Vector3D {1,0,0});
    Gamma.template insert<0,0>(ex,  1.0);
    Gamma(3,0) = -1.0;
    Gamma.template insert<6,0>(ex, -1.0);
    return Gamma;
  }
#endif

};
} // namespace OpenSees

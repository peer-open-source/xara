//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//----------------------------------------------------------------------------//
//
//                                 FEDEASLab
//       Finite Elements for Design Evaluation and Analysis of Structures
//
//----------------------------------------------------------------------------//
//
// Please cite the following resource in any derivative works:
//
// [1] Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//

#pragma once
#include <MatrixND.h>
#include <Matrix3D.h>
#include <Triad.h>
#include <Vector3D.h>
#include <Rotations.hpp>
#include "EuclidIsometry.h"

namespace OpenSees {

template <int nn>
class CrisfieldIsometry : public AlignedIsometry<nn> {
public:
  CrisfieldIsometry(const Vector3D& vecxz) : AlignedIsometry<nn>{vecxz} {}


  Matrix3D
  update_basis(const Matrix3D& RI, const Matrix3D& RJ, const Vector3D& dx) final
  {

    {
      Vector3D gammaw = CayleyFromVersor(
        VersorFromMatrix(RI*RJ.transpose())
        // qI.mult_conj(qJ)
      );

      gammaw *= 0.5;
      Rbar = CaySO3(gammaw)*RI;
  //  Qbar = VersorProduct(VersorFromMatrix(CaySO3(gammaw)), qI);
      Qbar = VersorFromMatrix(Rbar);
      Triad r{Rbar};
      r1 = r[1];
      r2 = r[2];
      r3 = r[3];
    }

    //
    // Compute the base vectors e2, e3
    //
    {
      // 'rotate' the mean rotation matrix Rbar on to e1 to
      // obtain e2 and e3 (using the 'mid-point' procedure)
      //
      // Vector3D e1, e2, e3;
      e[0]  = dx;
      e[0] /= dx.norm();

      // e2 = r2 - (e1 + r1)*((r2^e1)*0.5);
  
      Vector3D tmp;
      tmp  = e[0];
      tmp += r1;//Qbar.rotate(E1);
  
      e[1] = tmp;
      {
        // const Vector3D r2 = Qbar.rotate(E2);
        e[1] *= 0.5*r2.dot(e[0]);
        e[1].addVector(-1.0,  r2, 1.0);
      }
  
      // e3 = r3 - (e1 + r1)*((r3^e1)*0.5);
      e[2] = tmp;
      {
        // const Vector3D r3 = Qbar.rotate(E3);
        e[2] *= r3.dot(e[0])*0.5;
        e[2].addVector(-1.0,  r3, 1.0);
      }
    }

    Matrix3D E;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        E(i,j) = e[j][i];
    return E;
  }


  MatrixND<3,6> 
  getRotationGradient(int node) final {
    MatrixND<3,6> Gb{};

    constexpr Vector3D axis{1, 0, 0};
    constexpr Matrix3D ix = Hat(axis);

    double Ln = this->getLength();
    double n = 0;
    // constexpr Matrix3D ioi = axis.bun(axis);

    if (node == 0) {
      Gb.template insert<0,0>( ix, -1.0/Ln);
      // Gb.template insert<0,3>(ioi,  1.0-n);
      Gb(0,3) =   1.0; // - n;
      Gb(0,4) =    -n;
      Gb(0,2) =  n/Ln;
    }
    else if (node == nn-1) {
      Gb.template insert<0,0>( ix,   1.0/Ln);
      Gb(0,2) = -n/Ln;
    }
    
    return Gb;
  }

  //
  //
  //
  constexpr const Vector3D& 
  getBasisE1() const noexcept
  {
    return e[0];
  }

  constexpr const Vector3D&
  getBasisE2() const noexcept
  {
    return e[1];
  }

  constexpr const Vector3D&
  getBasisE3() const noexcept
  {
    return e[2];
  }

  const Versor&
  getReference()
  {
    return Qbar;
  }

  static inline void
  getLMatrix(const Matrix3D& A,
             const Vector3D& e1, const Vector3D& r1, const Vector3D &ri, MatrixND<12,3>& L)
  {
    Matrix3D L1, L2;
    Matrix3D rie1r1;
    Matrix3D e1e1r1;

    const double rie1 = ri.dot(e1);

    for (int k = 0; k < 3; k++) {
      const double e1r1k = (e1[k] + r1[k]);
      for (int j = 0; j < 3; j++) {
        rie1r1(j,k) = ri[j]*e1r1k;
        e1e1r1(j,k) = e1[j]*e1r1k;
      }
    }

    // L1  = ri'*e1 * A/2 + A*ri*(e1 + r1)'/2;
    L1.zero();
    L1.addMatrix(A, rie1*0.5);
    L1.addMatrixProduct(A, rie1r1, 0.5);

    // L2  = Sri/2 - ri'*e1*S(r1)/4 - Sri*e1*(e1 + r1)'/4;
    L2.zero();
    L2.addSpin(ri, 0.5);
    L2.addSpin(r1, -rie1/4.0);
    L2.addSpinMatrixProduct(ri, e1e1r1, -0.25);

    // L = [L1
    //      L2
    //     -L1
    //      L2];

    L.zero();
    L.assemble(L1, 0, 0,  1.0);
    L.assemble(L2, 3, 0,  1.0);
    L.assemble(L1, 6, 0, -1.0);
    L.assemble(L2, 9, 0,  1.0);

  }

  static inline const MatrixND<12,12> &
  getKs2Matrix(Matrix3D& A, 
               const Vector3D& e1, const Vector3D& r1, const double Ln, const Vector3D &ri, const Vector3D &z)
  {
    static MatrixND<12,12> ks2;

    //  Ksigma2 = [ K11   K12 -K11   K12
    //              K12'  K22 -K12'  K22
    //             -K11  -K12  K11  -K12
    //              K12'  K22 -K12'  K22];

    // U = (-1/2)*A*z*ri'*A + ri'*e1*A*z*e1'/(2*Ln)+...
    //      z'*(e1+r1)*A*ri*e1'/(2*Ln);

    const double zte1  =  z.dot(e1);
    const double ztr1  =  z.dot(r1);

    Matrix3D zrit, ze1t;
    Matrix rizt(3,3), rie1t(3,3);
    Matrix3D e1zt;

    //  const Matrix3D e1zt = e1.bun(z);

    // Chrystal's looping order
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        zrit(i,j)  = z[i]*ri[j];
        rizt(i,j)  = ri[i]*z[j];
        ze1t(i,j)  = z[i]*e1[j];
        e1zt(i,j)  = e1[i]*z[j];
        rie1t(i,j) = ri[i]*e1[j];
      }
    }


    const double C1 = ri.dot(e1)/(2.0*Ln);

    Matrix3D U;
    U.addMatrixTripleProduct(0.0, A, zrit, -0.5);
    U.addMatrixProduct(A, ze1t,   C1);
    U.addMatrixProduct(A, rie1t, (zte1 + ztr1)/(2.0*Ln));

    Matrix3D ks{};
    static Matrix3D m1;

    // K11 = U + U' + ri'*e1*(2*(e1'*z)+z'*r1)*A/(2*Ln);
    ks.addMatrix(U, 1.0);

    // Add matrix U transpose
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        ks(i,j) += U(j,i);

    ks.addMatrix(A, (2.0*zte1 + ztr1)*C1);

    ks2.zero();
    ks2.assemble(ks, 0, 0,  1.0);
    ks2.assemble(ks, 0, 6, -1.0);
    ks2.assemble(ks, 6, 0, -1.0);
    ks2.assemble(ks, 6, 6,  1.0);

    // K12 = (1/4)*(-A*z*e1'*Sri - A*ri*z'*Sr1 - z'*(e1+r1)*A*Sri);
    m1.zero();
    m1.addMatrixProduct(A, ze1t, -1.0);
    ks.zero();
    ks.addMatrixSpinProduct(m1, ri, 0.25);

    m1.zero();
    m1.addMatrixProduct(A, rizt, -1.0);
    ks.addMatrixSpinProduct(m1, r1, 0.25);
    ks.addMatrixSpinProduct(A, ri, -0.25*(zte1+ztr1));

    ks2.assemble(ks, 0, 3,  1.0);
    ks2.assemble(ks, 0, 9,  1.0);
    ks2.assemble(ks, 6, 3, -1.0);
    ks2.assemble(ks, 6, 9, -1.0);

    ks2.assembleTranspose(ks, 3, 0,  1.0);
    ks2.assembleTranspose(ks, 3, 6, -1.0);
    ks2.assembleTranspose(ks, 9, 0,  1.0);
    ks2.assembleTranspose(ks, 9, 6, -1.0);

    // K22 = (1/8)*((-ri'*e1)*Sz*Sr1 + Sr1*z*e1'*Sri + ...
    //       Sri*e1*z'*Sr1 - (e1+r1)'*z*S(e1)*Sri + 2*Sz*Sri);

    ks.zero();
    ks.addSpinProduct(z, r1, -0.125*ri.dot(e1));

    m1.zero();
    m1.addSpinMatrixProduct( r1, ze1t, 1.0);
    ks.addMatrixSpinProduct( m1, ri, 0.125);

    m1.zero();
    m1.addSpinMatrixProduct(ri, e1zt, 1.0);
    ks.addMatrixSpinProduct(m1, r1, 0.125);

    ks.addSpinProduct(e1, ri, -0.125*(zte1 + ztr1));
    ks.addSpinProduct( z, ri, 0.25);

    // Ksigma2 = [ K11   K12 -K11   K12;
    //             K12t  K22 -K12t  K22;
    //            -K11  -K12  K11  -K12;
    //             K12t  K22 -K12t  K22];

    ks2.assemble(ks, 3, 3, 1.0);
    ks2.assemble(ks, 3, 9, 1.0);
    ks2.assemble(ks, 9, 3, 1.0);
    ks2.assemble(ks, 9, 9, 1.0);

    return ks2;
  }

private:
    Versor Qbar;
    Matrix3D Rbar;
    Vector3D r1, r2, r3;
    Vector3D e[3];
};
}
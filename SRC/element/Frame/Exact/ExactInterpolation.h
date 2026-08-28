#pragma once
#include <cstddef>
#include <MatrixND.h>
#include <VectorND.h>
#include <Vector3D.h>
#include <Matrix3D.h>

namespace OpenSees {


static inline Matrix3D
Transpose3(const Matrix3D& A) noexcept
{
  Matrix3D AT{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      AT(i,j) = A(j,i);
  return AT;
}

template<int nsr>
static inline void
RotationOperator(MatrixND<nsr,nsr>& A, const Matrix3D& R)
{
  A.zero();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      A(i,j)     = R(i,j);
      A(i+3,j+3) = R(i,j);
    }
  }
  for (int i = 6; i < nsr; ++i)
    A(i,i) = 1.0;
}

template<int n>
static inline Vector3D
SubVector3(const VectorND<n>& v, int first) noexcept
{
  return {v[first], v[first+1], v[first+2]};
}

template<int nsr>
static inline VectorND<nsr>
ApplyOperator(const MatrixND<nsr,nsr>& A, const VectorND<nsr>& x)
{
  VectorND<nsr> y{0.0};
  for (int i = 0; i < nsr; ++i)
    for (int j = 0; j < nsr; ++j)
      y[i] += A(i,j)*x[j];
  return y;
}

template<int nsr>
static inline MatrixND<nsr,nsr>
PushTangent(const MatrixND<nsr,nsr>& A, const MatrixND<nsr,nsr>& Ks)
{
  MatrixND<nsr,nsr> ks{};
  for (int i = 0; i < nsr; ++i)
    for (int j = 0; j < nsr; ++j)
      for (int a = 0; a < nsr; ++a)
        for (int b = 0; b < nsr; ++b)
          ks(i,j) += A(i,a)*Ks(a,b)*A(j,b);
  return ks;
}


template<int nr, int nc>
static inline void
AddMatrixTriple(MatrixND<nc,nc>& Kblk,
                const MatrixND<nr,nc>& A,
                const MatrixND<nr,nr>& M,
                const MatrixND<nr,nc>& B,
                double scale)
{
  for (int i = 0; i < nc; ++i) {
    for (int j = 0; j < nc; ++j) {
      double kij = 0.0;
      for (int a = 0; a < nr; ++a) {
        const double Aai = A(a,i);
        if (Aai == 0.0)
          continue;

        double Mbj = 0.0;
        for (int b = 0; b < nr; ++b)
          Mbj += M(a,b)*B(b,j);

        kij += Aai*Mbj;
      }
      Kblk(i,j) += scale*kij;
    }
  }
}

//
// Logarithmic (Ibrahimbegovic)
//

template<std::size_t nen, int nwm>
static void
B_log(MatrixND<6+2*nwm,6+nwm>& B,
      double shape[2][nen],
      const Vector3D& dx,
      const Vector3D& th,
      const Vector3D& dth,
      int n)
{
  B.zero();

  for (int i = 0; i < 3; ++i)
    B(i,i) = shape[1][n];

  for (int i = 0; i < nwm; ++i) {
    B(6+i,     6+i) = shape[1][n];
    B(6+nwm+i, 6+i) = shape[0][n];
  }

  const Matrix3D T  = TExpSO3(th);
  const Matrix3D dR = ExpSO3(th);
  const Matrix3D Xi = dTanSO3(th, dth, 'L');

  B.assemble(Hat(dx)*T, 0, 3, shape[0][n]);
  B.assemble(T,         3, 3, shape[1][n]);
  B.assemble(dR*Xi,     3, 3, shape[0][n]);
}

template<std::size_t nen, int nwm>
static void
G_log(MatrixND<6+nwm,6+nwm>& G,
      const VectorND<6+2*nwm>& s,
      const VectorND<6+2*nwm>& S,
      const Vector3D& dx,
      double shape[2][nen],
      int i,
      int j,
      const Vector3D& th,
      const Vector3D& dth)
{
  G.zero();

  const Vector3D nsp = SubVector3(s, 0);
  const Vector3D msp = SubVector3(s, 3);
  const Vector3D Mmt = SubVector3(S, 3);

  const Matrix3D sn = Hat(nsp);
  const Matrix3D sm = Hat(msp);

  const Matrix3D T   = TExpSO3(th);
  const Matrix3D TT  = T.transpose();
  const Matrix3D snT = sn*T;
  const Matrix3D smT = sm*T;

  G.assemble(snT,   0, 3, -shape[1][i]*shape[0][j]);
  G.assemble(T^sn, 3, 0,  shape[1][j]*shape[0][i]);

  G.assemble(TT*Hat(dx)*sn*T,            3, 3,  shape[0][i]*shape[0][j]);
  G.assemble(TT*smT,                     3, 3, -shape[1][i]*shape[0][j]);
  G.assemble(dTanSO3(th, sn*dx, 'L'),    3, 3,  shape[0][i]*shape[0][j]);
  G.assemble(dTanSO3(th, msp,   'L'),    3, 3,  shape[1][i]*shape[0][j]);
  G.assemble(Transpose3(dTanSO3(th, Mmt, 'R')),
                                      3, 3,  shape[0][i]*shape[1][j]);
  G.assemble(ddTanSO3(th, dth, Mmt), 3, 3,  shape[0][i]*shape[0][j]);
}

} // namespace OpenSees
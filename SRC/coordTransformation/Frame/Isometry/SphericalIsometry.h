//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//----------------------------------------------------------------------------//
//
// Please cite the following resource in any derivative works:
//
// [1] Perez, C.M., and Filippou F.C.. (2024) 
//     "On Nonlinear Geometric Transformations of Finite Elements" 
//     Int. J. Numer. Meth. Engrg.; https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//
//
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
#pragma once
#include <cmath>
#include <array>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <MatrixND.h>
#include "EuclidIsometry.h"

namespace OpenSees {

template <int nn>
class SphericalIsometry : public Isometry<nn> {
public:
  SphericalIsometry(const Vector3D& vecxz);

  int initialize(std::array<Node*,nn>& nodes) final;
  int update(std::array<Node*,nn>& nodes) noexcept final;
  int update(const Matrix3D& RI,
             const Matrix3D& RJ,
             const Vector3D& dx, 
             std::array<Node*,nn>& nodes) final {return -1;}

  double getLength() const final { return Ln; }
  const Matrix3D& getRotation() const final;
  Vector3D getPosition() final;
  Vector3D getLocation() final;
  Matrix3D getInitialRotation() const final;
  Vector3D getPositionVariation(int ndf, double* du) final;
  Matrix3D getRotationDelta() final;
  MatrixND<3,6> getRotationGradient(int node) final;
  MatrixND<6*nn,6*nn> getHessian(const VectorND<6>& pw) final;

private:
  double L;
  double Ln;
  enum {pres, init};
  int m_I,m_J;
  // std::array<Node*,nn> nodes;
  Vector3D vz;   // vector in the x-z plane
  Vector3D dX;   // deformed length vector
  Vector3D c[2]; // current position of the center
  Matrix3D R[2]; // rotation matrices at the current and previous time step
  Vector3D Xc;   // center of the isometry
  Vector3D theta;
  double cg;
  double angle;
};


template <int nn>
SphericalIsometry<nn>::SphericalIsometry(const Vector3D& vecxz)
: vz(vecxz), Xc{},
  c{}, R{}, dX{}, Ln(0.0), L(0.0), 
  cg{0.25},
  m_I((nn+1)/2-1), m_J((nn+2)/2-1)
{

}

template <int nn>
int
SphericalIsometry<nn>::initialize(std::array<Node*,nn>& nodes) 
{
  for (int i=0; i<nn; i++)
    nodes[i]->getTrialRotation();

  const Vector &XI = nodes[m_I]->getCrds();
  const Vector &XJ = nodes[m_J]->getCrds();

  for (int i=0; i<3; i++)
    dX[i] = XJ[i] - XI[i];

  L = dX.norm();
  Ln = L;
  Vector3D e1 = dX/L;

  //
  Vector3D e2 = vz.cross(e1);

  const double ynorm = e2.norm();

  if (ynorm == 0.0)
    return -1;

  e2 /= ynorm;

  Vector3D e3 = e1.cross(e2);

  e2 = e3.cross(e1);

  for (int i = 0; i < 3; i++) {
    R[init](i,0) = e1[i];
    R[init](i,1) = e2[i];
    R[init](i,2) = e3[i];
  }
  R[pres].zero();
  R[pres].addDiagonal(1.0);

  //
  //
  
  const Vector& XC = nodes[m_I]->getCrds();
  Xc = Vector3D {XC[0], XC[1], XC[2]};
  // Cbar
  // c[init] = R[init]^Xc;
  return 0;
}

template <int nn>
const Matrix3D&
SphericalIsometry<nn>::getRotation() const {
  return R[pres];
}

template <int nn>
Matrix3D 
SphericalIsometry<nn>::getInitialRotation() const {
  return R[init];
}

template <int nn>
Vector3D
SphericalIsometry<nn>::getPosition() {
  // Return Delta c
  return {};
  // Vector3D Dc =  c[pres] - (R[init]^Xc) ; // (R[pres]^c[init]);
  // return Dc;
}

template <int nn>
Vector3D
SphericalIsometry<nn>::getLocation() {
  return c[pres];
}

template <int nn>
Matrix3D
SphericalIsometry<nn>::getRotationDelta() {
  return R[pres] - R[init];
}

template <int nn>
int
SphericalIsometry<nn>::update(std::array<Node*,nn>& nodes) 
noexcept
{
  Matrix3D RI = MatrixFromVersor(nodes[m_I]->getTrialRotation());
  Matrix3D RJ = MatrixFromVersor(nodes[m_J]->getTrialRotation());
  theta = LogSO3(RI^RJ);
  angle = theta.norm();
  R[pres] = RI*ExpSO3(theta*0.5);
  if (angle > 1e-14)
    cg = std::tan(angle*0.25)/angle;
  else
    cg = 0.25;

  return 0;
}


template <int nn>
Vector3D
SphericalIsometry<nn>::getPositionVariation(int ndf, double* du) 
{
  return Vector3D{};
}

template <int nn>
MatrixND<3,6>
SphericalIsometry<nn>::getRotationGradient(int node) {
  MatrixND<3,6> Gb{};
  Matrix3D WR{};
  if (node == m_I) {
    WR.addDiagonal(0.5);
    WR.addSpin(theta,  0.5*cg);
  }
  else if (node == m_J) {
    WR.addDiagonal(0.5);
    WR.addSpin(theta, -0.5*cg);
  }
  Gb.template insert<0,3>(WR, 1.0);
  return Gb;
}

// MatrixND<3,6>
// SphericalIsometry<nn>::getTranslationGradient(int node) final {
//   return {};
// }

template <int nn>
MatrixND<6*nn,6*nn>
SphericalIsometry<nn>::getHessian(const VectorND<6>& pw) {
  MatrixND<6*nn,6*nn> H{};
  const Matrix3D M = Hat(Vector3D{pw[3], pw[4], pw[5]});
  const double cm = (angle > 1e-14) ? 
                    (0.125/std::pow(std::cos(0.25*angle), 2))
                   : 0.125;
  for (int i=0; i<nn; i++) {
    for (int j=0; j<nn; j++) {
      if (i == j && (i == m_I || i == m_J))
        H.assemble(M, 6*i+3, 6*j+3, -cm);
      // else if ((i == m_I) && (j == m_J))
      //   H.assemble(M, 6*i+3, 6*j+3,  cm);
      else if ((i == m_I && j == m_J) ||
              (i == m_J && j == m_I)) {
        H.assemble(M, 6*i+3, 6*j+3,  cm);
      }
    }
  }
  return H;
}

} // namespace OpenSees

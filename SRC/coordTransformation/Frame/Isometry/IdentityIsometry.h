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
// References:
//
// [2] Filippou, F.C. (1998)
//     "FEDEASLab: Finite Elements for Design Evaluation and Analysis of Structures"
//
//
// Written: Claudio M. Perez, 
//          Filip C. Filippou
//          University of California, Berkeley
//
// Developed with FEDEASLab [2].
//
#pragma once
#include <array>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <MatrixND.h>
#include "EuclidIsometry.h"

namespace OpenSees {

template <int nn>
class IdentityIsometry : public Isometry<nn> {
public:
  IdentityIsometry(const Vector3D& vecxz);

  int initialize(std::array<Node*,nn>& nodes) final;
  int update(std::array<Node*,nn>& nodes) noexcept final;
  int update(const Matrix3D& RI,
             const Matrix3D& RJ,
             const Vector3D& dx, 
             std::array<Node*,nn>& nodes) final {return -1;}

  double getLength() const final { return L; }
  const Matrix3D& getRotation() const final;
  Vector3D getPosition() final;
  Vector3D getLocation() final;
  Matrix3D getInitialRotation() const final;
  Vector3D getPositionVariation(int ndf, double* du) final;
  Matrix3D getRotationDelta() final;
  MatrixND<3,6> getRotationGradient(int node) final;

  int
  setOffsets(std::array<Vector3D, nn>* offsets) override {
    this->offsets = offsets;
    return 0;
  }

private:
  double L;
  int m_I,m_J;
  Vector3D vz;   // vector in the x-z plane
  Vector3D dX;   // deformed length vector
  Vector3D c; // current position of the center
  Matrix3D R; // rotation matrices at the current and previous time step
  Vector3D Xc;   // center of the isometry
  Vector3D theta;
  double angle;
  std::array<Vector3D, nn>* offsets = nullptr; // offsets
};


template <int nn>
IdentityIsometry<nn>::IdentityIsometry(const Vector3D& vecxz)
: vz(vecxz), Xc{},
  c{}, R{}, dX{}, L(0.0), 
  m_I((nn+1)/2-1), m_J((nn+2)/2-1)
{

}

template <int nn>
int
IdentityIsometry<nn>::initialize(std::array<Node*,nn>& nodes) 
{
  for (int i=0; i<nn; i++)
    nodes[i]->getTrialRotation();

  const Vector &XI = nodes[m_I]->getCrds();
  const Vector &XJ = nodes[m_J]->getCrds();

  for (int i=0; i<3; i++)
    dX[i] = XJ[i] - XI[i];

  if (offsets != nullptr) [[unlikely]] {
    dX.addVector(1.0, (*offsets)[nn-1],  1.0);
    dX.addVector(1.0, (*offsets)[   0], -1.0);
  }

  L = dX.norm();
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
    R(i,0) = e1[i];
    R(i,1) = e2[i];
    R(i,2) = e3[i];
  }

  //
  //
  //
  // const Vector& XC = nodes[m_I]->getCrds();
  // Xc = Vector3D {XC[0], XC[1], XC[2]};
  // // Cbar
  // c[init] = R[init]^Xc;
  return 0;
}

template <int nn>
const Matrix3D&
IdentityIsometry<nn>::getRotation() const 
{
  return R;//Eye3;
}

template <int nn>
Matrix3D 
IdentityIsometry<nn>::getInitialRotation() const
{
  return R;
}

template <int nn>
Vector3D
IdentityIsometry<nn>::getPosition() 
{
  // Return Delta c
  return {};
}

template <int nn>
Vector3D
IdentityIsometry<nn>::getLocation() 
{
  return Vector3D{};
}

template <int nn>
Matrix3D
IdentityIsometry<nn>::getRotationDelta() 
{
  return Matrix3D{};//Eye3 - R; //
}

template <int nn>
int
IdentityIsometry<nn>::update(std::array<Node*,nn>& nodes) 
noexcept
{
  return 0;
}


template <int nn>
Vector3D
IdentityIsometry<nn>::getPositionVariation(int ndf, double* du) 
{
  return Vector3D{};
}

template <int nn>
MatrixND<3,6>
IdentityIsometry<nn>::getRotationGradient(int node) 
{
  return MatrixND<3,6> {};
}

} // namespace OpenSees

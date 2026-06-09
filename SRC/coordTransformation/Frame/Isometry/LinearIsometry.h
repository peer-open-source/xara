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
#pragma once
#include <Vector3D.h>
#include <Matrix3D.h>
#include <MatrixND.h>
#include "EuclidIsometry.h"

class Node;

namespace OpenSees {

template <int nn>
class LinearIsometry : public Isometry<nn>
{
public:
  LinearIsometry(const Vector3D& vecxz)
  : vz(vecxz), c{}, R{}, dX{}, Xc{}, offsets{nullptr}
  {
  }

  int
  setOffsets(std::array<Vector3D, nn>* offsets) final {
    this->offsets = offsets;
    return 0;
  }


  virtual int 
  initialize(std::array<Node*,nn>& nodes) final
  {
    for (int i=0; i<nn; i++)
      nodes[i]->getTrialRotation();

    const Vector &XI = nodes[   0]->getCrds();
    const Vector &XJ = nodes[nn-1]->getCrds();

    for (int i=0; i<3; i++)
      dX[i] = XJ[i] - XI[i];

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
      R[init](i,0) = e1[i];
      R[init](i,1) = e2[i];
      R[init](i,2) = e3[i];
    }

    const Vector& XC = nodes[ic]->getCrds();
    Xc = Vector3D {XC[0], XC[1], XC[2]};
    // Cbar
    c[init] = R[init]^Xc;

    return 0;
  }

  virtual
  int
  update(std::array<Node*,nn>& nodes) noexcept final 
  {
    return 0;
  }

  int
  update(const Matrix3D& RI, const Matrix3D& RJ, const Vector3D& dx, std::array<Node*,nn>& nodes) 
  noexcept final 
  {
    return 0;
  }


  MatrixND<3,6> 
  getRotationGradient(int node) final {
    MatrixND<3,6> Gb{};

    constexpr Matrix3D ix = Hat(Vector3D {1, 0, 0});

    if (node == 0) {
      Gb.template insert<0,0>( ix, -1.0/L);
      Gb(0,3) =   0.5;
    }
    else if (node == nn-1) {
      Gb.template insert<0,0>( ix,  1.0/L);
      Gb(0,3) =   0.5;
    }
    return Gb;
  }


  virtual Vector3D
  getRotationVariation(int ndf, double* du) final {
    // psi_r = omega
    return Vector3D {};
  }

  virtual Matrix3D
  getRotationSensitivity(std::array<Node*,nn> nodes) final {
    Matrix3D dR{};
    return dR;
  }

  double
  getLength() const override {
    return L;
  }

  const Matrix3D&
  getRotation() const final {
    return R[init];
  }

  Matrix3D
  getInitialRotation() const {
    return R[init];
  }

  virtual Matrix3D 
  getRotationDelta() final {
    return Matrix3D{};
  }

  Vector3D
  getLocation() {
    return c[init];
  }

  Vector3D
  getPosition() override {
    // Return Delta c
    return  c[init] - (R[init]^Xc) ;
  }

  Vector3D 
  getPositionVariation(int ndf, double* du) override {
    return Vector3D {};
  }

protected:
  constexpr static int ic = 0; // std::floor(0.5*(nn+1));
  enum {pres, init};
  double L;
  Vector3D vz, dX, Xc;
  Matrix3D R[2];
  Vector3D c[2];
  std::array<Vector3D, nn>* offsets = nullptr; // offsets
};
} // namespace OpenSees
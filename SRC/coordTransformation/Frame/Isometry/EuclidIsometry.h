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
// Written: Claudio M. Perez
//
// [1] Perez, C.M., and Filippou F.C. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024;
//     https://doi.org/10.1002/nme.7506
//
#pragma once
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <Versor.h>
#include <VectorND.h>
#include <Vector3D.h>
#include <MatrixND.h>
#include <Matrix3D.h>

namespace OpenSees {

template <int nn>
class Isometry
{
public:
  virtual int initialize(std::array<Node*,nn>& nodes) =0;
  // Used by EuclidFrameTransf
  virtual int update(std::array<Node*,nn>& nodes) noexcept =0;
  // Used by SouzaFrameTransf
  virtual int update(const Matrix3D& RI,
                     const Matrix3D& RJ,
                     const Vector3D& dx, 
                     std::array<Node*,nn>& nodes) =0;

  virtual int 
  setOffsets(std::array<Vector3D, nn>* offsets) {
    return -1;
  }

  virtual double    getLength() const =0;

  virtual const Matrix3D&  getRotation() const =0;
  virtual       Vector3D   getPosition() =0;
  virtual       Vector3D   getLocation() =0;
  virtual Matrix3D  getInitialRotation() const =0;

  virtual Vector3D  getPositionVariation(int ndf, double* du) =0; 
  virtual MatrixND<6*nn,6*nn> getRotationJacobian(const VectorND<6*nn>&pl) {
    // from obsolete interface, deprecated.
    return MatrixND<6*nn,6*nn> {};
  }
  virtual Matrix3D  getRotationDelta() =0;

  // compute derivative of the rotation matrix with respect to
  // nodal displacements and rotations.
  virtual MatrixND<3,6> getRotationGradient(int node) =0;
  virtual MatrixND<3,6> getTranslationGradient(int node) const {
    return MatrixND<3,6>{};
  }

  MatrixND<6,6>
  getIsometryGradient(int node) //const
  {
    MatrixND<6,6> G{};
    G.template insert<0,0>(getTranslationGradient(node));
    G.template insert<3,0>(getRotationGradient(node));
    return G;
  }

  virtual Vector3D
  getRotationVariation(int ndf, double* du) {
    // psi_r = omega
    Vector3D w{};
    for (int i=0; i<nn; i++) {
      auto Wi = this->getRotationGradient(i);
      for (int j=0; j<3; j++)
        for (int k=0; k<6; k++)
          w[j] += Wi(j,k) * du[ndf*i + k];
    }
    return w;
  }

  virtual MatrixND<6*nn,6*nn>
  getHessian(const VectorND<6>& pw) 
  {
    return MatrixND<6*nn,6*nn>{};
  }

  // compute derivative of the rotation matrix with respect to 
  // reference nodal coordinates.
  // each element of ix is 1-3 for x,y,z sensitivity, or 0 for no sensitivity.
  virtual Matrix3D
  getRotationSensitivity(std::array<Node*,nn> nodes) {
    Matrix3D dR{};
    return dR;
  }


  template <int ndf>
  inline MatrixND<nn*ndf,nn*ndf>
  getProjection(FrameTransform<nn,ndf>& t)
  {
    MatrixND<nn*ndf,nn*ndf> A{};
    A.addDiagonal(1.0);

    for (int a = 0; a < nn; ++a) {
      const Matrix3D Xa = Hat(t.getNodeLocation(a));

      for (int b = 0; b < nn; ++b) {
        MatrixND<3,ndf> Gv{};
        MatrixND<3,ndf> Gw{};

        Gv.template insert<0,0>(this->getTranslationGradient(b), 1.0);
        Gw.template insert<0,0>(this->getRotationGradient(b),    1.0);

        A.assemble(Gv,    a*ndf,   b*ndf, -1.0);
        A.assemble(Xa*Gw, a*ndf,   b*ndf,  1.0);
        A.assemble(Gw,    a*ndf+3, b*ndf, -1.0);
      }
    }

    return A;
  }
};


template <int nn>
class AlignedIsometry : public Isometry<nn>
{
public:

  AlignedIsometry(const Vector3D& vecxz)
  : vz(vecxz), Xc{}, c{}, R{}
  {
  }

  int
  setOffsets(std::array<Vector3D, nn>* offsets) override {
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

    if (offsets != nullptr) [[unlikely]] {
      dX.addVector(1.0, (*offsets)[nn-1],  1.0);
      dX.addVector(1.0, (*offsets)[   0], -1.0);
    }

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

    const Vector& XC = nodes[ic]->getCrds();
    Xc = Vector3D {XC[0], XC[1], XC[2]};
    if (offsets != nullptr) [[unlikely]] {
      Xc.addVector(1.0, (*offsets)[ic],  1.0);
    }
    // Cbar
    c[init] = R[init]^Xc;

    return this->update(nodes);
  }


  virtual
  int update(std::array<Node*,nn>& nodes) noexcept final 
  {
    const Matrix3D RI = MatrixFromVersor(nodes[0]->getTrialRotation());
    const Matrix3D RJ = MatrixFromVersor(nodes[nn-1]->getTrialRotation());

    Vector3D dx = dX;
    //
    // Update position
    //
    {
      const Vector& uI = nodes[   0]->getTrialDisp();
      const Vector& uJ = nodes[nn-1]->getTrialDisp();
      for (int k = 0; k < 3; k++)
        dx[k] += uJ(k) - uI(k);

      if (offsets != nullptr) [[unlikely]] {
        dx.addVector(1.0, (*offsets)[   0],  1.0);
        dx.addVector(1.0, RI*((*offsets)[0]), -1.0);
        dx.addVector(1.0, (*offsets)[nn-1], -1.0);
        dx.addVector(1.0, RJ*((*offsets)[nn-1]), 1.0);
      }
    }

    return this->update(RI, RJ, dx, nodes);
  }

  
  int
  update(const Matrix3D& RI, const Matrix3D& RJ, const Vector3D& dx, std::array<Node*,nn>& nodes) override 
  {
    // Calculate the deformed length
    Ln = dx.norm();

    if (Ln == 0.0) [[unlikely]] {
      opserr << "\nElement deformed length is 0.0\n";
      return -2;
    }

    //
    //
    //
    R[pres] = this->update_basis(RI, RJ, dx);
    
    //
    //
    //
    const Vector& UC = nodes[ic]->getTrialDisp();
    Vector3D uc {UC[0], UC[1], UC[2]};
    if (offsets != nullptr) {
      uc.addVector(1.0, (*offsets)[ic], -1.0);
      uc.addVector(1.0, nodes[ic]->getTrialRotation().rotate((*offsets)[ic]), 1.0);
    }
    // cbar
    c[pres] = R[pres]^(Xc + uc);
    return 0;
  }

  virtual Vector3D
  getRotationVariation(int ndf, double* du) final {
    // psi_r = omega
    Vector3D w{};
    for (int i=0; i<nn; i++) {
      auto Wi = this->getRotationGradient(i);
      for (int j=0; j<3; j++)
        for (int k=0; k<6; k++)
          w[j] += Wi(j,k) * du[ndf*i + k];
    }
    return w;
  }


  // MatrixND<3,6>
  // getTranslationGradient(int node) const override
  // {
  //   MatrixND<3,6> Gv{};

  //   if (node == ic) {
  //     Gv(0,0) = 1.0;
  //     Gv(1,1) = 1.0;
  //     Gv(2,2) = 1.0;
  //   }

  //   return Gv;
  // }

  double
  getLength() const override {
    return Ln;
  }

  const Matrix3D&
  getRotation() const final {
    return R[pres];
  }

  Matrix3D
  getInitialRotation() const override {
    return R[init];
  }

  Matrix3D 
  getRotationDelta() final {
    return R[pres] - R[init];
  }


  Vector3D
  getLocation() {
    return c[pres];
  }

  Vector3D
  getPosition() override {
    // Return Delta c
    Vector3D Dc =  c[pres] - (R[init]^Xc) ; // (R[pres]^c[init]);
    return Dc;
  }

  Vector3D 
  getPositionVariation(int ndf, double* du) override {
    return Vector3D {du[ndf*ic+0], du[ndf*ic+1], du[ndf*ic+2]};
  }

  virtual
  Matrix3D
  update_basis(const Matrix3D& RI, const Matrix3D& RJ, const Vector3D& dx) noexcept = 0;

protected:
  constexpr static int ic = 0; // std::floor(0.5*(nn+1));
  enum {pres, init};
  double L, Ln;
  Vector3D vz, dX, Xc;
  Matrix3D R[2];
  Vector3D c[2];
  Matrix3D dR;
  std::array<Vector3D, nn>* offsets = nullptr; // offsets
};

} // namespace OpenSees

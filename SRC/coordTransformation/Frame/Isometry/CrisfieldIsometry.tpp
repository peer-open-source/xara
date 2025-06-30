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
#include <array>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>
#include <Vector3D.h>
#include <MatrixND.h>
#include <Rotations.hpp>
#include "EuclidIsometry.h"

template <int nn>
class CrisfieldIsometry : public Isometry {
public:
  CrisfieldIsometry(std::array<Node*,nn>& nodes, const Vector3D& vecxz);
  
  int initialize() final;
  int update() final;
  double getLength() const final;
  Matrix3D getRotation() const final;
    Vector3D getPosition() final;
    Vector3D getPositionVariation(int ndf, double* du) final;
    Matrix3D getRotationDelta() final;
    MatrixND<3,6> getRotationGradient(int node) final;
    void setOffsets(std::array<Vector3D, nn>* offsets) {
        this->offsets = offsets;
    }

private:
  double L; // deformed length
  double Ln; // initial length
  std::array<Node*,nn> nodes;
  Vector3D vz;   // vector in the x-z plane
  Vector3D dX;   // deformed length vector
  Vector3D c[2]; // current position of the center
  Matrix3D R[2]; // rotation matrices at the current and previous time step
  Vector3D Xc;   // center of the isometry
  std::array<Vector3D, nn>* offsets; // rigid joint offsets
  int offset_flags; // flags for offsets
  constexpr static int ic = 0; // index of the center node
};


template <int nn>
CrisfieldIsometry<nn>::CrisfieldIsometry(std::array<Node*,nn>& nodes, const Vector3D& vecxz)
: nodes(nodes), vz(vecxz), Xc{},
  c{}, R{}, dX{}, Ln(0.0), L(0.0), 
  offsets(nullptr), offset_flags(0)
{
  for (int i = 0; i < nn; ++i) {
    nodes[i]->getTrialRotation();
  }
}

template <int nn>
int
CrisfieldIsometry<nn>::initialize()
{
  for (int i=0; i<nn; i++)
    nodes[i]->getTrialRotation();

  const Vector &XI = nodes[   0]->getCrds();
  const Vector &XJ = nodes[nn-1]->getCrds();

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

  Xc = nodes[ic]->getCrds();
  c[init] = R[init]^(Xc);

  this->update();
  return 0;
}

template <int nn>
int
CrisfieldIsometry<nn>::update()
{
  Vector3D e1 = dX;
  {
    //
    // Update state
    //
    const Vector& uI = nodes[   0]->getTrialDisp();
    const Vector& uJ = nodes[nn-1]->getTrialDisp();
    for (int k = 0; k < 3; k++) 
      e1[k] += uJ(k) - uI(k);

    if (offsets != nullptr) [[unlikely]] {
      e1.addVector(1.0, (*offsets)[   0],  1.0);
      e1.addVector(1.0, nodes[0]->getTrialRotation().rotate((*offsets)[0]), -1.0);
      e1.addVector(1.0, (*offsets)[nn-1], -1.0);
      e1.addVector(1.0, nodes[nn-1]->getTrialRotation().rotate((*offsets)[nn-1]), 1.0);
    }

    // Calculate the deformed length
    Ln = e1.norm();
    if (Ln == 0.0) [[unlikely]] {
        opserr << "\nCrisfieldIsometry: deformed length is 0.0\n";
        return -2;
    }
    e1 /= Ln;
  }

  {
    Versor q0 = VersorFromMatrix(R[init]);
    Versor qI = nodes[0]->getTrialRotation()*q0;
    Versor qJ = nodes[nn-1]->getTrialRotation()*q0;
    Vector3D gammaw = CayleyFromVersor(qJ.mult_conj(qI));

    gammaw *= 0.5;

    //  Qbar = VersorProduct(VersorFromMatrix(CaySO3(gammaw)), qI);
    Matrix3D Rbar = CaySO3(gammaw)*MatrixFromVersor(qI); // *q0);
    Vector3D v { Rbar(0,0), Rbar(1,0), Rbar(2,0) };
    double dot = v.dot(e1);
    if (std::fabs(std::fabs(dot)-1.0) < 1.0e-10) {
      R[pres] = Rbar;
    }
    else {
      v  = v.cross(e1);
      double scale = std::acos(dot)/v.norm();
      v *= scale; // ::acos(r1.dot(e1));

      R[pres] = ExpSO3(v)*Rbar;

      Vector3D r1 { R[pres](0,0), R[pres](1,0), R[pres](2,0) };
      Vector3D r2 { R[pres](0,1), R[pres](1,1), R[pres](2,1) };
      Vector3D r3 { R[pres](0,2), R[pres](1,2), R[pres](2,2) };
    }
  }
}
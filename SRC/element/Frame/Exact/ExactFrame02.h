
//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2026, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// ExactFrame02
//
// Geometrically exact 3D frame element using the logarithmic rotation
// parameterizations.
//
//===----------------------------------------------------------------------===//
#pragma once
#include <set>
#include <array>
#include <Element.h>
#include <VectorND.h>
#include <MatrixND.h>

#include <FrameSection.h>
#include <FrameTransform.h>
#include <Frame/FiniteElement.h>
#include <GroupSO3.h>
#include <Rotations.h>

class CrdTransf;

namespace OpenSees {

class FrameLoad;

template<std::size_t nen, int nwm=0>
class ExactFrame02:
  public FiniteElement<nen, 3, 6+nwm>
{
public:

  ExactFrame02(int tag,
               std::array<int,nen>& nodes,
               FrameSection *section[nen-1],
               CrdTransf& transf,
               Rotations::Parameters param=Rotations::Parameters::Incr);

  ~ExactFrame02();

  const char*
  getClassType() const final
  {
    return "ExactFrame02";
  }

  int setNodes() override;
  const Vector &getResistingForce() final;
  const Matrix &getTangentStiff() override;
  const Matrix &getMass() override;
  const Matrix &getInitialStiff() override;

  int update() final;
  int addLoad(ElementalLoad* , double scale) final;

  int revertToStart() override;
  int revertToLastCommit() override;
  int commitState() final {

    if (this->Element::commitState() != 0)
      opserr << "ExactFrame02::commitState() - failed in base class";

    past = pres;

    for (GaussPoint& point : pres) {
      if (point.material->commitState() != 0)
        return -1;
    }
    return 0;
  }

  int setParameter(const char **argv, int argc, Parameter &) final;
  int updateParameter(int parameterID, Information &) final;
  int activateParameter(int parameterID) final;
  const Vector& getResistingForceSensitivity(int grad) final;

  Response *setResponse(const char **argv, int argc, OPS_Stream &s) final;
  int getResponse(int responseID, Information &) final;

  int sendSelf(int cTag, Channel&) override;
  int recvSelf(int cTag, Channel&, FEM_ObjectBroker&) override;

  void Print(OPS_Stream& s, int flag) override;

private:
  constexpr static int
        nsr = 6+2*nwm,
        ndm = 3,
        nip = nen-1,
        ndf = 6+nwm;

  enum Respond: int {
    GlobalForce = 1,
    BasicPlasticDeformation = 4,
    LocalForce  = 2,
    BasicForce  = 7,
    BasicStiff  = 19,
  };

  static constexpr FrameStressLayout scheme = {
    FrameStress::N,
    FrameStress::Vy,
    FrameStress::Vz,
    FrameStress::T,
    FrameStress::My,
    FrameStress::Mz,
    FrameStress::Bimoment,
    FrameStress::Bishear,
  };

  struct GaussPoint {
    double point,
           weight;
    FrameSection* material;

    double   shape[2][nen];
    Matrix3D rotation;
    Vector3D curvature;
  };

  double xn[nen];
  double jxs;

  Matrix3D R0;

  std::set<FrameLoad*> frame_loads;

  std::array<GaussPoint,nip> pres;
  std::array<GaussPoint,nip> past;
  std::array<Matrix3D,nip> reference_rotation;
  std::array<Vector3D,nip> reference_curvature;

  CrdTransf* transform;
  const Rotations::Parameters parameterization;
  BeamIntegration* stencil;

  VectorND<nen*ndf> p;
  MatrixND<nen*ndf,nen*ndf> K;

  int parameterID;
};

} // namespace OpenSees
#include "ExactFrame02.tpp"

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
//===----------------------------------------------------------------------===//
//
// Claudio M. Perez
//
#pragma once
#include <set>
#include <array>
#include <Element.h>
#include <VectorND.h>
#include <MatrixND.h>

#include <Frame/FrameMass.hpp>
#include <FrameSection.h>
#include <FrameTransform.h>
#include <Frame/FiniteElement.h>
#include <GroupSO3.h>

class CrdTransf;

namespace OpenSees {

class FrameLoad;

template<std::size_t nen, int nwm=0>
class ExactFrame3d: 
  public FiniteElement<nen, 3, 6+nwm>
{
public:
  enum Logarithm {
    None,
    LogLeft,
    LogRight
  };

  ExactFrame3d(int tag,
               std::array<int,nen>& nodes,
               FrameSection *section[nen-1], 
               CrdTransf& transf
  );

  ~ExactFrame3d();

  const char*
  getClassType() const final
  {
    return "ExactFrame3d";
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
      opserr << "ExactFrame3d::commitState () - failed in base class";

    past = pres;

    for (GaussPoint& point : pres) {
      if (point.material->commitState() != 0)
        return -1;
    }
    return 0;
  }

  // Element: Parameters
  int setParameter(const char **argv, int argc, Parameter &) final;
  int updateParameter(int parameterID, Information &) final;
  int activateParameter(int parameterID) final;
  // Element: Sensitivity
  const Vector& getResistingForceSensitivity(int grad) final;

  Response *setResponse(const char **argv, int argc, OPS_Stream &s) final;
  virtual int getResponse(int responseID, Information &) final;

  // MovableObject
  int sendSelf(int cTag, Channel&) override;
  int recvSelf(int cTag, Channel&, FEM_ObjectBroker&) override;

  // TaggedObject
  void Print(OPS_Stream& s, int flag) override;

  private:
    //
    // Constexpr
    //
    constexpr static int 
          nsr = 6+2*nwm,        // Number of section resultants
          ndm = 3,              // Dimension of the problem (3D)
          nip = nen-1,          // Number of integration points
          ndf = 6+nwm;          // Degrees of freedom per node

    enum Respond: int {
      GlobalForce = 1,
      BasicPlasticDeformation = 4,
      LocalForce  = 2,
      BasicForce  = 7,
      BasicStiff  =19,
    };
    
    // Layout of stress resultants
    static constexpr FrameStressLayout scheme = {
      FrameStress::N,
      FrameStress::Vy,
      FrameStress::Vz,
      FrameStress::T,
      FrameStress::My,
      FrameStress::Mz,
      FrameStress::Bimoment,
  //  FrameStress::By,
  //  FrameStress::Bz,
      FrameStress::Bishear,
  //  FrameStress::Qy,
  //  FrameStress::Qz
    };

    //
    //
    //
    struct GaussPoint {
      double point,
             weight;
      FrameSection* material;

      double   shape[2][nen];
      Matrix3D rotation;
      Vector3D curvature;
    };
    // Node locations in local (scalar) coordinate
    double xn[nen];
    double jxs;

    Matrix3D R0;

    std::set<FrameLoad*> frame_loads;

    std::array<GaussPoint,nip> pres;
    std::array<GaussPoint,nip> past;
    CrdTransf*       transform;
    Logarithm               logarithm;
    BeamIntegration*        stencil;

    //
    //
    VectorND<nen*ndf> p;
    MatrixND<nen*ndf,nen*ndf> K;
    // Parameters
    int parameterID;
};

} // namespace OpenSees
#include "ExactFrame3d.tpp"

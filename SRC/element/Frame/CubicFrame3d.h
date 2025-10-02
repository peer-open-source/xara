//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
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
#ifndef CubicFrame3d_h
#define CubicFrame3d_h

#include <array>
#include <vector>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <VectorND.h>
#include <FrameSection.h>

#include <ID.h>

class Node;
class CrdTransf;
class BeamIntegration;
class Response;
class BasicFrame3d;

namespace OpenSees {

template <bool shear, int nwm>
class CubicFrame3d : public Element {
public:
  CubicFrame3d(int tag, 
               std::array<int, 2>&,
               std::vector<FrameSection*>&, 
               BeamIntegration&,
               CrdTransf&, 
               double rho);
  CubicFrame3d();
  ~CubicFrame3d();


  const char*
  getClassType() const override
  {
    return "CubicFrame3d";
  }

  int getNumExternalNodes() const override;
  const ID& getExternalNodes() override;
  Node** getNodePtrs() override;

  int getNumDOF() override;
  void setDomain(Domain*) override;

  int update() override;
  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;

  // public methods to obtain stiffness, mass, damping and residual
  const Matrix& getTangentStiff() final;
  const Matrix& getInitialStiff() final;
  const Matrix& getMass() final;

  void zeroLoad() final;
  int addLoad(ElementalLoad*, double loadFactor) final;
  int addInertiaLoadToUnbalance(const Vector& accel) final;

  const Vector& getResistingForce() final;
  const Vector& getResistingForceIncInertia() final;


  Response* setResponse(const char** argv, int argc, OPS_Stream& s) final;
  int getResponse(int responseID, Information& eleInfo) final;

  // Element: Parameters
  int setParameter(const char** argv, int argc, Parameter&) final;
  int updateParameter(int parameterID, Information&) final;
  int activateParameter(int parameterID) final;

  // Element: Sensitivity
  const Vector& getResistingForceSensitivity(int gradNumber) final;
  const Matrix& getInitialStiffSensitivity(int gradNumber) final;
  const Matrix& getMassSensitivity(int gradNumber) final;
  int commitSensitivity(int gradNumber, int numGrads) final;

  // MovableObject
  int sendSelf(int commitTag, Channel&) override;
  int recvSelf(int commitTag, Channel&, FEM_ObjectBroker&) override;

  // TaggedObject
  void Print(OPS_Stream& s, int flag) final;

protected:

private:
  void getBasicStiff(Matrix& kb, int initial = 0);
private:
  constexpr static int 
        nsr = shear? 6 : 4,   // number of section resultants
        ndm = 3,              // dimension of the problem (3D)
        nq  = 6+nwm*2,        // number of element dof's in the basic system
        NDF = 6+nwm,          //
        NEN = 2,              // number of element nodes
        maxNumSections = 20;


  int numSections;
  FrameSection** theSections;       // the materials
  CrdTransf* theCoordTransf;
  BeamIntegration* beamInt;

  double xi[maxNumSections];
  double wt[maxNumSections];
  double phizs[maxNumSections]; // Shear term 12EIz/(GA L^2)
  double phiys[maxNumSections]; // Shear term 12EIy/(GA L^2)

  ID connectedExternalNodes;    // Tags of nodes

  Node* theNodes[NEN];


  Vector Q;       // Applied nodal loads
  VectorND<nq> q;

  double density;             // Mass density per unit length
  double total_mass,
         twist_mass;
  int    mass_flag;
  bool   use_density;
  bool   mass_initialized;

  BasicFrame3d* loads;

  int parameterID;

  // Matrix K; // Element stiffness, damping, and mass Matrix
  MatrixND<NEN*NDF,NEN*NDF> K;
  VectorND<NDF*NEN> P; // Element resisting force vector
  Matrix K_wrap;
  Vector P_wrap;

  static constexpr FrameStressLayout scheme = {
      FrameStress::N, 
      FrameStress::T, 
      FrameStress::My, 
      FrameStress::Mz,
      FrameStress::Vy, 
      FrameStress::Vz,
  };
};
}
#include <CubicFrame3d.tpp>
#endif

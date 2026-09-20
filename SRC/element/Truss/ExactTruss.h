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
//
// [1] Crisfield, M.A. "Nonlinear Finite Element Analysis of Solids and Structures", 
//     Vol. 1, 1991, J.T. Wiley.
//
#pragma once
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <Matrix3D.h>
#include <FrameSection.h>
#include <Node.h>


class ExactTruss : public Element {
public:
  ExactTruss(int tag, int dim, 
             int Nd1, int Nd2, 
             FrameSection&,
             double rho,
             int doRayleighDamping, 
             int cMass, 
             int strain);

  ~ExactTruss();

  const char*
  getClassType() const
  {
    return "ExactTruss";
  }

  // public methods to obtain information about dof & connectivity
  int getNumExternalNodes() const override;
  const ID& getExternalNodes();
  Node** getNodePtrs();

  int getNumDOF() override;
  void setDomain(Domain* theDomain);

  // public methods to set the state of the element
  int commitState() final;
  int revertToLastCommit() final;
  int revertToStart() final;
  int update() final;

  // public methods to obtain stiffness, mass, damping and residual information
  const Matrix& getTangentStiff() final;
  const Matrix& getInitialStiff() final;
  const Matrix& getDamp() final;
  const Matrix& getMass() final;

  void zeroLoad() final;
  int addLoad(ElementalLoad* theLoad, double loadFactor) final;

  const Vector& getResistingForce() final;
  const Vector& getResistingForceIncInertia() final;

  // Sensitivity 
  int setParameter(const char** argv, int argc, Parameter& param) final;
  int updateParameter(int parameterID, Information& info) final;
  int activateParameter(int param) final;

  // public methods for element output
  void Print(OPS_Stream& s, int flag) final;

  Response* setResponse(const char** argv, int argc, OPS_Stream& s) final;
  int getResponse(int responseID, Information& );


private:
  enum class StrainType {
    Linear,
    Green,
    Corotational,
    Logarithmic
  } strain_type = StrainType::Corotational;

  double computeCurrentStrain();

  double locate(Vector3D&dX, Vector3D&dx, Vector3D&du) {
    dX.zero();
    dx.zero();
    du.zero();
    const Vector& uj = theNodes[1]->getTrialDisp();
    const Vector& ui = theNodes[0]->getTrialDisp();
    const Vector& Xj = theNodes[1]->getCrds();
    const Vector& Xi = theNodes[0]->getCrds();
    for (int i = 0; i < numDIM; i++) {
      du[i] = uj(i) - ui(i);
      dX[i] = Xj(i) - Xi(i);
      dx[i] = dX[i] + du[i];
    }

    return dx.norm();
  }

  double interpolate(Vector3D&B, double&Ln) {
    Vector3D dX, dx, du;
    Ln = this->locate(dX, dx, du);

    switch (strain_type) {
      case StrainType::Linear:
        B = dX/Lo;
        return 1.0;
      case StrainType::Green:
        B = (dX + du)/(Lo*Lo);
        return Ln/Lo;
      case StrainType::Corotational:
        B = (dX + du)/(Lo*Ln);
        return 1.0;
      case StrainType::Logarithmic:
        B = (Lo/Ln)*(dX + du)/(Lo*Ln);
        return 1.0;
    }
  }

  // Layout of stress resultants
  static constexpr FrameStressLayout section_layout = {
    FrameStress::N,
  };

  // private attributes
  FrameSection* theSection; // material (section) response
  ID connectedExternalNodes;           // contains the tags of the end nodes
  int numDOF;                          // number of dof for element
  int numDIM;                          // number of dimensions

  double Lo;             // initial length of truss
  double rho;            // mass density per unit length
  int doRayleighDamping; // flag to include Rayleigh damping
  int cMass;             // consistent mass flag

  Node* theNodes[2];

  Matrix3D R; // Rotation matrix
  int parameterID;

  Vector* theLoad;   // pointer to the load vector P
  Matrix* theMatrix; // pointer to objects matrix (a class wide Matrix)
  Vector* theVector; // pointer to objects vector (a class wide Vector)

  static Matrix M2;
  static Matrix M4;
  static Matrix M6;
  static Matrix M12;

  static Vector V2;
  static Vector V4;
  static Vector V6;
  static Vector V12;
};

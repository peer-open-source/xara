//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
//
#pragma once

#define ND_TAG_J2BeamThread3d 92516

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Information.h>
#include <Parameter.h>
#include <Matrix3D.h>
#include <VectorND.h>
#include <Vector3D.h>
#include <domain/DomainStatus.h>
class Response;

namespace OpenSees {

class J2BeamThread3d : public NDMaterial {
public:
  J2BeamThread3d(int tag, 
                 double E, double nu, 
                 double sigY, double Hi, double Hk,
                 double density
                );
  J2BeamThread3d();
  ~J2BeamThread3d();

  const char* getClassType() const override {return "J2BeamThread3d";}

  int setTrialStrain(const Vector&) override;
  int setTrialStrainIncr(const Vector& v) override;
  const Matrix& getTangent() override;
  const Matrix& getInitialTangent() override;
  const Vector& getStress() override;
  const Vector& getStrain() override;
  bool threadSafe() const override {return true;}
  double getRho() override { return density; }


  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;

  NDMaterial* getCopy() override;
  NDMaterial* getCopy(const char* type) override;
  const char* getType() const override;
  int getOrder() const;

  // Response *setResponse (const char **argv, int argc, OPS_Stream &);
  // int getResponse (int responseID, Information &);

  int sendSelf(int commitTag, Channel& ) override;
  int recvSelf(int commitTag, Channel& , FEM_ObjectBroker&) override;

  void Print(OPS_Stream& s, int flag) override;

  int setParameter(const char** argv, int argc, Parameter&) override;
  int updateParameter(int parameterID, Information&) override;
  int activateParameter(int paramID);

  const Vector& getStressSensitivity(int gradIndex, bool conditional);
  int commitSensitivity(const Vector& depsdh, int gradIndex, int numGrads);


private:
  DomainStatus returnMap() noexcept;
  double E;
  double nu;
  double sigmaY;
  double Hiso;
  double Hkin;
  double density;

  int parameterID;
  Matrix* SHVs;

  Matrix3D D;     // tangent
  Vector3D Tepsilon;     // Trial strains
  Vector3D sigma;
  Vector s_wrapper;
  Vector e_wrapper;
  Matrix D_wrapper;

  Vector3D  xsi_n1;     // converged relative (shifted) stress
  double    q_n1;       // converged equivalent stress  (scalar)

  double alphan;
  double alphan1;

  // enum class Request : int {
  //   PlasticStrain = 1,
  // };

  double epsPn[3];
  double epsPn1[3];
  // Vector3D xsi;

  double dg_n1;
};
} // namespace OpenSees

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
//===----------------------------------------------------------------------===//
//
// Description: Elastic isotropic model where stress components 22, 33, and 23
// are condensed out.
//
// Written: CMP
// Created: Feb 2026
//
#ifndef ElasticIsotropicBeamThread_h
#define ElasticIsotropicBeamThread_h

#include <ElasticIsotropicMaterial.h>
#include <memory>
#include <Matrix.h>
#include <Vector.h>
#include <memory>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <ID.h>
#include <Information.h>
#include <Parameter.h>

class MaterialBuilder;

namespace OpenSees {

class ElasticIsotropicBeamThread : public NDMaterial
{
  public:
    ElasticIsotropicBeamThread(ElasticIsotropicMaterial&);
  private:
    ElasticIsotropicBeamThread(const ElasticIsotropicBeamThread &);
  public:
    ~ElasticIsotropicBeamThread();
    const char *getClassType() const override {
      return "ElasticIsotropic<BeamFiber>";
    }

    int setTrialStrain(const Vector &v) override;
    int setTrialStrainIncr (const Vector &v) override;
    const Matrix &getTangent() override;
    const Matrix &getInitialTangent() override;
    const Vector &getStress() override;
    const Vector &getStrain() override;
    bool threadSafe() const override {return true;}
    double getRho() override {return rho;}
        
    int commitState() override;
    int revertToLastCommit() override;
    int revertToStart() override;
    
    NDMaterial *getCopy() override;
    NDMaterial *getCopy (const char *type) override;
    const char *getType() const override;
    int getOrder() const override;

    int setParameter(const char **argv, int argc, Parameter &param) override;
    const Vector& getStressSensitivity(int gradIndex, bool conditional) override;

    int sendSelf(int commitTag, Channel &) override {return -1;};
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override {return -1;}

  private:
    Vector3D stress;
    Vector3D Tepsilon;		// Trial strains
    std::shared_ptr<Matrix3D> D;	// Elastic constants
    Vector retStrain;
    Vector retStress;
    Matrix retTangent;
    double rho;
    int parameterID;
    MaterialBuilder& builder;
};
}

#endif

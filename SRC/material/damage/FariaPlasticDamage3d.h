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
//
//
#pragma once
//
// Written: Thanh Do
// Created: 07/16
//
#include <NDMaterial.h>
#include <MatrixND.h>
#include <ID.h>
class Vector;
class Matrix;
class Channel;

class FariaPlasticDamage3d : public NDMaterial
{
  public:
  FariaPlasticDamage3d(int tag, 
                        double E, 
                        double nu, 
                        double ft,
                        double fc, 
                        double beta, 
                        double Ap, 
                        double An, 
                        double Bn,
                        double density);
  FariaPlasticDamage3d();
  ~FariaPlasticDamage3d();

  const char *getClassType() const override {return "FariaPlasticDamage";}

  int setTrialStrain(const Vector &v) override;
  int setTrialStrain(const Vector &v, const Vector &r);
  int setTrialStrainIncr(const Vector &v) override;
  int setTrialStrainIncr(const Vector &v, const Vector &r);
  const Matrix &getTangent();
  const Matrix &getInitialTangent();
  
  const Vector &getStress();
  const Vector &getStrain();
  
  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;

  NDMaterial *getCopy(const char *type);
  NDMaterial *getCopy() override;
  const char *getType() const override;
  int getOrder() const override;

  int sendSelf(int commitTag, Channel &) override;  
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;    
  void Print(OPS_Stream &s, int flag) override;       


  private:
    // parameters
    double E;     // elastic modulus
    double nu;    // Poisson ratio 
    double ft;    // tensile yield strength
    double Fc;    // compressive yield strength
    double beta;  // plastic deformation rate
    double Ap;    // damage parameter
    double An;    // damage parameter
    double Bn;    // damage parameter

    // current state variables
    double rp;    // positive damage threshold
    double rn;    // negative damage threshold
    double dp;    // positive damage variable
    double dn;    // negative damage variable

    OpenSees::VectorND<6> eps;   // strain
    OpenSees::VectorND<6> sig;   // stress
    OpenSees::VectorND<6> sige;  // effective stress
    OpenSees::VectorND<6> eps_p; // plastic strain

    // committed state variables
    double rpCommit; 
    double rnCommit; 
    double dpCommit; 
    double dnCommit; 

    OpenSees::VectorND<6> epsCommit;
    OpenSees::VectorND<6> sigCommit;
    OpenSees::VectorND<6> sigeCommit;  
    OpenSees::VectorND<6> eps_pCommit;

    // tangent matrices
    OpenSees::MatrixND<6,6> Ce, C, Ccommit; 
    Matrix retTangent, retInitialTangent;
    Vector retStress, retStrain;

    double density;
};


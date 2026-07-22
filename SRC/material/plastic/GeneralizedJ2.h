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
#pragma once
//
#include <NDMaterial.h>
#include <MatrixND.h>
#include <Matrix.h>
#include <Vector.h>
class Channel;

namespace OpenSees {

class GeneralizedJ2 : public NDMaterial
{
public:
  enum class HRule { LP, GP /*, NLK */ };

  GeneralizedJ2(int tag, 
                double E, 
                double nu, 
                double Fy,
                double Fs, 
                double Hiso, 
                double Hkin, 
                double delta,
                double density,
                HRule rule = HRule::GP);
  ~GeneralizedJ2();

  const char *getClassType() const override {return "GeneralizedJ2";}
  const char *getType() const override {return "ThreeDimensional";}
  bool threadSafe() const override {return true;}

  int setTrialStrain (const Vector &v) override;
  int setTrialStrainIncr (const Vector &v) override;

  const Matrix &getTangent() override;
  const Matrix &getInitialTangent() override;
  
  const Vector &getStress() override;
  const Vector &getStrain() override;
  
  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;

  NDMaterial *getCopy(const char *type) override;
  NDMaterial *getCopy() override;
  int getOrder() const override {return 6;}

  int sendSelf(int commitTag, Channel &) override;  
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;    
  void Print(OPS_Stream &, int flag) override;       

private:
  // parameters
  double E;       // elastic modulus
  double nu;      // Poisson ratio
  double Fy;      // yield stress
  double Fs;      // saturated yield stress ( = Fy + phi )
  double Hiso;    // isotropic hardening parameter (Hi)
  double Hkin;    // kinematic hardening parameter (Hk)
  double delta;   // saturation hardening parameter (Auricchio–Taylor delta)

  // hardening rule selector (LP or GP; NLK needs Hnl)
  HRule rule;

  // history (past/present)
  struct {
    VectorND<6> sig;    // Cauchy stress
    VectorND<6> eps_p;  // plastic strain
    VectorND<9> sig_b;  // back-stress (deviatoric, Voigt)
    double sigdnrm;     // || dev(s) - back-stress || (J2 metric)
    double alpha;       // isotropic hardening scalar
  } past, pres;

  // Total strain, elastic tangent
  VectorND<6>      eps;   // total strain
  MatrixND<6,6>    Ce;    // elastic stiffness
  MatrixND<6,6>    C;     // current consistent tangent

  // wrappers for OpenSees API returns
  Matrix retTangent, retInitialTangent;
  Vector retStress,  retStrain;

  double density;

private:

  // hardening rules
  struct Hardening {
    double Hi, Hk, G, delta, phi; // phi = Fs - Fy
    HRule rule;
  };
  struct HState {
    double A2;           // sig_nrm - past.sigdnrm
    double sn;           // ||trial_dev - b||
    VectorND<9> se;      // trial deviatoric stress
    VectorND<9> sb;      // back-stress
    double alpha;        // isotropic hardening scalar
  };


  int updateState(const VectorND<6> &eps);

  static int hard_LP(const Hardening &hd,
                     double yf_tr,
                     double &Dlam, double &xi_p, double &Tlam) noexcept;

  static int hard_GP(const Hardening &hd, 
                     double yf_tr, 
                     const HState &hs,
                     double &Dlam, double &xi_p, double &Tlam) noexcept;

};

} // namespace OpenSees

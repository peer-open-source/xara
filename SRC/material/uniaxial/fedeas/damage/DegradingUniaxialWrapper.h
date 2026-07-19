/* ============================================================================
 * FEDEASLab - Release 5.1, July 2020
 * Matlab Finite Elements for Design, Evaluation and Analysis of Structures
 * Professor Filip C. Filippou (filippou@berkeley.edu)
 * Department of Civil and Environmental Engineering, UC Berkeley
 * Copyright(c) 1998-2021. The Regents of the University of California. 
 * All Rights Reserved.
 * ============================================================================
 * Adapted from MATLAB source code by Claudio Perez                     07/2021
 */
#ifndef DegradingUniaxialWrapper_h
#define DegradingUniaxialWrapper_h

#include <Parsing.h>
#include <UniaxialMaterial.h>
#include <VectorND.h>
#include "evolution.h"
using namespace OpenSees;

/* Fracture data. Probably rename too */
struct DmgFrac {
  bool Activ;
  double psiF, psiU;
};

/* Damage parameters; eventually rename */
#define MAX_EVOL_PARAMS 5

struct DmgEvol {
  int nParam = 0;
  double theParam[MAX_EVOL_PARAMS];
};



struct DamageIndex {
  bool valid = false;

  double Cd0=0, 
         Cd1=0, 
         Cwc=0,
         Ccd=1,
         psi_d0=0, 
         psi_d1=0;

  enum class Type {
    None, mBeta, oBeta, Log, Weibull, Bilinear, Trilinear
  } type = Type::None;

  DmgEvol evol;  /* Damage evolution data */
  DmgFrac frac;  /* Fracture data */

  struct DmgResp update(double psi_tild) const {
    switch (type) {
      case Type::mBeta:
        return dmglib_MBeta(psi_tild, evol, frac);
      case Type::oBeta:
        return dmglib_OBeta(psi_tild, evol, frac);
      case Type::Log:
        return dmglib_Logn(psi_tild, evol, frac);
      case Type::Weibull:
        return dmglib_Wbl(psi_tild, evol, frac);
      case Type::Bilinear:
        return dmglib_Bilin(psi_tild, evol, frac);
      case Type::Trilinear:
        return dmglib_Trilin(psi_tild, evol, frac);
      case Type::None:
      default:
        return dmglib_none(psi_tild, evol, frac);
    }
  }
};


class DegradingUniaxialWrapper : public UniaxialMaterial {
public:
  
  struct Data {
    double tol = 1e-10;
    DamageIndex idx[2];
  };

  DegradingUniaxialWrapper(int tag, UniaxialMaterial &, Data data);
  ~DegradingUniaxialWrapper();

  const char *
  getClassType() const {return "UniaxialDamage";}

  int    setTrialStrain(double strain, double strainRate = 0.0);
  int    setTrialStrain(double strain, double temperature, double strainRate);
  double getStrain();
  double getStrainRate();
  double getStress();
  double getTangent();
  double getDampTangent();
  double getInitialTangent();

  int commitState();
  int revertToLastCommit();
  int revertToStart();

  UniaxialMaterial *getCopy();

  int sendSelf(int commitTag, Channel &);
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

  void Print(OPS_Stream &s, int flag);

  int setCoupling(double);
  int setParameter(const char **argv, int argc, Parameter &);
  int updateParameter(int parameterID, Information &);

  double getStressSensitivity(int gradIndex, bool conditional);
  double getStrainSensitivity(int gradIndex);
  double getInitialTangentSensitivity(int gradIndex);
  double getDampTangentSensitivity(int gradIndex);
  double getRhoSensitivity(int gradIndex);
  int commitSensitivity(double strainGradient, int gradIndex, int numGrads);
  
//  int setDamageWrapper(Tcl_Interp*, std::string);

private:
  UniaxialMaterial *theMaterial;

  Data data;

  //
  // State variables
  //
  double m_stress,
         m_tangent, 
         m_rate_tol=1e-6;

  struct UniaxialState {
  /* This struct defines the interface between
   * the public wrapper and it's external 
   * implementation. */
   double  e, ep, De, se, kt, ke;
  };


  struct DamageState {
    struct Step {
      double
            strain  = 0.0,
            base_stress[2] = {0.0, 0.0};
      VectorND<2>
            d       = {{0.0, 0.0}},
            vpEx    = {{0.0, 0.0}},
            psi     = {{0.0, 0.0}},
            psiEx   = {{0.0, 0.0}};
    } test, past;
  } dstate;

  int applyDamage(
        const Data&     data,
        DamageState& instance,
        const double inputs[3],
        double response[2]);

  // Enums used for readable array indices.
  enum {PAST, PRES};
  // enum {ax, mi, mj};
  // Rows of Cin matrix; indicates moment sign, axial sign
  enum {pp, pn, np, nn};

  // UniaxialState past,pres;
};

#endif // DegradingUniaxialWrapper_H


//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// ============================================================================
// FEDEASLab - Release 5.1, July 2020
// Matlab Finite Elements for Design, Evaluation and Analysis of Structures
// Professor Filip C. Filippou (filippou@berkeley.edu)
// Department of Civil and Environmental Engineering, UC Berkeley
// Copyright(c) 1998-2021. The Regents of the University of California. 
// All Rights Reserved.
// ============================================================================
// Adapted from MATLAB source code by Claudio Perez                     07/2021
//
#include <functional> // std::hash
#include <string>
#include <cmath>

#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Logging.h>
#include <ModelRegistry.h>
#include <VectorND.h>
using namespace OpenSees;

#include "DegradingUniaxialWrapper.h"

#define WRAPPER_CMD "UniaxialDamage"
static const std::hash<std::string> hasher;
static const int MatTag = hasher(WRAPPER_CMD);


DegradingUniaxialWrapper::DegradingUniaxialWrapper(int tag,
                                                   UniaxialMaterial &material,
                                                   Data data)
  : UniaxialMaterial(tag, MatTag)
  , theMaterial(nullptr)
  , m_stress(0.0)
  , data(data)
{
  theMaterial = material.getCopy();
  m_tangent = theMaterial->getInitialTangent();
}


DegradingUniaxialWrapper::~DegradingUniaxialWrapper()
{
  if (theMaterial)
    delete theMaterial;
}


int
DegradingUniaxialWrapper::setTrialStrain(double strain, double temp,
                                         double strainRate)
{

  theMaterial->setTrialStrain(strain, temp, strainRate);

  if (true) { //  && abs(strain_incr) > m_rate_tol){

    double inputs[] = {
      strain,
      theMaterial->getStress(),
      theMaterial->getTangent(),
    };

    double outputs[2];

    this->applyDamage(data, dstate, inputs, outputs);

    this->m_stress  = outputs[0];
    this->m_tangent = outputs[1];

  } else {
    this->m_stress = theMaterial->getStress();
    this->m_tangent = theMaterial->getTangent();
  }
  return 0;
}

int
DegradingUniaxialWrapper::setTrialStrain(double trialStrain, double strainRate)
{
  return this->setTrialStrain(trialStrain, 0.0, strainRate);
}

double
DegradingUniaxialWrapper::getTangent()
{
  return this->m_tangent;
}

double
DegradingUniaxialWrapper::getStress()
{
  return m_stress;
}

double
DegradingUniaxialWrapper::getInitialTangent()
{
  return theMaterial->getInitialTangent();
}

double
DegradingUniaxialWrapper::getDampTangent()
{
  return theMaterial->getDampTangent();
}

double
DegradingUniaxialWrapper::getStrain()
{
  return theMaterial->getStrain();
}

double
DegradingUniaxialWrapper::getStrainRate()
{
  return theMaterial->getStrainRate();
}

int
DegradingUniaxialWrapper::commitState()
{
  dstate.past = dstate.test;
  return theMaterial->commitState();
}

int
DegradingUniaxialWrapper::revertToLastCommit()
{
  return theMaterial->revertToLastCommit();
}

int
DegradingUniaxialWrapper::revertToStart()
{
  return theMaterial->revertToStart();
}


UniaxialMaterial *
DegradingUniaxialWrapper::getCopy()
{
  return  new DegradingUniaxialWrapper(this->getTag(), *theMaterial, data);
}


int
DegradingUniaxialWrapper::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
DegradingUniaxialWrapper::recvSelf(int cTag, Channel &theChannel,
                                   FEM_ObjectBroker &theBroker)
{
  return -1;
}


void
DegradingUniaxialWrapper::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << WRAPPER_CMD ", tag: " << this->getTag() << "\n";
    s << "  material: " << theMaterial->getTag() << "\n";
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";
    s << "\"material\": " << theMaterial->getTag() << ", ";
    for (int i=0; i<2; ++i) {
      if (i == 0)
        s << "\"pos\": {";
      else
        s << "\"neg\": {";
      s << "\"Cwc\": " << data.idx[i].Cwc << ", ";
      s << "\"Ccd\": " << data.idx[i].Ccd << "}";
      if (i == 0)
        s << ", ";
    }
    s << "}";
  }
}

int
DegradingUniaxialWrapper::setParameter(const char **argv, int argc,
                                       Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}

int
DegradingUniaxialWrapper::updateParameter(int parameterID, Information &info)
{
  return 0;
}

double
DegradingUniaxialWrapper::getStressSensitivity(int gradIndex, bool conditional)
{
  return theMaterial->getStressSensitivity(gradIndex, conditional);
}

double
DegradingUniaxialWrapper::getStrainSensitivity(int gradIndex)
{
  return theMaterial->getStrainSensitivity(gradIndex);
}

double
DegradingUniaxialWrapper::getInitialTangentSensitivity(int gradIndex)
{
  return theMaterial->getInitialTangentSensitivity(gradIndex);
}

double
DegradingUniaxialWrapper::getDampTangentSensitivity(int gradIndex)
{
  return theMaterial->getDampTangentSensitivity(gradIndex);
}

double
DegradingUniaxialWrapper::getRhoSensitivity(int gradIndex)
{
  return theMaterial->getRhoSensitivity(gradIndex);
}

int
DegradingUniaxialWrapper::commitSensitivity(double strainGradient,
                                            int gradIndex, int numGrads)
{
  return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
}


static inline double
sgn(double val) {
  return (double(0) < val) - (val < double(0));
}



int
DegradingUniaxialWrapper::applyDamage(
               const Data&  data,
               DamageState& state,
               const double inputs[3],
               double response[2])
{

  enum {pos,   neg};
  static constexpr double gtol = 1e-8;

  // Workspace variables
  double  Dpsi   [2],
          CwcIF  [2],
          Hshe   [2],
          DdDpsiM[2],
          DPsi   [2];

  // Dereference/cast input data structures
  DamageState::Step& past = state.past;
  DamageState::Step& pres = state.test;

  // extract effective state variables
  pres.strain        = inputs[0];
  double base_stress = inputs[1];

  // extract past history variables
  VectorND<2> d     = past.d;
  VectorND<2> vpEx  = past.vpEx;
  VectorND<2> psiP  = past.psi;
  VectorND<2> psiEx = past.psiEx;
  // plastic deformation increments
  double Dvp = pres.strain - past.strain;
  double vp[2] = {pres.strain, past.strain};
  
  //
  // State update algorithm
  //

  // splitting of effective hinge forces into +ve and -ve components (N, Mi, Mj)
  pres.base_stress[pos]  = base_stress > 0.0 ? base_stress : 0.0;
  pres.base_stress[neg]  = base_stress < 0.0 ? base_stress : 0.0;

  // positive loading
  if ((Dvp > 0.0) && (vp[PRES] > vpEx[pos])) {
    CwcIF[pos] = 1.0;
    if (vp[PAST] < vpEx[pos]) {
      Dvp = vp[PRES] - vpEx[pos] + data.idx[pos].Cwc*(vpEx[pos] - vp[PAST]);
    }
    vpEx[pos] = vp[PRES];

  } else {
    CwcIF[pos] = data.idx[pos].Cwc;
  }

  // negative loading
  if ((Dvp < 0.0) && (vp[PRES] < vpEx[neg])) {
    CwcIF[neg] = 1.0;
    if (vp[PAST] > vpEx[neg])
      Dvp = vp[PRES] - vpEx[neg] + data.idx[neg].Cwc * (vpEx[neg] - vp[PAST]);
    vpEx[neg] = vp[PRES];
  } else {
    CwcIF[neg] = data.idx[neg].Cwc;
  }

  // energy dissipation
  // initialize the energy increment Dpsi
  Dpsi[pos] = ((pres.base_stress[pos] + past.base_stress[pos])*0.5*Dvp);
  Dpsi[neg] = ((pres.base_stress[neg] + past.base_stress[neg])*0.5*Dvp);

  // calculate DPsi
  // effect of coupling between +ve and -ve force (coupling damage coefficient Ccd)
  DPsi[pos] = CwcIF[pos]*Dpsi[pos] + CwcIF[neg]*data.idx[pos].Ccd*Dpsi[neg];
  DPsi[neg] = CwcIF[neg]*Dpsi[neg] + CwcIF[pos]*data.idx[neg].Ccd*Dpsi[pos];

  
  // update energy psi
  VectorND<2> psi{};
  for (int m=0; m < 2; m++) {
    psi[m] = psiP[m] + DPsi[m];
    double g  = psi[m] - psiEx[m];
#if 1
    if (g > gtol) {
        psiEx[m] = psi[m];

        double x = (psi[m] - data.idx[m].psi_d0)
                    / (data.idx[m].psi_d1 - data.idx[m].psi_d0);

        DmgResp dmgresp = data.idx[m].update(x);

        d[m] = std::min(1.0 - data.tol, dmgresp.y);

        DdDpsiM[m] =
            dmgresp.dydx
          / (data.idx[m].psi_d1 - data.idx[m].psi_d0);
    }
#else
    double en_psi = g > data.tol;
    psiEx[m]  = psiEx[m] + en_psi*(psi[m] - psiEx[m]);

    // damage indices
    double psi_tild = (psiEx[m] - data.idx[m].psi_d0) / (data.idx[m].psi_d1 - data.idx[m].psi_d0);
    double DdDx;
    if (std::fabs(psi_tild) < data.tol) {
      DdDx = 0.0;
    }
    else {
      DmgResp dmgresp = data.idx[m].update(psi_tild);
      DdDx = dmgresp.dydx;
      d[m] = std::min(1.-data.tol, dmgresp.y);
    }
    DdDpsiM[m] = en_psi[m]*DdDx / (data.idx[m].psi_d1 - data.idx[m].psi_d0);
  #endif
  }

  // Update stress
  response[0] = (1. -d[pos])*pres.base_stress[pos] + (1.-d[neg])*pres.base_stress[neg];

  // tangent stiffness under damage
  // Heavyside function for derivative of positive/negative forces
  Hshe[pos] = (1. + sgn(base_stress))*0.5; 
  Hshe[neg] = (1. - sgn(base_stress))*0.5;

  //
  // Tangent
  //
  double kt = inputs[2];

#if 0
  // first contribution ka
  double ka = kt*((1.-d[pos])*Hshe[pos] + (1.-d[neg])*Hshe[neg]);

  // second contribution kb: derivative of +ve psi and -ve psi wrt Dv; TODO: optimize
  double DvpDv  = inputs[2];
  double DpsipDv = 0.5*(kt*Dvp*Hshe[pos] + DvpDv*(pres.base_stress[pos] + past.base_stress[pos]));
  double DpsinDv = 0.5*(kt*Dvp*Hshe[neg] + DvpDv*(pres.base_stress[neg] + past.base_stress[neg]));
  double DpsiDv  = DpsipDv + DpsinDv;

  // derivative of +ve Psi wrt Dv = DPsipDv; derivative of -ve Psi wrt Dv = DPsinDv
  // for the axial contribution DPsipDv(1,:) and DPsinDv(1,:)
  double DPsipDv = CwcIF[pos]*DpsipDv + CwcIF[neg]*data.idx[pos].Ccd*DpsinDv;
  double DPsinDv = CwcIF[neg]*DpsinDv + CwcIF[pos]*data.idx[neg].Ccd*DpsipDv;

  // scaling of DPsiDv by DdDpsiM to get DdpDv and DdnDv
  double DdpDv = DPsipDv * DdDpsiM[pos];
  double DdnDv = DPsinDv * DdDpsiM[neg];

  // for second stiffness contribution scale damage index derivatives
  // by +ve, -ve hinges forces
  kb = - DdpDv*pres.base_stress[pos] - DdnDv*pres.base_stress[neg];

  response[1] = ka + kb;
#else 
  {

    double DseAvg_p =
        0.5*(pres.base_stress[pos] + past.base_stress[pos]
          + kt*Hshe[pos]*Dvp);

    double DseAvg_n =
        0.5*(pres.base_stress[neg] + past.base_stress[neg]
          + kt*Hshe[neg]*Dvp);

    double DPsipDv =
        CwcIF[pos]*DseAvg_p
      + data.idx[pos].Ccd*CwcIF[neg]*DseAvg_n;

    double DPsinDv =
        data.idx[neg].Ccd*CwcIF[pos]*DseAvg_p
      + CwcIF[neg]*DseAvg_n;

    double DdpDv = DdDpsiM[pos]*DPsipDv;
    double DdnDv = DdDpsiM[neg]*DPsinDv;

    response[1] =
        kt*((1.0 - d[pos])*Hshe[pos] + (1.0 - d[neg])*Hshe[neg])
      - DdpDv*pres.base_stress[pos]
      - DdnDv*pres.base_stress[neg];
  }
#endif


  // Update current damage history variables
  pres.d     = d;
  pres.vpEx  = vpEx;
  pres.psi   = psi;
  pres.psiEx = psiEx;

  state.past = past;
  state.test = pres;
 
  return 0;
}

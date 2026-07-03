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
// Written in Matlab: Thanh Do
// Created: 07/26
//
#include <cmath>
#include "FariaPlasticDamage3d.h"
#include <Channel.h>
#include <MatrixND.h>
#include <Voight.hpp>
#include <MaterialResponse.h>

#define MND
using namespace OpenSees;


FariaPlasticDamage3d::FariaPlasticDamage3d(int tag,
                                           double E, 
                                           double nu, 
                                           double Ft,
                                           double Fc, 
                                           double beta, 
                                           double Ap,
                                           double An, 
                                           double Bn,
                                           double rho)
: NDMaterial(tag,ND_TAG_PlasticDamageConcrete3d),
  E(E), nu(nu),
  ft(std::fabs(Ft)), Fc(std::fabs(Fc)), 
  beta(beta), Ap(Ap), An(An), Bn(Bn),
  retTangent(C),
  retInitialTangent(Ce),
  retStress(sig),
  retStrain(eps),
  density(rho)
{
  this->revertToStart();
  this->commitState();
}

FariaPlasticDamage3d::FariaPlasticDamage3d()
  : NDMaterial(0, ND_TAG_PlasticDamageConcrete3d)
  , retTangent(C)
  , retInitialTangent(Ce)
  , retStress(sig)
  , retStrain(eps)
  , density(0.0)
{

}

FariaPlasticDamage3d::~FariaPlasticDamage3d()
{

}


// Stress invariants: octahedral normal and shear
static inline void
StrsInvar(const VectorND<6> &sig, double &sigoct, double &tauoct)
{
  // normal stress
  sigoct = (sig(0) + sig(1) + sig(2))/3.0;
  // shear stress
  double J2 = (std::pow((sig(0) - sig(1)),2) 
            +  std::pow((sig(1) - sig(2)),2)
            +  std::pow((sig(2) - sig(0)),2))/6.0
            +  std::pow(sig(3),2)
            +  std::pow(sig(4),2)
            +  std::pow(sig(5),2);
  tauoct = std::sqrt((2./3.)*J2);
}


//
// Stress decomposition function: algebraic approach
//
#if 1
#include <concrete/StrsDec.cpp>
#else
template <typename T> static inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
static inline void
StrsDecA(const Vector &sig, VectorND<6> &sigpos, //VectorND<6> &signeg, 
                            MatrixND<6,6> *Qpos) //,   MatrixND<6,6> *Qneg)
{
#if 1
  for (int i=0; i<6; i++) {
    // positive and negative stress tensors
    sigpos(i) = (sig(i) + std::fabs(sig(i))) / 2.0;
  }

  // projection tensors
  if (Qpos == nullptr)
    return;

  for (int i=0; i<6; i++) {
    (*Qpos)(i,i) = (1.0 + sgn(sig(i))) / 2.0;
  }
  return;
#else
  for (int i=0; i<6; i++) {
    if (sig(i) > 1e-8) {
      sigpos(i) = sig(i);
      signeg(i) = 0.;
    } else if (sig(i)  < -1.e-8) {
      sigpos(i) = 0.;
      signeg(i) = sig(i);
    } else {
      sigpos(i) = sig(i)/2.;
      signeg(i) = sig(i)/2.;
    }
  }
  if (Qpos == nullptr || Qneg == nullptr)
    return;
  for (int i=0; i<6; i++) {
    if (sig(i) > 1e-8) {
      (*Qpos)(i,i) = 1;
      (*Qneg)(i,i) = 0;
    } else if (sig(i)  < -1.e-8) {
      (*Qpos)(i,i) = 0;
      (*Qneg)(i,i) = 1;
    } else {
      (*Qpos)(i,i) = 0.5;
      (*Qneg)(i,i) = 0.5;
    }
  }
#endif
}
#endif


// Compute dn = 1 - (rn0/rn)*(1-An) - An*exp(Bn*(1 - rn/rn0))
static double
negative_damage(double rn0, double rn, double An, double Bn)
{
  // 1) Compute alpha = rn0/rn
  double alpha = rn0 / rn;

  // 2) Compute x = Bn * (1 - rn/rn0) without cancellation
  double inv_alpha = 1.0 / alpha;
  double x = Bn * (1.0 - inv_alpha);

  double exm1 = std::expm1(x);

  // 4) Combine terms. Use fma to reduce rounding:
  //    term1 = (1 - alpha)*(1 - An)
  double term1 = std::fma(-(1.0 - An), alpha, (1.0 - An));  
  //    which is = (1 - An) - alpha*(1 - An)
  //
  double term2 = An * exm1;

  // 5) dn = term1 - term2
  double dn = std::fma(-An, exm1, term1);
  return dn;
}

// Compute dp = 1 - (rp0/rp)*exp(Ap*(1 - rp/rp0))
// with improved precision when rp is near rp0
static double
positive_damage(double rp0, double rp, double Ap)
{
  // 1) alpha = rp0 / rp
  double alpha = rp0 / rp;

  // 2) x = Ap * (1 - rp/rp0) = Ap * (1 - 1/alpha)
  double inv_alpha = 1.0 / alpha;
  double x = Ap * (1.0 - inv_alpha);

  // 3) exm1 = exp(x) - 1
  double exm1 = std::expm1(x);

  // 4)
  double term1 = 1.0 - alpha;

  // 5) dp = term1 - alpha * exm1
  //    use fma to compute -(alpha * exm1) + term1 in one rounding
  double dp = std::fma(-alpha, exm1, term1);

  return dp;
}


static double
negative_equivalent_stress(const VectorND<6> &signeg, double k)
{
  double sigoct, tauoct;
  StrsInvar(signeg, sigoct, tauoct);                  //  find octahedral stresses
  double eta = k*sigoct + tauoct;
  if (eta <= 0.0)
    return 0.0;
  return std::sqrt(std::sqrt(3.0)*eta);
}


static double
negative_surface(const VectorND<6> &signeg, double rn, double k)
{
  return negative_equivalent_stress(signeg, k) - rn;    //  negative equivalent stress
}

int
FariaPlasticDamage3d::setTrialStrain(const Vector &strain)
{
  bool plasticCorrection = false;
  bool positiveDamageActive = false;
  bool negativeDamageActive = false;


  double F2c = 1.16*Fc; // f2c: biaxial compressive strength
  double k = std::sqrt(2.0)*(F2c - Fc)/(2.0*F2c - Fc);
  // initial damage threshold
  double rp0 = ft/std::sqrt(E);
  double rn0 = std::sqrt((-k + std::sqrt(2.0))*Fc/std::sqrt(3.0));

  constexpr static double toln = 1e-5,
                          tolp = 1e-5;
  constexpr static double dn_max = 1.0 - 1e-5,
                          dp_max = 1.0 - 1e-5;

  // retrieve history variables
  eps_p = eps_pCommit;
  VectorND<6> sigeP = sigeCommit;
  rp = rpCommit;
  rn = rnCommit;
  dp = dpCommit;
  dn = dnCommit;

  // elastic compliance
  MatrixND<6,6> Se{};
  Ce.invert(Se);

  // current strain
  for (int i=0; i<6; i++)
    eps(i) = strain(i);

  // incremental strain
  VectorND<6> Depse_tr = eps - eps_p; // ?[cmp]
  VectorND<6> Deps = eps - epsCommit;


  //
  // PLASTIC part
  //

  // elastic trial stress
  VectorND<6> sige_tr  =  sigeCommit;
  sige_tr.addMatrixVector(Ce, Deps, 1.0);

  // decomposition of trial stress tensor
  VectorND<6> sigpos{}, signeg = sige_tr;
  StrsDecA(sige_tr, sigpos, nullptr);
  signeg -= sigpos;

  MatrixND<6,6> Cbar{};
  if (negative_surface(signeg, rn, k) <= (toln*rn0)) {
    // elastic state, accept trial response
    sige = sige_tr;
    Cbar = Ce;
  }
  else {
    // Correction

    //  norm of trial effective stress
    double nrm = std::sqrt(
                 sige_tr(0)*sige_tr(0)
               + sige_tr(1)*sige_tr(1)
               + sige_tr(2)*sige_tr(2)
               + 2.0*sige_tr(3)*sige_tr(3)
               + 2.0*sige_tr(4)*sige_tr(4)
               + 2.0*sige_tr(5)*sige_tr(5));

    // normalized trial effective stress
    VectorND<6> L_tr = sige_tr/nrm;

    double L_trDotDeps = L_tr.dot(Deps);

    // plastic strain increment
    VectorND<6> Deps_p = Depse_tr*(beta*E/nrm*L_trDotDeps);

    double lam  = 1.0 - beta*E/nrm * L_trDotDeps;

    sige.addVector(0.0, sige_tr, lam);                 //  corrected effective stress

    // check damage
    signeg = sige;
    StrsDecA(sige, sigpos, nullptr);                   //  decompose the effective stress  
    signeg -= sigpos;

    if ((negative_surface(signeg, rn, k) <= toln*rn0) || (L_trDotDeps <= 0.0)) {
      //  no damage or sige and eps in different direction
      sige = sige_tr;
      Cbar = Ce;
    }

    else {
      plasticCorrection = true;

      //  update plastic strain
      eps_p += Deps_p;

      VectorND<6> L_tr_temp {};
      // tangent in effective space, Cbar
      for (int i=0; i<3; i++) 
        L_tr_temp(i) = L_tr(i);
      for (int i=3; i<6; i++)
        L_tr_temp(i) = 2.0*L_tr(i);

      double Dlam_Dnrm = 2.0*beta*E/std::pow(nrm,3)*sige_tr.dot(Deps);

      VectorND<6> Dlam_Deps = L_tr;
      Dlam_Deps *= -beta*E/nrm;       // Dlam_Deps = -beta*E/nrm*L_tr;
      // Dlam_Deps = Dlam_Dnrm * Ce * Dnrm_Dsig + Ce*Dlam_Dsig + Dlam_Deps; 
      Dlam_Deps.addMatrixVector(1.0, Ce, L_tr_temp,         Dlam_Dnrm);
      Dlam_Deps.addMatrixVector(1.0, Ce,      Deps, -beta*E/(nrm*nrm));

      Cbar.addMatrix(Ce, lam);
      Cbar.addTensorProduct(sige_tr, Dlam_Deps, 1.0);
    }
  }


  //
  // DAMAGE part
  //

  // decompose into positive and negative effective stress tensor
  MatrixND<6,6> Qpos{}, Qneg = Voight::IImix;
  signeg = sige;
  StrsDecA(sige, sigpos, &Qpos);    // decompose the effective stress
  signeg -= sigpos;                 // signeg = sige - sigpos
  Qneg   -= Qpos;

  // positive damage
  double taup,
         Ddp_Drp = 0.;
  {
  // calculate equivalent stresses
    VectorND<6> tmp = Se*sigpos; // {}
    // Ce.solve(sigpos, tmp);
    taup = std::sqrt(tmp.dot(sigpos));      // positive equivalent stress


    if ((taup - rp) <= (tolp*rp0)) {               // no positive damage
      Ddp_Drp = 0;
    }
    else {                                         // positive damage evolves
      positiveDamageActive = true;
      rp = taup;                                   // update rp = max(taup, rp)
      // dp = 1. - rp0/rp * std::exp(Ap*(1. - rp/rp0));
      dp = positive_damage(rp0, rp, Ap);          // update dp
      Ddp_Drp =  (Ap*rp + rp0)/(rp*rp) * std::exp(Ap*(1.0 - rp/rp0));
      if (dp < 0) {
        dp = 0;
        Ddp_Drp = 0;
      }       
      // dp = dp*dp_max;                             // cap the damage variable 
      // Ddp_Drp = Ddp_Drp*dp_max;
      if (dp > dp_max) {
        dp = dp_max; 
        Ddp_Drp = 0;
      }
    }
  }

  // negative damage
  double Ddn_Drn = 0;
  double gn = negative_surface(signeg, rn, k); // negative equivalent stress
  {
    if (gn <= toln*rn0) {                    // no negative damage
      Ddn_Drn = 0;
    }
    else {                                  // negative damage evolves
      negativeDamageActive = true;
      // rn = taun;                         // update rn
      rn += gn;
      // dn = 1. - rn0/rn*(1-An) - An*std::exp(Bn*(1. - rn/rn0));
      dn = negative_damage(rn0, rn, An, Bn); // update dn
      Ddn_Drn = rn0/(rn*rn)*(1.0-An) + An*Bn/rn0*std::exp(Bn*(1.0 - rn/rn0));
      if (dn < 0) {
        dn = 0;
        Ddn_Drn = 0;
      }
      // cap the damage variable
      if (dn > dn_max) {
        dn = dn_max;
        Ddn_Drn = 0;
      }
    }
  }

  // stress update
  sig.addVector(0.0, sigpos, 1.0-dp);
  sig.addVector(1.0, signeg, 1.0-dn);

  //
  // TANGENT
  //
  {
    C.zero();

    // Negative part
    {
      MatrixND<6,6> Dsigneg_Deps = Qneg*Cbar;

      VectorND<6> Dtaun_Dsigneg{};
      double sigoct, tauoct;
      StrsInvar(signeg, sigoct, tauoct);             // find octahedral stresses
      double eta = k*sigoct + tauoct;
      double taun = eta <= 0.0 ? 0.0 : std::sqrt(std::sqrt(3.0)*eta);   // negative equivalent stress
      if (taun <= toln) {
        Dtaun_Dsigneg.zero();
      }
      else {
        // norm of deviatoric stress
        VectorND<6> s = Voight::IIdevMix*signeg;            // deviatoric stress

        double nrms = std::sqrt(s[0]*s[0] + s[1]*s[1] + s[2]*s[2]
                         + 2.0*(s[3]*s[3] + s[4]*s[4] + s[5]*s[5]));
        VectorND<6> n = s;
        double sqrt_eta = std::sqrt(eta);
        double Dtaun_Dsigoct = std::pow(3.0,0.25) * k/2.0/sqrt_eta;
        Dtaun_Dsigneg.addVector(0.0, Voight::ivol, Dtaun_Dsigoct/3.0); // = Dtaun_Dsigoct * Dsigoct_Dsigneg
        if (std::abs(nrms) <= 1e-8) //toln) 
          n.zero();
        else {
          n/=nrms;
          double Dtaun_Dtauoct = std::pow(3.0,0.25) / 2.0/sqrt_eta;
          Dtaun_Dsigneg.addVector(1.0,    n, Dtaun_Dtauoct/std::sqrt(3.0)); // += Dtaun_Dtauoct * Dtauoct_Dsigneg
        }

        for (int i=3; i<6; i++)
          Dtaun_Dsigneg(i) *= 2.0;
      }

      VectorND<6> Ddn_Deps = Dsigneg_Deps ^ Dtaun_Dsigneg;
      Ddn_Deps *= Ddn_Drn;
      C.addMatrix(Dsigneg_Deps , (1.0 - dn));
      C.addTensorProduct(signeg, Ddn_Deps, -1.0);
    }

    // Positive part
    {
      MatrixND<6,6> Dsigpos_Deps = Qpos*Cbar;
      VectorND<6> Dtaup_Dsigpos{};
      if (taup <= tolp) {
        Dtaup_Dsigpos.zero();
      }
      else  {
        Dtaup_Dsigpos = Se*sigpos;
        // Ce.solve(sigpos, Dtaup_Dsigpos);
        Dtaup_Dsigpos/=taup;
      }

      // Ddp_Deps = Ddp_Drp * Dsigpos_Deps' * Dtaup_Dsigpos;
      VectorND<6> Ddp_Deps{};
      Ddp_Deps.addMatrixTransposeVector(0.0, Dsigpos_Deps, Dtaup_Dsigpos, Ddp_Drp);

      C.addMatrix(Dsigpos_Deps , (1.0 - dp));
      C.addTensorProduct(sigpos, Ddp_Deps, -1.0);
    }
  }

  state[0] = (plasticCorrection || positiveDamageActive || negativeDamageActive)
      ? State::Plastic
      : State::Elastic;
  state[1] = positiveDamageActive 
           ? State::Plastic : State::Elastic; // positive damage state
  state[2] = (negativeDamageActive || plasticCorrection)
           ? State::Plastic : State::Elastic; // negative damage state
  return 0;
}


int
FariaPlasticDamage3d::setTrialStrain(Vector const&v1, Vector const&v2)
{
  return this->setTrialStrain(v1);
}


int
FariaPlasticDamage3d::setTrialStrainIncr(const Vector &strain)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

int
FariaPlasticDamage3d::setTrialStrainIncr(const Vector &strain, const Vector &rate)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

const Matrix&
FariaPlasticDamage3d::getTangent()
{
  return retTangent;
}

const Matrix&
FariaPlasticDamage3d::getInitialTangent()
{
  return retInitialTangent;
}

const Vector&
FariaPlasticDamage3d::getStress()
{
  return retStress;
}

const Vector&
FariaPlasticDamage3d::getStrain()
{
  return retStrain;
}

int
FariaPlasticDamage3d::commitState()
{
  rpCommit    = rp;
  rnCommit    = rn;
  dpCommit    = dp;
  dnCommit    = dn;

  epsCommit   = eps;
  sigCommit   = sig;
  sigeCommit  = sige;
  eps_pCommit = eps_p;

  Ccommit = C;
  return 0;
}

int
FariaPlasticDamage3d::revertToLastCommit()
{
  rp = rpCommit;
  rn = rnCommit;
  dp = dpCommit;
  dn = dnCommit;

  eps = epsCommit;
  sig = sigCommit;
  sige = sigeCommit;
  eps_p = eps_pCommit;
  C  = Ccommit;
  return 0;
}

int
FariaPlasticDamage3d::revertToStart()
{
  eps.zero();
  sig.zero();
  sige.zero();
  eps_p.zero();

  double G  = E/2./(1. +    nu);     // Shear modulus
  double K  = E/3./(1. - 2.*nu);     // Bulk  modulus
  Ce.zero();
  Ce.addMatrix(Voight::IIvol, K);
  Ce.addMatrix(Voight::IIdevCon, 2.*G);
  // Ce.addMatrix(IIdevMix, 2*G);

  C = Ce;
  Ccommit = Ce;

  double F2c = 1.16*Fc;
  double k = std::sqrt(2.0)*(F2c - Fc)/(2.0*F2c - Fc);

  // initial damage threshold
  double rp0 = ft/std::sqrt(E);
  double rn0 = std::sqrt((-k+std::sqrt(2.0))*Fc/std::sqrt(3.0));

  rp = rp0;
  rn = rn0;
  dp = 0.;
  dn = 0.;
  rpCommit = rp;
  rnCommit = rn;
  dpCommit = dp;
  dnCommit = dn;
  return 0;
}

NDMaterial*
FariaPlasticDamage3d::getCopy(const char *type)
{
  if (strcmp(type,"ThreeDimensional") == 0) {
    return this->getCopy();
  }
  else {
    return NDMaterial::getCopy(type);
  }
}

NDMaterial*
FariaPlasticDamage3d::getCopy()
{
  auto *copy = new FariaPlasticDamage3d(this->getTag(),
                                        E, nu,
                                        ft,
                                        Fc,
                                        beta,
                                        Ap,
                                        An,
                                        Bn,
                                        density);

  copy->rp = rp;
  copy->rn = rn;
  copy->dp = dp;
  copy->dn = dn;
  copy->eps = eps;
  copy->sig = sig;
  copy->sige = sige;
  copy->eps_p = eps_p;

  copy->rpCommit = rpCommit;
  copy->rnCommit = rnCommit;
  copy->dpCommit = dpCommit;
  copy->dnCommit = dnCommit;
  copy->epsCommit = epsCommit;
  copy->sigCommit = sigCommit;
  copy->sigeCommit = sigeCommit;
  copy->eps_pCommit = eps_pCommit;
  copy->Ce = Ce;
  copy->C = C;
  copy->Ccommit = Ccommit;

  return copy;
}

const char*
FariaPlasticDamage3d::getType() const
{
  return "ThreeDimensional";
}

int
FariaPlasticDamage3d::getOrder() const
{
  return 6;
}


int 
FariaPlasticDamage3d::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(1+8+4 + 5*6);

  data(0) = this->getTag();

  data(1) = E;
  data(2) = nu;
  data(3) = ft;
  data(4) = Fc;
  data(5) = beta;
  data(6) = Ap;
  data(7) = An;
  data(8) = Bn;

  data(9)  = rpCommit;
  data(10) = rnCommit;
  data(11) = dpCommit;
  data(12) = dnCommit;  

  for (int i = 0; i < 6; i++) {
    data(13+i) = epsCommit(i);
    data(13+6+i) = sigCommit(i);
    data(13+12+i) = sigeCommit(i);
    data(13+18+i) = eps_pCommit(i);
    // data(13+24+i) = sigePCommit(i);
  }

  int res = 0;
  int dbTag = this->getDbTag();

  res = theChannel.sendVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FariaPlasticDamage3d::sendSelf -- could not send Vector\n";
    return res;
  }
  Matrix Ccm(Ccommit);
  res = theChannel.sendMatrix(dbTag, commitTag, Ccm);
  if (res < 0) {
    opserr << "FariaPlasticDamage3d::sendSelf -- could not send Ccommit matrix\n";
    return res;
  }  
  
  return res;
}

int 
FariaPlasticDamage3d::recvSelf(int commitTag, Channel &theChannel, 
				  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dbTag = this->getDbTag();

  static Vector data(43);  
  res = theChannel.recvVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FariaPlasticDamage3d::recvSelf -- could not receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));

  E = data(1);
  nu = data(2);
  ft = data(3);
  Fc = data(4);
  beta = data(5);
  Ap = data(6);
  An = data(7);
  Bn = data(8);

  rpCommit = data(9);
  rnCommit = data(10);
  dpCommit = data(11);
  dnCommit = data(12);

  for (int i = 0; i < 6; i++) {
    epsCommit(i) = data(13+i);
    sigCommit(i) = data(13+6+i);
    sigeCommit(i) = data(13+12+i);
    eps_pCommit(i) = data(13+18+i);
  }
  
  Matrix Ccm(Ccommit);
  res = theChannel.recvMatrix(dbTag, commitTag, Ccm);
  if (res < 0) {
    opserr << "FariaPlasticDamage3d::recvSelf -- could not receive Ccommit matrix\n";
    return res;
  }  

  double G  = E/2./(1. +    nu);     // Shear modulus
  double K  = E/3./(1. - 2.*nu);     // Bulk  modulus
  Ce.zero();
  Ce.addMatrix(Voight::IIvol, K);
  Ce.addMatrix(Voight::IIdevCon, 2.*G);

  rp = rpCommit;
  rn = rnCommit;
  dp = dpCommit;
  dn = dnCommit;
  eps = epsCommit;
  sig = sigCommit;
  sige = sigeCommit;
  eps_p = eps_pCommit;
  C = Ccommit;
  
  return res;
}

void 
FariaPlasticDamage3d::Print(OPS_Stream &s, int flag) 
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";
    s << "\"E\": " << E << ", ";
    s << "\"nu\": " << nu << ", ";
    s << "\"Ft\": " << ft << ", ";
    s << "\"Fc\": " << Fc << ", ";
    s << "\"beta\": " << beta << ", ";
    s << "\"Ap\": " << Ap << ", ";
    s << "\"An\": " << An << ", ";
    s << "\"Bn\": " << Bn;
    s << "}";
  }

  else {
    opserr << this->getType() << ": " << this->getTag() << "\n";
    opserr << "stress: " << Vector(retStress) << "\n";
    opserr << "strain: " << Vector(retStrain) << "\n";
    opserr << "tangent: " << Matrix(retTangent) << "\n";
  }
}       




Response*
FariaPlasticDamage3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = nullptr;
  const char *matType = this->getType();

  output.tag("NdMaterialOutput");
  output.attr("matType",this->getClassType());
  output.attr("matTag",this->getTag());

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) {
    const Vector &res = this->getStress();
    int size = res.Size();
    
    if ((strcmp(matType,"PlaneStress") == 0 && size == 3) ||
        (strcmp(matType,"PlaneStrain") == 0 && size == 3)) {
      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");
    } else if (strcmp(matType,"ThreeDimensional") == 0 && size == 6) {
      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma33");
      output.tag("ResponseType","sigma12");
      output.tag("ResponseType","sigma23");
      output.tag("ResponseType","sigma13");
    } else {
      for (int i=0; i<size; i++) 
      output.tag("ResponseType","UnknownStress");
    }
    theResponse =  new MaterialResponse(this, 1, this->getStress());

  } 
  else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
    const Vector &res = this->getStrain();
    int size = res.Size();
    if ((strcmp(matType,"PlaneStress") == 0 && size == 3) ||
        (strcmp(matType,"PlaneStrain") == 0 && size == 3)) {
      output.tag("ResponseType","eta11");
      output.tag("ResponseType","eta22");
      output.tag("ResponseType","eta12");
    } else if (strcmp(matType,"ThreeDimensional") == 0 && size == 6) {
      output.tag("ResponseType","eps11");
      output.tag("ResponseType","eps22");
      output.tag("ResponseType","eps33");
      output.tag("ResponseType","eps12");
      output.tag("ResponseType","eps23");
      output.tag("ResponseType","eps13");
    } else {
      for (int i=0; i<size; i++) 
        output.tag("ResponseType","UnknownStrain");
    }      
    theResponse =  new MaterialResponse(this, 2, this->getStrain());
  }
  else if (strcmp(argv[0], "isElastic") == 0) {
    double data[1];
    theResponse = new MaterialResponse(this, 20, Vector(&data[0], 1));
  }
  else if (strcmp(argv[0], "damage") == 0) {
    double data[2];
    theResponse = new MaterialResponse(this, 21, Vector(&data[0], 2));
  }
  output.endTag(); // NdMaterialOutput

  return theResponse;
}


int
FariaPlasticDamage3d::getResponse(int responseID, Information &info)
{
  switch (responseID) {
    case 1:
      return info.setVector(this->getStress());

    case 2:
      return info.setVector(this->getStrain());

    case 4:
      return info.setMatrix(this->getTangent());

    case 20: {
      double data[3];
      data[0] = (state[0] == State::Elastic) ? 1.0 : 0.0;
      data[1] = (state[1] == State::Elastic) ? 1.0 : 0.0;
      data[2] = (state[2] == State::Elastic) ? 1.0 : 0.0;
      return info.setVector(Vector(&data[0], 3));
    }

    case 21: {
      double data[2];
      data[0] = dp;
      data[1] = dn;
      return info.setVector(Vector(&data[0], 2));
    }

    default: {
      return -1;
    }
  }
}

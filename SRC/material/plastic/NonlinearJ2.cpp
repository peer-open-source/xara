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
// General isotropic J2 plasticity with nonlinear isotropic/kinematic hardening
//
// Details and references are provided in the header file NonlinearJ2.h.
//
//
// material J2 tag E nu 
//      -F {Fy Fi Fs}  -b {a b} 
//      -C {C1 C2 ...} -g {g1 g2 ...}
//      -YFtol tol -maxIter n -density rho
//===----------------------------------------------------------------------===//
//
// Written: Claudio M. Perez
//
#include "NonlinearJ2.h"
#include <VectorND.h>
#include <Voigt.hpp>
#include <MatrixND.h>
#include <Logging.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <cmath>
#include <limits>
#include <cstring>
// #include "hardening/FlowStress.h"

using namespace OpenSees;

static const double EPS_DBL = std::numeric_limits<double>::epsilon();

NonlinearJ2::NonlinearJ2(int tag,
                     double E_, double nu_,
                     double Fy_,
                     double density,
                     double Hiso,
                     double a, double DInf,
                     double b, double QInf,
                     const std::vector<double> &Ck,
                     const std::vector<double> &gammak,
                     double YFtol,
                     int MaxIter)
: NDMaterial(tag, 0),
  E(E_), nu(nu_), Fy(Fy_),
  Hiso_(Hiso), Q_{DInf, QInf}, b_{a, b},
  newton_tolerance(YFtol), MaxIter_(MaxIter),
  density_(density),
  retInitialTangent(Ce), retTangent(C),
  retStress(pres.sig), retStrain(eps)
{
  G = E/2./(1. +    nu);
  int nc = Ck.size();
  for (int i=0;i<nc;i++) {
    Ck_.push_back((2.0/3.0)*Ck[i]);
    gammak_.push_back(gammak[i]);
  }
  this->revertToStart();
  this->commitState();
}

NonlinearJ2::~NonlinearJ2() = default;

//
// Core state determination
//
int
NonlinearJ2::updateState(const VectorND<6> &eps)
{

  // const double G  = E/2./(1. +    nu);  // Shear modulus
  const double K  = E/3./(1. - 2.*nu);  // Bulk  modulus

  // Deviatoric strain
  const VectorND<6> eps_dev = Voigt::Dev(eps);
  const double      eps_vol = Voigt::Trace(eps);

  // Trial deviatoric stress s_tr = 2G * (eps_dev - eps_p_d)
  const VectorND<6> eps_p_d = Voigt::Dev(past.eps_p);
  VectorND<9> s_tr{};
  s_tr.addVector(0.0, Voigt::ExpandVector(eps_dev),  2.0*G);
  s_tr.addVector(1.0, Voigt::ExpandVector(eps_p_d), -2.0*G);

  // Sum of back-stresses at last commit
  VectorND<9> Xn{};
  for (const auto &x : past.sig_b)
    Xn.addVector(1.0, x, 1.0);

  // Yield function at trial
  const double sig_nrm = metric(s_tr - Xn);

  const double f_tr = sig_nrm - Isotropic::Y(*this, past, 0.0)*phi_n*mises_scale;

  double g_tol, f_tol;
  this->scale_tolerance(sig_nrm, g_tol, f_tol);


  // Elastic admissible?
  if (f_tr <= f_tol) {
    // Stress = P*s_tr + K*eps_v*ivol
    pres.sig = Voigt::ReduceVector(s_tr);
    Voigt::AddVol(pres.sig, K*eps_vol);
    // pres.sig = P * s_tr;
    // pres.sig.addVector(1.0, Voigt::ivol,  K*eps_vol);

    pres.eps             = eps;
    pres.eps_p           = past.eps_p;
    pres.mises_strain    = past.mises_strain;
    pres.sig_b           = past.sig_b;

    // Elastic tangent
    C.zero();
    C.addMatrix(Voigt::IoI,          K);
    C.addMatrix(Voigt::IIdevCon, 2.0*G);

    // Currently not computing energy
    if (false) {
      // Energy computations (use avg stress with Deps)
      const VectorND<6> sigav = (pres.sig + past.sig) * 0.5;
      const VectorND<6> Deps  = pres.eps - past.eps;
      const VectorND<6> Deps_p = pres.eps_p - past.eps_p;
      pres.Estr = past.Estr + sigav.dot(Deps);
      pres.Ehst = past.Ehst + sigav.dot(Deps_p);
    }
    return 0;
  }

  //
  // Plastic correction, Newton on lambda
  //
  double lamda   = 0.0;
  double g, Dg, phi_a=0.0;
  int iter = 0;
  VectorND<9> m;
  this->newton_update(s_tr,lamda,   m,g,Dg,phi_a);

  // Solve consistency condition
  while ((std::fabs(g) > g_tol) && (iter < MaxIter_)) {

    Dg = std::abs(Dg)>EPS_DBL ? Dg : ((Dg>=0)?EPS_DBL:-EPS_DBL);
    // Newton step
    lamda -= g / Dg;
    if (lamda < 0.0)
      lamda = 0.0;

    // Recompute Z, n and residual
    this->newton_update(s_tr,lamda,  m,g,Dg,phi_a);

    iter++;
  }

  if (iter == MaxIter_ && std::abs(g) > g_tol) {
    // Currently cannot print since opserr is not thread safe.
    return -1;
  }

  //
  // Plastic updates at converged lamda
  //

  // eps_p += I^{-1} * (P' * lamda*n)
  // const VectorND<6> n6 = P * m;
  const VectorND<6> n6 = Voigt::ReduceVector(m);
  // apply I^{-1} (double the shear components)
  VectorND<6> deps_p_Iinv = n6;
  for (int i=3;i<6;i++)
    deps_p_Iinv[i] *= 2.0;

  pres.eps_p = past.eps_p;
  pres.eps_p.addVector(1.0, deps_p_Iinv, lamda);


  // Update isotropic variables
  const double d_mises_strain = lamda*flow_rate*(phi_n*mises_scale);
  Isotropic::update(*this, past, d_mises_strain, pres);

  // Update back-stress components
  Kinematic::update(*this, past, d_mises_strain, m, pres);

  // Assemble full stress = P*s_dev + K*eps_v*ivol

  // Correct deviatoric stress: s = s_tr - 2G*lamda*n
  pres.sig.addMatrixVector(0.0, P, s_tr, 1.0);
  pres.sig.addMatrixVector(1.0, P, m, -2.0*G*lamda);
  pres.sig.addVector(1.0, Voigt::ivol,  K*eps_vol);

  pres.eps = eps;

  //
  // Algorithmic tangent
  //

  // Effective shear modulus factor
  const double theta_phi = lamda*phi_m / phi_a;

  C.zero();
  // 1. bulk term
  C.addMatrix(Voigt::IoI,      K);
  // 2. scale contravariant deviatoric identity
  C.addMatrix(Voigt::IIdevCon,  2.0*G*(1.0 - 2.0*G*theta_phi));

  const double theta_lam = -2.0*G/(flow_rate*Dg);
  C.addTensorProduct(n6, n6,  -2.0*G*(theta_lam - 2.0*G*theta_phi*phi_nnmm));

  // term with back-stress coupling
  Kinematic::addTangent(*this, past, d_mises_strain, m, -2.0*G*theta_phi*theta_lam, C);


  // Energy
  if (false) {
    const VectorND<6> sigav = (pres.sig + past.sig) * 0.5;
    pres.Estr = past.Estr + sigav.dot(pres.eps - past.eps);
    pres.Ehst = past.Ehst + sigav.dot(pres.eps_p - past.eps_p);
  }

  return 0;
}


void 
NonlinearJ2::newton_update(const VectorND<9>& s_tr, 
                          double lambda,
                          // outputs:
                          VectorND<9>& m,
                          double& g, 
                          double& Dg,
                          double& phi_z) const noexcept 
{
  double d_mises_strain = (flow_rate*lambda)*phi_n*mises_scale;

  // Partial stress, Z
  m = s_tr;
  m -= Kinematic::A(*this, past, d_mises_strain);

  phi_z  = metric(m); // =  phi(Z)

  // Newton residual function
  g = phi_z - 2.0*G*d_mises_strain*(phi_n/mises_scale)
    - Kinematic::H(*this, past, d_mises_strain)*(phi_n/mises_scale)
    - Isotropic::Y(*this, past, d_mises_strain)*phi_n*mises_scale;

  if (phi_z < 1e-16) {
    m.zero();
    Dg = 0.0;
    return;
  }

  m    *= phi_m/phi_z; // m = Z * phi(m)/phi(a)

  // Derivative of newton residual
  // Note: d(norm(Z))/d(lamda) == n . X' == dX(n)
  Dg = Kinematic::dX(*this, past, d_mises_strain, m)*mises_rate/flow_rate
      - Kinematic::dH(*this, past, d_mises_strain)*mises_rate*(phi_n/mises_scale)
      - 2.0*G*mises_rate*(phi_n/mises_scale)
      - Isotropic::H(*this, past, d_mises_strain)*mises_scale*phi_n*mises_rate;
}

void
NonlinearJ2::scale_tolerance(double sig_nrm, double& gtol, double& ftol) const noexcept
{
  const double yieldRadius = Isotropic::Y(*this, past, 0.0) * phi_n * mises_scale;

  const double yieldScale =
      std::max(
          1.0,
          std::max(sig_nrm, std::abs(yieldRadius))
      );

  double round_tol = 256.0
          * std::numeric_limits<double>::epsilon()
          * yieldScale;

  gtol = std::max(
          newton_tolerance,
          round_tol
      );
  ftol = std::max(1e-10*yieldScale, round_tol);
}

//
// NDMaterial API
//

int
NonlinearJ2::setTrialStrain(const Vector &v)
{
  for (int i=0;i<6;i++)
    eps[i] = v(i);
  return updateState(eps);
}

int
NonlinearJ2::setTrialStrainIncr(const Vector &dv)
{
  for (int i=0;i<6;i++)
    eps[i] += dv(i);
  return updateState(eps);
}


int 
NonlinearJ2::commitState()
{
  past = pres;
  return 0;
}


int 
NonlinearJ2::revertToLastCommit()
{
  pres = past;
  C = Ce;
  return 0;
}


int
NonlinearJ2::revertToStart()
{
  // Zero present/past
  eps.zero();

  pres.sig.zero();
  pres.eps.zero();
  pres.eps_p.zero();
  pres.sig_b.assign(Ck_.size(), VectorND<9>{});
  pres.mises_strain = 0.0;
  pres.Estr  = 0.0;
  pres.Ehst  = 0.0;

  past = pres;

  // Elastic stiffness
  const double G  = E/2./(1. +    nu);     // Shear modulus
  const double K  = E/3./(1. - 2.*nu);     // Bulk  modulus
  Ce.zero();
  Ce.addMatrix(Voigt::IoI,          K);
  Ce.addMatrix(Voigt::IIdevCon, 2.0*G);
  C = Ce;

  return 0;
}


const Matrix &
NonlinearJ2::getTangent() 
{
  // return wrapper around C
  return retTangent;
}

const Matrix &
NonlinearJ2::getInitialTangent() 
{
  return retInitialTangent; 
}

const Vector &
NonlinearJ2::getStress() 
{ 
  return retStress;
}

const Vector &
NonlinearJ2::getStrain() 
{ 
  return retStrain;
}


NDMaterial *
NonlinearJ2::getCopy()
{
  auto *m = new NonlinearJ2(this->getTag(),
                          E, nu, 
                          Fy, density_, 
                          Hiso_, b_[0], Q_[0], b_[1], Q_[1],
                          Ck_, gammak_, 
                          newton_tolerance, MaxIter_);
  // constructor assumed Ck without 2/3 and rescales. 
  // We have to undo that here to avoid double rescaling
  for (size_t i=0;i<Ck_.size();i++)
    m->Ck_[i] = Ck_[i];
  m->Ce   = this->Ce;
  m->C    = this->C;
  m->eps  = this->eps;
  m->past = this->past;
  m->pres = this->pres;
  return m;
}


NDMaterial *
NonlinearJ2::getCopy(const char *type)
{
  if (type && std::strcmp(type, "ThreeDimensional") == 0)
    return this->getCopy();

  return NDMaterial::getCopy(type);
}



void
NonlinearJ2::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";
    s << "\"E\": " << E << ", ";
    s << "\"nu\": " << nu << ", ";
    s << "\"Fy\": " << Fy << ", ";
    s << "\"Q\": [" << Q_[0] << ", " << Q_[1] << "], ";
    s << "\"b\": [" << b_[0] << ", " << b_[1] << "], ";
    s << "\"C\": [";
    for (size_t i=0;i<Ck_.size();i++) {
      s << Ck_[i]*(3.0/2.0);
      if (i < Ck_.size()-1) s << ", ";
    }
    s << "], ";
    s << "\"gamma\": [";
    for (size_t i=0;i<gammak_.size();i++) {
      s << gammak_[i];
      if (i < gammak_.size()-1) s << ", ";
    }
    s << "], ";
    s << "\"Hiso\": " << Hiso_ << ", ";
    s << "\"gtol\": " << newton_tolerance << ", ";
    s << "\"MaxIter\": " << MaxIter_ << ", ";
    s << "\"density\": " << density_;
    s << "}";
  }
  else {
    s << "NonlinearJ2 tag: " << this->getTag()
      << " E: " << E << " nu: " << nu
      << " fy: " << Fy
      << " m(backstresses): " << int(Ck_.size())
      << " tol: " << newton_tolerance << " iters: " << MaxIter_
      << " rho: " << density_ << "\n";
  }
}

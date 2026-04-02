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
// Written: Claudio M. Perez
// Adapted from the MATLAB implementation by Nikolay Velkov.
//
// Developed with FEDEASLab [1].
//
// References:
//
//  [1] Filippou, F.C. (1998)
//     "FEDEASLab: Finite Elements for Design Evaluation and Analysis of Structures"
//  [2] Hartloper, Sousa, Lignos (2021), ASCE, JSE, Vol. 147, Issue 4
//      "Constitutive Modeling of Structural Steels: Nonlinear Isotropic/Kinematic Hardening"
//      Material Model and Its Calibration
//  [3] Suchocki (2022), Acta Mechanica, Vol 233, Issue 1, pp. 83-120
//      On finite element implementation of cyclic elastoplasticity: theory, coding, and exemplary problems
//  [4] Simo, Hughes (1998), Computational Inelasticity, Springer
//
//
// material J2 tag E nu 
//      -F {Fy Fi Fs}  -b {a b} 
//      -C {C1 C2 ...} -g {g1 g2 ...}
//      -YFtol tol -maxIter n -density rho
//
#include "NonlinearJ2.h"
#include <VectorND.h>
#include <Voight.hpp>
#include <MatrixND.h>
#include <Logging.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <cmath>
#include <limits>
#include <cstring>
#include "hardening/FlowStress.h"

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
  E(E_), nu(nu_), fy(Fy_),
  Hiso_(Hiso),
  a_(a), DInf_(DInf), b_(b), QInf_(QInf),
  YFtol_(YFtol*Fy_), MaxIter_(MaxIter),
  density_(density),
  retInitialTangent(Ce), retTangent(C),
  retStress(pres.sig), retStrain(eps)
{
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
NonlinearJ2::updateState()
{

  const double G  = E/2./(1. +    nu);  // Shear modulus
  const double K  = E/3./(1. - 2.*nu);  // Bulk  modulus

  // Deviatoric and volumetric strain
  const VectorND<6> eps_dev = Voight::Dev(eps);
  const double      eps_vol = Voight::Trace(eps);

  // Trial deviatoric stress s_tr = 2G * (eps_d - eps_p_d)
  const VectorND<6> eps_p_d = Voight::Dev(past.eps_p);
  VectorND<9> s_tr{};
  s_tr.addVector(0.0, P^eps_dev,  2.0*G);
  s_tr.addVector(1.0, P^eps_p_d, -2.0*G);

  // Sum of back-stresses at last commit
  VectorND<9> Xn{};
  for (const auto &x : past.sig_b)
    Xn.addVector(1.0, x, 1.0);

  // Trial yield function
  const double sig_nrm = metric(s_tr - Xn);

  const double yf_tr   = sig_nrm - Isotropic::Y(*this, past, 0.0)*phi_n*mises_scale;

  // Elastic admissible?
  if (yf_tr <= YFtol_) {
    // Stress = P*s_tr + K*eps_v*ivol
    pres.sig = P * s_tr;
    pres.sig.addVector(1.0, Voight::ivol,  K*eps_vol);

    pres.eps      = eps;
    pres.eps_p    = past.eps_p;
    pres.gamma    = past.gamma;
    pres.sig_b    = past.sig_b;

    // Elastic tangent
    C.zero();
    C.addMatrix(Voight::IoI,      K);
    C.addMatrix(Voight::IIdevCon, 2.0*G);

    // Energies (use avg stress with Deps)
    if (false) {
      const VectorND<6> sigav = (pres.sig + past.sig) * 0.5;
      const VectorND<6> Deps  = pres.eps - past.eps;
      const VectorND<6> Deps_p = pres.eps_p - past.eps_p;
      pres.Estr = past.Estr + sigav.dot(Deps);
      pres.Ehst = past.Ehst + sigav.dot(Deps_p);
    }
    return 0;
  }

  //
  // Plastic correction (Newton on lamda)
  //
  double lamda   = 0.0;
  double g, Dg, phi_a=0.0;
  int iter = 0;
  VectorND<9> m;

  this->newton_update(s_tr,lamda,  m,g,Dg,phi_a);

  if constexpr (false) {

    // Build a left bracket at lambda=0 for bisection fallback
    double g_lo, Dg_lo, phi_lo;
    VectorND<9> m_lo;
    this->newton_update(s_tr, 0.0, m_lo, g_lo, Dg_lo, phi_lo);
    double lam_lo = 0.0;

    // Right bracket (unknown at start)
    double lam_hi = std::numeric_limits<double>::quiet_NaN();
    double g_hi   = std::numeric_limits<double>::quiet_NaN();
    bool   have_hi = false;

    const int    MaxBacktrack = 6;       // try at most 6 halvings
    const double eta          = 1e-4;    // sufficient decrease: |g_new| <= (1-eta)|g|
    const double EPSgprime    = 1e-14;   // slope floor to avoid divide-by-~0

    while ((std::abs(g) > YFtol_) && (iter < MaxIter_)) {

      // Safeguard slope sign/magnitude (prefer negative; flip tiny values to a usable magnitude)
      if (std::abs(Dg) < EPSgprime) Dg = (Dg >= 0.0 ? -EPSgprime : -EPSgprime); // prefer descent

      // Newton trial (project to λ >= 0)
      double delta = - g / Dg;
      double alpha = 1.0;
      double lam_trial, g_trial, Dg_trial, phi_trial;
      VectorND<9> m_trial;

      bool accepted = false;
      for (int bt = 0; bt <= MaxBacktrack; ++bt) {
        lam_trial = lamda + alpha * delta;
        if (lam_trial < 0.0)
          lam_trial = 0.0;  // projection

        this->newton_update(s_tr, lam_trial, m_trial, g_trial, Dg_trial, phi_trial);

        const bool finite = std::isfinite(g_trial) && std::isfinite(Dg_trial);
        const bool decr   = std::abs(g_trial) <= (1.0 - eta) * std::abs(g);

        if (finite && decr) {
          accepted = true;
          break;  // accept damped step
        }
        alpha *= 0.5; // backtrack
        // opserr << "Backtracking step " << bt+1 << ", alpha = " << alpha << ", |g| = " << std::abs(g_trial) << "\n";
      }

      if (!accepted) {
        opserr << "Newton step rejected\n";
        // Fallback: try a bisection step if we have a right bracket (g changes sign)
        if (!have_hi && g_trial < 0.0) {
          have_hi = true;
          lam_hi = lam_trial;
          g_hi = g_trial;
        }
        if (have_hi) {
          lam_trial = 0.5 * (lam_lo + lam_hi);
          this->newton_update(s_tr, lam_trial, m_trial, g_trial, Dg_trial, phi_trial);
        } else {
          // No bracket yet: take a conservative half step toward zero
          lam_trial = 0.5 * std::max(0.0, lamda + delta);
          this->newton_update(s_tr, lam_trial, m_trial, g_trial, Dg_trial, phi_trial);
        }
      }

      // Maintain bracket (g usually decreases with lambda; track sign)
      if (g_trial < 0.0) {
        have_hi = true;
        lam_hi = lam_trial;
        g_hi = g_trial;
      } else {
        lam_lo  = lam_trial; 
        g_lo = g_trial;
      }

      // Update iterate
      lamda = lam_trial;
      m     = m_trial;
      g     = g_trial;
      Dg    = (std::abs(Dg_trial) >= EPSgprime)
            ? Dg_trial
            : Dg_trial - std::copysign(EPSgprime, Dg_trial);

      ++iter;
    }
  } 
  else {
    while ((std::abs(g) > YFtol_) && (iter < MaxIter_)) {

      Dg = std::abs(Dg)>EPS_DBL ? Dg : ((Dg>=0)?EPS_DBL:-EPS_DBL);
      // Newton step
      lamda -= g / Dg;
      if (lamda < 0.0)
        lamda = 0.0;

      // Recompute Z, n and residual
      this->newton_update(s_tr,lamda,  m,g,Dg,phi_a);

      iter++;
    }
  }

  if (iter == MaxIter_ && std::abs(g) > YFtol_) {
    opserr << "Material failed to converge after ";
    opserr << MaxIter_ << " iterations, |g| = " 
    << std::abs(g) 
    << " > " << YFtol_
    << "\n";
    return -1;
  }

  //
  // Plastic updates at converged lamda
  //

  // eps_p += I^{-1} * (P' * lamda*n)
  const VectorND<6> n6 = P * m;
  // apply I^{-1} (double the shear components)
  VectorND<6> deps_p_Iinv = n6;
  for (int i=3;i<6;i++)
    deps_p_Iinv[i] *= 2.0;

  pres.eps_p = past.eps_p;
  pres.eps_p.addVector(1.0, deps_p_Iinv, lamda);


  // Update isotropic variable
  const double d_mises_strain = lamda*flow_rate*(phi_n*mises_scale);
  Isotropic::update(*this, past, d_mises_strain, pres);

  // Update back-stress components: 
  // X_k = w_k * ( X_k + (2/3) C_k * lamda * n )
  Kinematic::update(*this, past, d_mises_strain, m, pres);

  // Assemble full stress = P*s_dev + K*eps_v*ivol

  // Correct deviatoric stress: s = s_tr - 2G*lamda*n
  pres.sig.addMatrixVector(0.0, P, s_tr, 1.0);
  pres.sig.addMatrixVector(1.0, P, m, -2.0*G*lamda);
  pres.sig.addVector(1.0, Voight::ivol,  K*eps_vol);

  pres.eps = eps;

  //
  // Algorithmic tangent
  //

  // Effective shear modulus factor
  const double theta_phi = lamda*phi_m / phi_a;

  C.zero();
  // volumetric identity scaled by bulk modulus
  C.addMatrix(Voight::IoI,      K);
  // scale contravariant deviatoric identity
  C.addMatrix(Voight::IIdevCon,  2.0*G*(1.0 - 2.0*G*theta_phi));

  const double theta_lam = -2.0*G/(flow_rate*Dg);
  C.addTensorProduct(n6, n6,  -2.0*G*(theta_lam - 2.0*G*theta_phi*phi_nnmm));

  // term with back-stress coupling
  Kinematic::addTangent(*this, past, d_mises_strain, m, -2.0*G*theta_phi*theta_lam, C);


  // Energies
  if (false) {
    const VectorND<6> sigav = (pres.sig + past.sig) * 0.5;
    pres.Estr = past.Estr + sigav.dot(pres.eps - past.eps);
    pres.Ehst = past.Ehst + sigav.dot(pres.eps_p - past.eps_p);
  }

  return 0;
}

//
// NDMaterial API
//

int
NonlinearJ2::setTrialStrain(const Vector &v)
{
  for (int i=0;i<6;i++)
    eps[i] = v(i);
  return updateState();
}

int
NonlinearJ2::setTrialStrainIncr(const Vector &dv)
{
  for (int i=0;i<6;i++)
    eps[i] += dv(i);
  return updateState();
}

const Matrix &NonlinearJ2::getTangent()         { return retTangent; }
const Matrix &NonlinearJ2::getInitialTangent()  { return retInitialTangent; }
const Vector &NonlinearJ2::getStress()          { return retStress; }
const Vector &NonlinearJ2::getStrain()          { return retStrain; }


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
  pres.gamma = 0.0;
  pres.Estr  = 0.0;
  pres.Ehst  = 0.0;

  past = pres;

  // Elastic stiffness
  const double G  = E/2./(1. +    nu);     // Shear modulus
  const double K  = E/3./(1. - 2.*nu);     // Bulk  modulus
  Ce.zero();
  Ce.addMatrix(Voight::IoI,      K);
  Ce.addMatrix(Voight::IIdevCon, 2.0*G);
  C = Ce;

  return 0;
}


NDMaterial *
NonlinearJ2::getCopy()
{
  auto *m = new NonlinearJ2(this->getTag(), E, nu, 
                          fy, density_, 
                          Hiso_,
                          a_, DInf_, b_, QInf_,
                          Ck_, gammak_, YFtol_, MaxIter_);
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


int 
NonlinearJ2::sendSelf(int, Channel &)
{
  // Not implemented
  return -1;
}

int 
NonlinearJ2::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
  // Not implemented
  return -1;
}


void
NonlinearJ2::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": \"NonlinearJ2\", ";
    s << "\"tag\": " << this->getTag() << ", ";
    s << "\"E\": " << E << ", ";
    s << "\"nu\": " << nu << ", ";
    s << "\"fy\": " << fy << ", ";
    s << "\"Q\": [" << QInf_ << ", " << DInf_ << "], ";
    s << "\"b\": [" << b_ << ", " << a_ << "], ";
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
    s << "\"YFtol\": " << YFtol_ << ", ";
    s << "\"MaxIter\": " << MaxIter_ << ", ";
    s << "\"density\": " << density_;
    s << "}";
  } else {
    s << "NonlinearJ2 tag: " << this->getTag()
      << " E: " << E << " nu: " << nu
      << " fy: " << fy
      << " a: " << a_ << " DInf: " << DInf_
      << " b: " << b_ << " QInf: " << QInf_
      << " m(backstresses): " << int(Ck_.size())
      << " tol: " << YFtol_ << " iters: " << MaxIter_
      << " rho: " << density_ << "\n";
  }
}

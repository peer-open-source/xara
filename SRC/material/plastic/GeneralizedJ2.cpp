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
// cmp
//
// References:
//
//  [1] Filippou, F.C. (1998)
//     "FEDEASLab: Finite Elements for Design Evaluation and Analysis of Structures"
//
//  [2] Auricchio, F., and R.L. Taylor. 
//      "Two Material Models for Cyclic Plasticity: Nonlinear Kinematic Hardening and Generalized Plasticity." 
//      International Journal of Plasticity 11, no. 1 (1995): 65–98. 
//      https://doi.org/10.1016/0749-6419(94)00039-5.
//
//
// #include <matrix/routines/poly34.h>
#include <GeneralizedJ2.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <Voigt.hpp>
#include <cstring>
// #include <constmath/constmath.hpp>

#include "GeneralizedJ2.h"
#include <Logging.h>
#include <Channel.h>
#include <cmath>
#include <limits>

using namespace OpenSees;

static const double SQRT23 = std::sqrt(2.0/3.0);


//
// Hardening rules
//

int
GeneralizedJ2::hard_LP(const Hardening &hd, double yf_tr,
                       const TrialState &hs,
                       double &Dlam, double &xi_p, double &theta_phi, double &Tlam) noexcept
{
  const double H  = hd.Hi + hd.Hk;
  const double G  = hd.G;

  Dlam = yf_tr / (2.0*G + H);
  if (Dlam < 0.0)
    Dlam = 0.0;

  xi_p = 2.0*G / (2.0*G + H);
  Tlam = 1.0;
  theta_phi = Dlam/hs.sn;
  return 0;
}


int
GeneralizedJ2::hard_GP(const Hardening &hd,
                       double yf_tr, const TrialState &hs,
                       double &Dlam, double &xi_p, double &theta_phi, double &Tlam) noexcept
{
  const double H     = hd.Hk + hd.Hi;
  const double G     = hd.G;
  const double delta = hd.delta;
  const double phi   = hd.phi;   // phi = Fs - Fy
  const double G1 = G + H/2.0;

  const double A2 = hs.A2;       // = sig_nrm - past.sigdnrm
  const double A1 = yf_tr;
  const double A3 = delta - 2.0*G;
  const double A4 = (delta + H)*phi;

  // Quadratic a x^2 + b x + c = 0
  const double a = 2.0*G1*A3;
  const double b =  A4 - A1*(delta - 2.0*G) + (2.0*G + H)*A2;
  const double c = -A1*A2;

  const double discr = b*b - 4.0*a*c;
  if (discr < 0.0) {
    // Should not happen if inputs are consistent
    Dlam = 0.0;
    xi_p = 1.0;
    Tlam = 1.0;
    return -1;
  }

  const double sqd = std::sqrt(discr);
  const double x1  = (-b + sqd)/(2.0*a);
  const double x2  = (-b - sqd)/(2.0*a);

  // smallest positive root
  Dlam = std::numeric_limits<double>::infinity();
  if (x1 > 0.0 && x1 < Dlam) Dlam = x1;
  if (x2 > 0.0 && x2 < Dlam) Dlam = x2;
  if (!std::isfinite(Dlam)) Dlam = 0.0;

  // algorithmic tangent coefficient xi_p (Adiscr)
  // xi_p = TG*(B1 + B2) / ((TG + H)*B1 - (delta - 2.0*G)*B2 + B3)
  const double B1    = A2 + A3*Dlam;
  const double B2    = A1 - (2.0*G + H)*Dlam; 
  const double B3    = A4; 
  const double Denom = (2.0*G + H)*B1 - A3*B2 + B3;

  // Adiscr
  xi_p  = (Denom != 0.0) ? (2.0*G*(B1 + B2)/Denom) : 1.0;
  Tlam  = 1.0;
  theta_phi = Dlam/(hs.sn);
  return 0;
}


int
GeneralizedJ2::hard_AF(const Hardening &hd,
                       double yf_tr, 
                       const TrialState &hs,
                       double &Dlam, double &xi_p, 
                       double & theta_phi,
                       double &Tlam,
                       VectorND<9> &M, // flow direction
                       double & theta_nlk
                      ) noexcept
{
  return -1;
}

//
// Constructors / basic
//

GeneralizedJ2::GeneralizedJ2(int tag, 
                             double E_, 
                             double nu_, 
                             double Fy_,
                             double Fs_, 
                             double Hiso_, 
                             double Hkin_, 
                             double delta_,
                             double density_,
                             HRule rule_)
: NDMaterial(tag, 0),
  E(E_), nu(nu_), Fy(Fy_), Fs(Fs_),
  Hiso(Hiso_), Hkin(Hkin_), delta(delta_),
  rule(rule_),
  retTangent(C), retInitialTangent(Ce), 
  retStress(pres.sig), retStrain(eps),
  density(density_)
{
  this->revertToStart();
  this->commitState();
}


GeneralizedJ2::~GeneralizedJ2() = default;

//
// State update/return mapping
//

int
GeneralizedJ2::updateState(const VectorND<6> &eps)
{
  static constexpr MatrixND<6,9> P {
    // NOTE: this appears transposed because MatrixND is column-major
    1.0000,        0,        0,        0,        0,        0,
         0,   1.0000,        0,        0,        0,        0,
         0,        0,   1.0000,        0,        0,        0,
         0,        0,        0,   0.5000,        0,        0,
         0,        0,        0,   0.5000,        0,        0,
         0,        0,        0,        0,   0.5000,        0,
         0,        0,        0,        0,   0.5000,        0,
         0,        0,        0,        0,        0,   0.5000,
         0,        0,        0,        0,        0,   0.5000};

  // Elastic moduli and projections
  const double G  = E/2./(1. +    nu);     // Shear modulus
  const double K  = E/3./(1. - 2.*nu);     // Bulk  modulus
  const double TG  = 2.0*G;

  const double norm_unit = 1.0;
  const double flow_rate = SQRT23; // 1.0 for J2, sqrt(3/2) for Mises


  // deviatoric part of total strain
  const VectorND<6> eps_dev = Voigt::Dev(eps);
  const double      eps_vol = Voigt::Trace(eps); // volumetric strain

  // Trial deviatoric stress: s_tr = 2G * (eps_dev - eps_p_d)
  const VectorND<6> eps_p_d = Voigt::Dev(past.eps_p);
  VectorND<9> sig_dev_tr{};
  sig_dev_tr.addVector(0.0, P^eps_dev,  2.0*G);
  sig_dev_tr.addVector(1.0, P^eps_p_d, -2.0*G);

  // Yield function at trial
  const VectorND<9> eta = sig_dev_tr - past.sig_b;
  const double sig_nrm = eta.norm()*norm_unit;
  const double yf_tr   = sig_nrm - flow_rate*(Fy + Hiso*past.alpha);
  const double A2      = sig_nrm - past.sigdnrm;

  // Elastic (admissible) or plastic
  bool elastic = (yf_tr <= 0.0) || (A2 < 0.0);

  if (elastic) {
    // Elastic update
    // Full Cauchy stress = deviatoric + K*eps_v*IVOL

    // VectorND<6> sig = P*sig_dev_tr;
    // sig.addVector(1.0, Voigt::ivol, K*eps_vol);
    pres.sig = Voigt::ReduceVector(sig_dev_tr);
    Voigt::AddVol(pres.sig, K*eps_vol);

    // pres.sig     = sig;
    // pres.eps     = eps;
    pres.eps_p   = past.eps_p;
    pres.sig_b   = past.sig_b;
    pres.alpha   = past.alpha;
    pres.sigdnrm = (sig_dev_tr - pres.sig_b).norm()*norm_unit;

    // Tangent
    C.zero();
    C.addMatrix(Voigt::IoI,          K);
    C.addMatrix(Voigt::IIdevCon, 2.0*G);
    return 0;
  }

  //
  // Plastic correction
  //
  // Flow direction (stress space) n = eta / ||eta||
  VectorND<9> n{};
  if (sig_nrm > 0.0) {
    n = eta;
    n /= sig_nrm;
  }

  // Hardening rule selection
  const double phi = Fs - Fy; // phi is distance to asymptote
  const Hardening hd{Hiso, Hkin, G, delta, phi, rule};
  TrialState hs{};
  hs.A2    = A2;
  hs.sn    = sig_nrm;
  hs.se    = sig_dev_tr;
  hs.sb    = past.sig_b;
  hs.alpha = past.alpha;
  hs.flow_stress = (Fy + Hiso*past.alpha)*flow_rate;

  double Dlam=0.0, xi_p=1.0, Tlam=1.0;
  double theta_phi, theta_nlk=0;
  int status = -1;
  switch (rule) {
    case HRule::LP: status = hard_LP(hd, yf_tr, hs,          Dlam, xi_p, theta_phi, Tlam); break;
    case HRule::GP: status = hard_GP(hd, yf_tr, hs,          Dlam, xi_p, theta_phi, Tlam); break;
    case HRule::AF: status = hard_AF(hd, yf_tr, hs,          Dlam, xi_p, theta_phi, Tlam, n, theta_nlk); break;
    default:        status = hard_GP(hd, yf_tr, hs,          Dlam, xi_p, theta_phi, Tlam); break;
  }
  if (status != 0) {
    opserr << "Hardening rule failed; falling back to elastic update\n";
    // Something went wrong in the hardening; fall back to elastic
    pres = past;
    C.zero();
    C.addMatrix(Voigt::IoI,          K);
    C.addMatrix(Voigt::IIdevCon, 2.0*G);
    return -1;
  }

  // Plastic strain increment:
  // eps_p += I\(P'*(Dlam*n9))
  VectorND<6> deps = P*n;
  for (int i=3;i<6;i++)
    deps[i] *= 2.0; // I^{-1} effect
  pres.eps_p = past.eps_p;
  pres.eps_p.addVector(1.0, deps, Dlam);

  // Correct deviatoric stress: s = s_tr - 2G * Dlam * n
  VectorND<9> sig_dev = sig_dev_tr;
  sig_dev.addVector(1.0, n, -2.0*G*Dlam);

  // Update isotropic variable alpha
  pres.alpha = past.alpha + flow_rate*Dlam;

  // Update back-stress: b = (b + (2/3)Hk * Dlam * n) * Tlam
  pres.sig_b = past.sig_b;
  pres.sig_b.addVector(1.0, n, (2.0/3.0)*Hkin*Dlam);
  pres.sig_b *= Tlam;

  // Assemble full stress
  pres.sig = P*sig_dev;
  pres.sig.addVector(1.0,  Voigt::ivol,  K*eps_vol);

  // Store norm for next step
  pres.sigdnrm = (sig_dev - pres.sig_b).norm();

  //
  // Algorithmic tangent 
  //
  // (Auricchio–Taylor J2)
  // Ct = K*IoI + TG*(1 - xi_p)*(n o n) + TG*(1 - TG*Dlam/sig_nrm)*(Id - n o n)
  // Ct = K*IoI + TG*(1 - TG*Dlam/sig_nrm)*Id + TG*(1 - xi_p - (1 - TG*Dlam/sig_nrm))*(n o n)
  //    = K*IoI + TG*(1 - TG*Dlam/sig_nrm)*Id + TG*(TG*Dlam/sig_nrm - xi_p)*n o n
  //                      |       C      |               C             Ad
  C.zero();
  // 1. bulk term
  C.addMatrix(Voigt::IoI, K);

  // IIdevCon := IIcon - IoI/3
#if 1 // Original implementation
  const double a = 1.0 - TG*Dlam/sig_nrm;
  const double b = (1.0 - xi_p);

  // term: TG * a * (Id - n o n)
  MatrixND<6,6> nn{};
  nn.addTensorProduct(P*n,P*n, 1.0);
  MatrixND<6,6> Pnn = Voigt::IIdevCon;
  Pnn.addMatrix(nn, -1.0);
  C.addMatrix(Pnn, 2.0*G*a);

  // term: TG * b * (n o n)
  C.addMatrix(nn, 2.0*G*b);
#else
  const double a = 1.0 - TG*Dlam/sig_nrm;
  const double b = (1.0 - xi_p);

  // 2. scale contravariant deviatoric identity
  C.addMatrix(Voigt::IIdevCon,  2.0*G*(1.0 - 2.0*G*theta_phi));

  const double theta_lam = xi_p;
  const VectorND<6> Pn = P*n;
  C.addTensorProduct(Pn, Pn,  -2.0*G*(theta_lam - 2.0*G*theta_phi));

  if (theta_nlk != 0) {
    // TODO
  }

#endif
  return 0;
}

//
// NDMaterial API
//

int
GeneralizedJ2::setTrialStrain(const Vector &v)
{
  assert(v.Size() == 6);
  // copy into internal Voigt
  for (int i=0;i<6;i++)
    eps[i] = v(i);
  return updateState(eps);
}


int
GeneralizedJ2::setTrialStrainIncr(const Vector &dv)
{
  for (int i=0;i<6;i++)
    eps[i] += dv(i);
  return updateState(eps);
}

const Matrix&
GeneralizedJ2::getTangent()
{
  return retTangent;
}

const Matrix&
GeneralizedJ2::getInitialTangent()
{
  return retInitialTangent;
}

const Vector&
GeneralizedJ2::getStress()
{
  return retStress;
}

const Vector&
GeneralizedJ2::getStrain()
{
  return retStrain;
}

int
GeneralizedJ2::commitState()
{
  past = pres;
  return 0;
}

int
GeneralizedJ2::revertToLastCommit()
{
  pres = past;
  C = Ce;
  return 0;
}

int
GeneralizedJ2::revertToStart()
{
  // zero states
  eps.zero();
  pres.sig.zero();
  pres.eps_p.zero();
  pres.sig_b.zero();
  pres.sigdnrm = 0.0;
  pres.alpha   = 0.0;

  past = pres;

  // Build elastic matrices
  const double G  = E/2./(1. +    nu);     // Shear modulus
  const double K  = E/3./(1. - 2.*nu);     // Bulk  modulus

  Ce.zero();
  Ce.addMatrix(Voigt::IoI,          K);
  Ce.addMatrix(Voigt::IIdevCon, 2.0*G); // Idp
  C = Ce;
  return 0;
}


NDMaterial *
GeneralizedJ2::getCopy()
{
  auto *m = new GeneralizedJ2(this->getTag(), E, nu, Fy, Fs, Hiso, Hkin, delta, density, rule);
  // copy state
  m->eps       = this->eps;
  m->past      = this->past;
  m->pres      = this->pres;
  return m;
}


NDMaterial *
GeneralizedJ2::getCopy(const char *type)
{
  if (type && std::strcmp(type, "ThreeDimensional") == 0)
    return this->getCopy();

  return NDMaterial::getCopy(type);
}


int
GeneralizedJ2::sendSelf(int ctag, Channel &theChannel)
{
  return 0;
}

int
GeneralizedJ2::recvSelf(int ctag, Channel &theChannel, FEM_ObjectBroker &)
{
  return 0;
}

void
GeneralizedJ2::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"GeneralizedJ2\", ";
    s << "\"E\": " << E << ", ";
    s << "\"nu\": " << nu << ", ";
    s << "\"Fy\": " << Fy << ", ";
    s << "\"Fs\": " << Fs << ", ";
    s << "\"Hiso\": " << Hiso << ", ";
    s << "\"Hkin\": " << Hkin << ", ";
    s << "\"delta\": " << delta << ", ";
    s << "\"rule\": \"" << (rule==HRule::LP ? "LP" : "GP") << "\", ";
    s << "\"density\": " << density;
    s << "}";
    return;
  }
  else {
    s << "GeneralizedJ2 tag: " << this->getTag()
      << " E: " << E << " nu: " << nu
      << " Fy: " << Fy << " Fs: " << Fs
      << " Hiso: " << Hiso << " Hkin: " << Hkin
      << " delta: " << delta
      << " rule: " << (rule==HRule::LP ? "LP" : "GP")
      << " rho: " << density << "\n";
  }
}

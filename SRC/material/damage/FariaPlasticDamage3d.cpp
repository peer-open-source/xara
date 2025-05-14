/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Written in Matlab: Thanh Do
// Created: 07/16
//
#include "FariaPlasticDamage3d.h"
#include <Channel.h>
#include <MatrixND.h>
#include <cmath>
#include <Projector.hh>
#include <elementAPI.h>

#define MND
using namespace OpenSees;
// static Matrix Idp(6,6); 
// static Matrix I(6,6);
// static Matrix Id(6,6);


void * OPS_ADD_RUNTIME_VPV(OPS_NewPlasticDamageConcrete3d)
{
  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 5 || numArgs > 9) {
    opserr << "Want: nDMaterial PlasticDamageConcrete3d $tag $E $nu $ft $fc <$beta $Ap $An $Bn>\n";
    return 0;	
  }
  
  int iData[1];
  double dData[8];
  dData[4] = 0.6;  // beta
  dData[5] = 0.5;  // Ap
  dData[6] = 2.0;  // An
  dData[7] = 0.75; // Bn
  
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return 0;
  }
  
  numData = numArgs - 1;;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: " << iData[0] <<"\n";
    return 0;
  }  
  
  theMaterial = new PlasticDamageConcrete3d(iData[0], 
					    dData[0], dData[1], dData[2],dData[3], 
					    dData[4], dData[5], dData[6],dData[7]);

  return theMaterial;
}


PlasticDamageConcrete3d::PlasticDamageConcrete3d(int tag,
                                                  double E, 
                                                  double nu, 
                                                  double _ft,
                                                  double _fc, 
                                                  double _beta, 
                                                  double _Ap,
                                                  double _An, 
                                                  double _Bn)
  :NDMaterial(tag,ND_TAG_PlasticDamageConcrete3d),
 E(E), nu(nu), ft(_ft), fc(_fc), beta(_beta), Ap(_Ap), An(_An), Bn(_Bn)
 // eps(6), sig(6), sige(6), eps_p(6), sigeP(6),
 //epsCommit(6), sigCommit(6), sigeCommit(6), eps_pCommit(6) //, sigePCommit(6)
// , Ce(6,6), C(6,6), Ccommit(6,6)
{
  eps.zero();
  sig.zero();
  sige.zero();
  eps_p.zero();
  // Ce.zero();

  // additional material parameters
  double G  = E/2./(1+nu);       // shear modulus
  double K  = E/3./(1-2*nu);     // bulk  modulus

  Ce.zero();
  Ce.addMatrix(IIvol, K);
  Ce.addMatrix(IIdevCon, 2.*G);
  
  C = Ce;
  // Ccommit.zero();
  
  double f2c = 1.16*fc;
  double k = sqrt(2.0)*(f2c - fc)/(2.*f2c - fc);

  // initial damage threshold
  double rp0 = ft/sqrt(E);
  double rn0 = sqrt((-k+sqrt(2.0))*fc/sqrt(3.0));
      
  rp = rp0;
  rn = rn0;
  dp = 0.;
  dn = 0.;

  this->commitState();
}

PlasticDamageConcrete3d::PlasticDamageConcrete3d()
  :NDMaterial (0, ND_TAG_PlasticDamageConcrete3d)
  //  eps(6), sig(6), sige(6), eps_p(6), sigeP(6),
  //  epsCommit(6), sigCommit(6), sigeCommit(6), eps_pCommit(6)//, sigePCommit(6)
  //  Ce(6,6), C(6,6), Ccommit(6,6)
{

}

PlasticDamageConcrete3d::~PlasticDamageConcrete3d()
{

}

/*
%% Stress invariant function: octahedral normal and shear stresses
function [sigoct, tauoct] = StrsInvar (sig)

% normal stress
sigoct = (sig(1) + sig(2) + sig(3))/3;

% shear stress
J2 = ((sig(1) - sig(2))^2 + (sig(1) - sig(3))^2 + (sig(2) - sig(3))^2)/6 +...
     (sig(4))^2 + (sig(5))^2 + (sig(6))^2;
tauoct = (2/3*J2)^0.5;
end
*/

static inline void
StrsInvar(const Vector &sig, double &sigoct, double &tauoct)
{
  // normal stress
  sigoct = (sig(0) + sig(1) + sig(2))/3.;
  // shear stress
  double J2 = (pow((sig(0) - sig(1)),2) + pow((sig(0) - sig(2)),2) + pow((sig(1) - sig(2)),2))/6. + 
               pow(sig(3),2) + pow(sig(4),2) + pow(sig(5),2);
  tauoct = sqrt(2./3.*J2);  
}

/*
%% Stress decomposition function: algebraic approach
function [sigpos, signeg, Qpos, Qneg] = StrsDecA (sig)

% positive and negative stress tensors
sigpos = (sig + abs(sig))/2;
signeg = sig - sigpos;

% projection tensors
Qpos = diag((1+sign(sig))/2);
Qneg = eye(6) - Qpos;
end
*/
template <typename T> static inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
template <typename V6, typename M6>
static inline void
StrsDecA(const Vector &sig, V6 &sigpos, V6 &signeg, 
                            M6 &Qpos,   M6 &Qneg)
{
  for (int i=0; i<6; i++) {
    sigpos(i) = (sig(i) + std::fabs(sig(i))) / 2.;
    signeg(i) = sig(i) - sigpos(i);
    Qpos(i,i) = (1. + sgn(sig(i))) / 2.;
    Qneg(i,i) = 1. - Qpos(i,i);
  }
  return;
  for (int i=0; i<6; i++) {
    if (sig(i) > 1e-8) {
      sigpos(i) = sig(i);
      signeg(i) = 0.;
      Qpos(i,i) = 1;
      Qneg(i,i) = 0;
    } else if (sig(i)  < -1.e-8) {
      sigpos(i) = 0.;
      signeg(i) = sig(i);
      Qpos(i,i) = 0;
      Qneg(i,i) = 1;
    } else {
      sigpos(i) = sig(i)/2.;
      signeg(i) = sig(i)/2.;
      Qpos(i,i) = 0.5;
      Qneg(i,i) = 0.5;
    }
  }
}

int
PlasticDamageConcrete3d::setTrialStrain(Vector const&v1, Vector const&v2)
{
  return this->setTrialStrain(v1);
}


int
PlasticDamageConcrete3d::setTrialStrain(const Vector &strain)
{

  double f2c = 1.16*fc; // f2c: biaxial compressive strength
  double k = sqrt(2.0)*(f2c - fc)/(2.*f2c - fc);
  // initial damage threshold
  double rp0 = ft/sqrt(E);
  double rn0 = sqrt((-k+sqrt(2.0))*fc/sqrt(3.0));

  constexpr static double tol = 1.0e-5;

  // retrieve history variables
  eps_p = eps_pCommit;
  // sigeP = sigePCommit;
  VectorND<6> sigeP = sigeCommit;
  rp = rpCommit;
  rn = rnCommit;
  dp = dpCommit;
  dn = dnCommit;

  // current strain
  eps = strain;

  // static VectorND<6> Depse_tr; // (6);
  // static VectorND<6> Deps; // (6);
  // incremental strain
  VectorND<6> Depse_tr = eps - eps_p;
  VectorND<6> Deps = eps - epsCommit;

  //
  // PLASTIC part
  //

  // elastic trial stress
  // static VectorND<6> sige_tr;//(6);
  VectorND<6> sige_tr  = sigeCommit;
  sige_tr += Ce*Deps;

  // decomposition of trial stress tensor
  VectorND<6> sigpos{}; // (6);
  VectorND<6> signeg{}; // (6);
#ifdef MND
  MatrixND<6,6> Qpos{}, Qneg{};
#else 
  Matrix Qpos(6,6), Qneg(6,6);
  Qpos.Zero();
  Qneg.Zero();
#endif
  StrsDecA(sige_tr, sigpos, signeg, Qpos, Qneg);

  // compute stress invariants of the negative stress tensor
  double sigoct, tauoct;
  StrsInvar(signeg, sigoct, tauoct);


  // check tauneg
  double taun = std::sqrt(std::sqrt(3.0)*(k*sigoct + tauoct)); // negative equivalent stress

  // Correction
  
  MatrixND<6,6> Cbar{};
  if ((taun - rn) <= (tol*rn0)) { // elastic state, accept trial response
    sige = sige_tr;                                                    
    Cbar = Ce;                                                          
  }
  else {
    //  norm of trial effective stress
    double nrm = sqrt(pow(sige_tr(0),2) + pow(sige_tr(1),2) + pow(sige_tr(2),2) 
               + 2*pow(sige_tr(3),2) + 2*pow(sige_tr(4),2) + 2*pow(sige_tr(5),2));

    // static Vector L_tr(6);
    VectorND<6> L_tr = sige_tr; 
    L_tr/=nrm;    // normalized trial effective stress

    //  Deps_p = beta*E*(L_tr'*Deps)*Depse_tr/nrm;     // plastic strain increment
    double L_trDotDeps = L_tr.dot(Deps);

    // static Vector Deps_p(6);
    VectorND<6> Deps_p = Depse_tr;
    Deps_p *= beta*E*L_trDotDeps/nrm;

    double lam  = 1 - beta*E/nrm * L_trDotDeps;      //  scale factor

    sige  = sige_tr; 
    sige *= lam;                   //  corrected effective stress

    // check damage
    StrsDecA(sige, sigpos, signeg, Qpos, Qneg);         //  decompose the effective stress        
    StrsInvar(signeg, sigoct, tauoct);                  //  find octahedral stresses
    taun = sqrt(sqrt(3.0)*(k*sigoct + tauoct));          //  negative equivalent stress


    if ((taun - rn <= tol*rn0) || (L_trDotDeps <= 0)) { //  no damage or sige and eps in different direction
      sige = sige_tr;
      Cbar = Ce;
    }

    else {
      eps_p  = eps_p + Deps_p;                          //  update plastic strain

      VectorND<6> L_tr_temp{}; // (6);
      // tangent in effective space, Cbar
      //  L_tr_temp = [L_tr(1:3); 2*L_tr(4:6)];
      for (int i=0; i<3; i++) 
        L_tr_temp(i) = L_tr(i);
      for (int i=3; i<6; i++)
        L_tr_temp(i) = 2*L_tr(i);

      double Dlam_Dnrm = 2*beta*E/pow(nrm,3)*sige_tr.dot(Deps);

      // static VectorND<6> Dnrm_Dsig; // (6);
      // static VectorND<6> Dlam_Dsig; // (6);
      VectorND<6> Dnrm_Dsig = L_tr_temp;
      VectorND<6> Dlam_Dsig = Deps; 
      Dlam_Dsig *= -beta*E/(nrm*nrm); // Dlam_Dsig = -beta*E/(nrm*nrm)*Deps;
      // static Vector Dlam_Deps(6);
      VectorND<6> Dlam_Deps = L_tr;
      Dlam_Deps *= -beta*E/nrm;       // Dlam_Deps = -beta*E/nrm*L_tr;
      // Dlam_Deps = Dlam_Dnrm * Ce * Dnrm_Dsig + Ce*Dlam_Dsig + Dlam_Deps; 
      Dlam_Deps += Dlam_Dnrm * (Ce * Dnrm_Dsig);
      Dlam_Deps.addMatrixVector(1.0, Ce, Dlam_Dsig, 1.0);

      Cbar.addMatrix(Ce, lam);
      Cbar.addTensorProduct(sige_tr, Dlam_Deps, 1.0);
    }
  }


  //
  // DAMAGE part
  //
  Qpos.zero();
  Qneg.zero();
  // decompose into positive and negative effective stress tensor
  StrsDecA(sige, sigpos, signeg, Qpos, Qneg);    // decompose the effective stress  

  // calculate equivalent stresses
  VectorND<6> tmp;
  Ce.solve(sigpos, tmp);
  double taup = sqrt(tmp.dot(sigpos));                // positive equivalent stress

  StrsInvar(signeg, sigoct, tauoct);             // find octahedral stresses
  taun = sqrt((sqrt(3.)*(k*sigoct + tauoct)));   // negative equivalent stress

  double Ddp_Drp = 0.;

  // positive damage
  if ((taup - rp) <= (tol*rp0)) {                // no positive damage
    Ddp_Drp = 0;

  }
  else {                                       // positive damage evolves
    rp = taup;                                   // update rp = max(taup, rp)
    dp = 1 - rp0/rp * exp(Ap*(1 - rp/rp0));

    Ddp_Drp =  (Ap*rp + rp0)/(rp*rp) * exp(Ap*(1 - rp/rp0));               
    dp = dp*(1-tol);                             // cap the damage variable 
    Ddp_Drp = Ddp_Drp*(1-tol);        
    if (dp > 1-tol) {
      dp = 1- tol; Ddp_Drp = 0;
    }
  }

  // negative damage
  double Ddn_Drn = 0;
  if (taun - rn <= tol*rn0) {                    // no negative damage
    Ddn_Drn = 0;
  }  else {                                      // negative damage evolves
    rn = taun;                                   // update rn
    dn = 1 - rn0/rn*(1-An) - An*exp(Bn*(1 - rn/rn0));
    Ddn_Drn = rn0/(rn*rn)*(1-An) + An*Bn/rn0*exp(Bn*(1 - rn/rn0));
    dn = dn*(1-tol);                             // cap the damage variable
    Ddn_Drn = Ddn_Drn*(1-tol);
    if (dn > 1-tol) {
      dn = 1- tol; Ddn_Drn = 0;
    }
  }

  // stress update
  sig = (1-dp)*sigpos + (1-dn)*signeg;

  //
  // TANGENT
  //
#ifdef MND
  MatrixND<6,6> Dsigpos_Deps{}, Dsigneg_Deps{};
#else
  static Matrix Dsigpos_Deps(6,6), Dsigneg_Deps(6,6);
#endif
  Dsigpos_Deps = Qpos*Cbar;
  Dsigneg_Deps = Qneg*Cbar;
  
  // static Vector s(6);
  VectorND<6> s = IIdevMix*signeg;            // deviatoric stress


  // norm of deviatoric stress
  double nrms = sqrt( pow(s(0),2) + pow(s(1),2) + pow(s(2),2) +  
                    2*pow(s(3),2) + 2*pow(s(4),2) + 2*pow(s(5),2));


  static VectorND<6> n = s; // (6);
  if (nrms <= tol) 
    n.zero(); 
  else {
    n/=nrms;
  }

  VectorND<6> Dtaup_Dsigpos;
  VectorND<6> Dtaun_Dsigneg;

  if (taup <= tol) {
    Dtaup_Dsigpos.zero();   //  Dtaup_Dsigpos = zeros(6,1); 
  }
  else  {
    Ce.solve(sigpos, Dtaup_Dsigpos);  //  Dtaup_Dsigpos = (Ce\sigpos)/taup; end
    Dtaup_Dsigpos/=taup;
  }

  if (taun <= tol) {
    Dtaun_Dsigneg.zero();
  }
  else {
    double Dtaun_Dsigoct = pow(3,0.25) * k/2./sqrt(k*sigoct + tauoct);
    double Dtaun_Dtauoct = pow(3,0.25) / 2./sqrt(k*sigoct + tauoct);
    // static Vector Dsigoct_Dsigneg(6);
    VectorND<6> Dsigoct_Dsigneg = OpenSees::ivol;
    Dsigoct_Dsigneg /= 3.;
    // static Vector Dtauoct_Dsigneg(6);
    VectorND<6> Dtauoct_Dsigneg = n; Dtauoct_Dsigneg /= sqrt(3.0);
    Dtaun_Dsigneg = Dtaun_Dsigoct * Dsigoct_Dsigneg + Dtaun_Dtauoct * Dtauoct_Dsigneg;
  }

  // static Vector Dnrm_Deps(6);
  // static Vector Ddp_Deps(6);
  // static Vector Ddn_Deps(6);
  // Ddp_Deps = Ddp_Drp * Dsigpos_Deps' * Dtaup_Dsigpos;
  VectorND<6> Ddp_Deps = Dsigpos_Deps ^ Dtaup_Dsigpos;
  Ddp_Deps *= Ddp_Drp; 

  VectorND<6> Ddn_Deps = Dsigneg_Deps ^ Dtaun_Dsigneg;
  Ddn_Deps *= Ddn_Drn;
  
  C.zero();
  C.addMatrix(Dsigpos_Deps , (1-dp));
  C.addMatrix(Dsigneg_Deps , (1-dn));
  C.addTensorProduct(sigpos, Ddp_Deps, -1.0);
  C.addTensorProduct(signeg, Ddn_Deps, -1.0);
    // - (sigpos % Ddp_Deps) 
    // - (signeg % Ddn_Deps);  

  return 0;
}

int
PlasticDamageConcrete3d::setTrialStrainIncr(const Vector &strain)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

int
PlasticDamageConcrete3d::setTrialStrainIncr(const Vector &strain, const Vector &rate)
{
  eps += strain;
  this->setTrialStrain(eps);
  return 0;
}

const Matrix&
PlasticDamageConcrete3d::getTangent()
{
  static Matrix Wrapper(C);
  return Wrapper;
}

const Matrix&
PlasticDamageConcrete3d::getInitialTangent()
{
  static Matrix Wrapper(Ce);
  return Wrapper;
}

const Vector&
PlasticDamageConcrete3d::getStress()
{
  static Vector wrapper(sig);
  return wrapper;
}

const Vector&
PlasticDamageConcrete3d::getStrain()
{
  static Vector wrapper(eps);
  return wrapper;
}

int
PlasticDamageConcrete3d::commitState()
{
  rpCommit    = rp;
  rnCommit    = rn;
  dpCommit    = dp;
  dnCommit    = dn;
  epsCommit   = eps;
  sigCommit   = sig;
  sigeCommit  = sige;
  eps_pCommit = eps_p;
  // sigePCommit = sige;
  sigeCommit = sige;
  Ccommit = C;
  return 0;
}

int
PlasticDamageConcrete3d::revertToLastCommit()
{
  C  = Ccommit;
  rp = rpCommit;
  rn = rnCommit;
  dp = dpCommit;
  dn = dnCommit;
  eps = epsCommit;
  sig = sigCommit;
  sige = sigeCommit;
  eps_p = eps_pCommit;
  // sigeP = sigeCommit;
  return 0;
}

int
PlasticDamageConcrete3d::revertToStart()
{
  eps.zero();
  sig.zero();
  sige.zero();
  eps_p.zero();
  // sigeP.zero();
  double G  = E/2./(1+nu);       // shear modulus
  double K  = E/3./(1-2*nu);     // bulk  modulus
  Ce.zero();
  Ce.addMatrix(IIvol, K);
  Ce.addMatrix(IIdevCon, 2.*G);
  return 0;
}

NDMaterial*
PlasticDamageConcrete3d::getCopy(const char *type)
{
  if (strcmp(type,"ThreeDimensional") == 0 || strcmp(type,"3D") == 0) {
    PlasticDamageConcrete3d *theCopy =
      new PlasticDamageConcrete3d (this->getTag(), E, nu, ft, fc, beta, Ap, An, Bn);
  
    return theCopy;  
  } else {
    return NDMaterial::getCopy(type);
  }
}

NDMaterial*
PlasticDamageConcrete3d::getCopy()
{
  PlasticDamageConcrete3d *theCopy =
    new PlasticDamageConcrete3d (this->getTag(), E, nu, ft, fc, beta, Ap, An, Bn);
  
  return theCopy;
}

const char*
PlasticDamageConcrete3d::getType() const
{
  return "ThreeDimensional";
}

int
PlasticDamageConcrete3d::getOrder() const
{
  return 6;
}


int 
PlasticDamageConcrete3d::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(10);

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PlasticDamageConcrete3d::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
PlasticDamageConcrete3d::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PlasticDamageConcrete3d::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

void 
PlasticDamageConcrete3d::Print(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"type\": \"" << this->getClassType() << "\", ";
    s << "\"tag\": " << this->getTag() << ", ";
    s << "\"E\": " << E << ", ";
    s << "\"nu\": " << nu << ", ";
    s << "\"ft\": " << ft << ", ";
    s << "\"fc\": " << fc << ", ";
    s << "\"beta\": " << beta << ", ";
    s << "\"Ap\": " << Ap << ", ";
    s << "\"An\": " << An << ", ";
    s << "\"Bn\": " << Bn;
    s << "}";
  }
  else {
    opserr << "PlasticDamageConcrete3d: " << this->getTag() << "\n";
  }
}       

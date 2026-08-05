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
// - Isotropic hardening is Chaboche-Rousselier type and consits of a linear 
//   term and two Voce exponential terms.
// - Kinematic hardening employs a sum of Armstrong-Frederick back-stress terms,
//   as developed by Chaboche. With one term, this reproduces the linear Prager
//   type kinematic hardening when gamma=0.
//
// This model is similar to UVCmultiaxial [2], with the following distinctions:
// 1) UVCmultiaxial does not include the linear isotropic hardening term. This 
//    was incorporated in the present model to allow it to replicate standard
//    linearly hardening J2 plasticity when the Voce parameters are set to zero.
//    It is also able to reproduce the common "linear + one exponential" plasticity
//    model that is offered in, eg, FEAP, and employed by [4].
// 2) UVCmultiaxial [2] derives the consistency condition assuming the flow direction
//    coincides with the direction of the trial state. This is only valid for 
//    proportional loading when nonlinear kinematic hardening is employed. 
//    The present model does not make this assumption.
// 3) UVCmultiaxial [2] employs an exponential integration scheme for the back-stress
//    evolution that is derived by assuming the flow direction is constant over the increment.
//    This appears to be a reasonable procedure, and may be reproduced in the present model 
//    by setting:
//             bs_integration =  BackStressIntegration::Exponential;
//    However, the present model also supports a backward-Euler integration scheme
//    for the back-stress evolution which is more commonly employed in the literature. 
// 4) In the original implementation of UVCmultiaxial [2], a missplaced factor of sqrt(2/3)
//    in the kinematic hardening modulus produces an inconsistent tangent. This is not
//    the case in the present model, which produces a consistent tangent.
// 5) The original implementation of UVCmultiaxial [2] symmetrizes the tangent, which
//    becomes asymmetric when nonlinear kinematic hardening is employed. 
//    The present model does not symmetrize the tangent, and produces a consistent tangent.
//
// This model is also similar to that of [3], which supports equivalent hardening rules.
// In [3], backward-Euler integration is employed for the back-stress evolution (item 3 above). 
//
// The conventions adopted (scaling the yield surface, consistency parameter, flow direction, etc) 
// are more consistent with [2] and [4] than with [3]. These conventions are represented by
// compile-time parameters, but currently are:
// - The yield function takes the Frobenius norm. Alternatives would be the J2 invariant
// - The flow direction m is normalized to unit length.
// - The scalar unknown lambda solved for in the state determination is the 
//   Frobenius norm of the plastic strain increment
// 
// TODO:
// - Refactor Isotropic and Kinematic hardening classes, remove from header.
// - Store the back-stress components in 6D if possible, minimize use of 9D representation.
// - Support alternative nonlinear isotropic hardening rules:
//   - Mroz and Maciejewski
//
// References:
//
//  [1] Filippou, F.C. (1998)
//     "FEDEASLab: Finite Elements for Design Evaluation and Analysis of Structures"
//  [2] Hartloper, Sousa, Lignos (2021), ASCE, JSE, Vol. 147, Issue 4
//      "Constitutive Modeling of Structural Steels: Nonlinear Isotropic/Kinematic Hardening"
//      Material Model and Its Calibration
//  [3] Suchocki (2022), Acta Mechanica, Vol 233, Issue 1, pp. 83-120
//      "On finite element implementation of cyclic elastoplasticity: 
//      theory, coding, and exemplary problems"
//  [4] Simo, Hughes (1998), Computational Inelasticity, Springer
//
//===----------------------------------------------------------------------===//
//
// Written: Claudio M. Perez
//
#pragma once
#include <NDMaterial.h>
#include <Voight.hpp>
#include <MatrixND.h>
#include <VectorND.h>
#include <vector>
#include <string>

class Vector;
class Matrix;
class Channel;
class FEM_ObjectBroker;
class OPS_Stream;
//
namespace OpenSees {


class NonlinearJ2 : public NDMaterial
{
public:
  NonlinearJ2(int tag,
              double E, double nu, 
              double fy, double density,
              double Hiso,
              double a, double DInf,     // iso: reduction branch
              double b, double QInf,     // iso: cyclic hardening branch
              const std::vector<double> &Ck,        // kin: magnitudes
              const std::vector<double> &gammak,    // kin: saturation rates
              double YFtol = 1e-16,
              int    MaxIter = 15);

  ~NonlinearJ2();

  const char *getClassType() const override {return "NonlinearJ2";}
  const char *getType() const override {return "ThreeDimensional";}

  // Copies
  NDMaterial *getCopy() override;
  NDMaterial *getCopy(const char *type) override;

  int getOrder() const override { return 6; }
  double getRho() override { return density_; }

  // State setting
  int setTrialStrain (const Vector &v) override;
  int setTrialStrainIncr (const Vector &dv) override;

  // Accessors
  const Matrix &getTangent() override;
  const Matrix &getInitialTangent() override;
  const Vector &getStress() override;
  const Vector &getStrain() override;
  bool threadSafe() const override {return true;}

  // History management
  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;

  // TaggedObject interface
  void Print(OPS_Stream &s, int flag) override;

private:
  // Parameters
  double E, nu, Fy, G;
  // Mass density
  double density_;
  // Isotropic hardening
  double Q_[2], b_[2], Hiso_;
  // Kinematic hardening
  std::vector<double> Ck_, gammak_;
  // Solver control
  double newton_tolerance;
  int    MaxIter_;

  // Tangent tensors
  OpenSees::MatrixND<6,6> Ce, C;

  // Trial/total strain & wrappers for OpenSees API
  OpenSees::VectorND<6>   eps;   // current total strain
  Matrix                  retInitialTangent, retTangent;
  Vector                  retStress, retStrain;

  // History
  struct State {
    OpenSees::VectorND<6> sig;     // Cauchy stress
    OpenSees::VectorND<6> eps;     // total strain at commit
    OpenSees::VectorND<6> eps_p;   // plastic strain tensor
    std::vector<OpenSees::VectorND<9>> sig_b; // back-stress components (deviatoric)
    double mises_strain;           // mises equivalent plastic strain
    double Estr;                   // strain energy
    double Ehst;                   // hysteretic energy
  } past, pres;
  std::vector<double> wk_;

private:
  // Core update
  int updateState(const VectorND<6> &eps);
  // bool isLinearHardening() const {return (Ck_.size() == 1 && gammak_[0] == 0.0);}


  // 6x9 Voight mapping used like: (P^vec6) -> 9x1   and   (P*vec9) -> 6x1
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

  //
  // constants for different conventions in scaling the yield function, 
  // flow direction, and consistency parameter
  //
  static constexpr double SQRT23 = 0.81649658092772603; // sqrt(2/3)
  static constexpr double mises_scale = SQRT23;
  static constexpr double norm_unit = 1.0; // SQRT23; // sqrt(3/2) for J2 norm
  static constexpr double flow_rate = 1.0; // SQRT23;
  // derived
  static constexpr double phi_n = norm_unit;
  // d(mises_strain)/d(lamda) = flow_rate*norm_unit*mises_scale
  static constexpr double mises_rate = flow_rate*phi_n*mises_scale;
  static constexpr double phi_m = flow_rate*phi_n*phi_n;
  static constexpr double phi_nnmm = (phi_n*phi_n)/(phi_m*phi_m);


  static inline double metric(const VectorND<9>& v) noexcept {
    return v.norm()*phi_n;
  }

  inline void newton_update(const VectorND<9>& s_tr, 
                            double lambda,
                            // outputs:
                            VectorND<9>& m,
                            double& g, 
                            double& Dg,
                            double& phi_a) const noexcept;

  inline void scale_tolerance(double sig_nrm, double& gtol, double& ftol) const noexcept;


  //
  // Hardening Models (TODO: these need refactoring and organizing)
  //


  // Nonlinear Chaboche kinematic hardening
  struct Kinematic {

    // method for integrating backstress; currently set at compile time
    enum class BackStressIntegration {
      Exponential, // Hartloper version
      BackwardEuler
    };
    static constexpr BackStressIntegration bs_integration =  BackStressIntegration::Exponential; // BackStressIntegration::BackwardEuler;

    static inline VectorND<9> A(const NonlinearJ2 &m,
                                const NonlinearJ2::State &past,
                                double lamda) noexcept {
      // Return the *partial* back-stress x_phi = sum phi_k x_k
      // to form the partial stress Z = s^tr - x_phi
      VectorND<9> x_phi{};
      const size_t nc = past.sig_b.size();
      for (size_t i=0;i<nc;i++) {
        // const double phi = 1.0/(1.0 + m.gammak_[i]*lamda);
        
        double phi; //, dphi;
        switch (bs_integration) {
          case BackStressIntegration::BackwardEuler:
            phi = 1.0/(1.0 + m.gammak_[i]*lamda);
            // dphi = -m.gammak_[i]*phi*phi;
            break;
          case BackStressIntegration::Exponential:
            phi = std::exp(-m.gammak_[i]*lamda);
            // dphi = -m.gammak_[i]*phi;
            break;
        }
        x_phi.addVector(1.0, past.sig_b[i], phi);
      }
      return x_phi;
    }

    static double H(const NonlinearJ2& m, 
                    const NonlinearJ2::State &past,
                    double mises_strain) noexcept {
      double H = 0.0;
      const size_t nc = m.Ck_.size();
      for (size_t i=0;i<nc;i++) {
        if (m.Ck_[i] == 0.0 || m.gammak_[i] == 0.0)
          continue;
        
        double phi;
        switch (bs_integration) {
          case BackStressIntegration::BackwardEuler:
            phi = 1.0/(1.0 + m.gammak_[i]*mises_strain);
            break;
          case BackStressIntegration::Exponential:
            phi = std::exp(-m.gammak_[i]*mises_strain);
            break;
        }
        // const double phi = 1.0/(1.0 + m.gammak_[i]*mises_strain);
        //
        // NOTE: For backward Euler, this is equivalent to:
        // H += Ck[i]*phi*lambda;
        H += (1.0 - phi)*m.Ck_[i]/m.gammak_[i];
      }
      return H;
    }

    static inline double dH(const NonlinearJ2& m,
                            const NonlinearJ2::State &s,
                            double d_mises_strain) noexcept {
      double dH = 0.0;
      const size_t nc = m.Ck_.size();
      for (size_t i=0;i<nc;i++) {
        if (m.Ck_[i] == 0.0 || m.gammak_[i] == 0.0)
          continue;

        double phi, dphi;
        switch (bs_integration) {
          case BackStressIntegration::BackwardEuler:
            phi = 1.0/(1.0 + m.gammak_[i]*d_mises_strain);
            dphi = -m.gammak_[i]*phi*phi;
            break;
          case BackStressIntegration::Exponential:
            phi = std::exp(-m.gammak_[i]*d_mises_strain);
            dphi = -m.gammak_[i]*phi;
            break;
        }
        // const double phi = 1.0/(1.0 + m.gammak_[i]*d_mises_strain);
        // const double dphi = -m.gammak_[i]*phi*phi;
        //
        dH += -(m.Ck_[i]/m.gammak_[i])*dphi;
      }
      return dH;
    }

    static inline double dX(const NonlinearJ2 &m,
                            const NonlinearJ2::State &past,
                            double lamda,
                            const VectorND<9>& n) noexcept {
      // Return |x|' = n . x'
      double dX = 0.0;
      // const double gamma = past.mises_strain + SQRT23*lamda;
      const size_t nc = m.gammak_.size();
      for (size_t i=0;i<nc;i++) {
        double phi, dphi;
        switch (bs_integration) {
          case BackStressIntegration::BackwardEuler:
            phi = 1.0/(1.0 + m.gammak_[i]*lamda);
            dphi = -m.gammak_[i]*phi*phi;
            break;
          case BackStressIntegration::Exponential:
            phi = std::exp(-m.gammak_[i]*lamda);
            dphi = -m.gammak_[i]*phi;
            break;
        }
        //
        dX -= dphi * n.dot(past.sig_b[i]);
      }
      return dX;
    }

    static inline void addTangent(const NonlinearJ2 &m,
                                  const NonlinearJ2::State &past,
                                  double lamda,
                                  const VectorND<9>& n,
                                  double theta,
                                  // outputs:
                                  MatrixND<6,6,double> &C) noexcept {
      const size_t nc = past.sig_b.size();
      // const VectorND<6> Pn = P * n;
      const VectorND<6> Pn = Voight::ReduceVector(n);
      for (size_t i=0; i<nc; i++) {
        double phi, dphi;
        switch (bs_integration) {
          case BackStressIntegration::BackwardEuler:
            phi = 1.0/(1.0 + m.gammak_[i]*lamda);
            dphi = -m.gammak_[i]*phi*phi;
            break;
          case BackStressIntegration::Exponential:
            phi = std::exp(-m.gammak_[i]*lamda);
            dphi = -m.gammak_[i]*phi;
            break;
        }
        //
        // VectorND<6> Px = P * past.sig_b[i];
        const VectorND<6> Px = Voight::ReduceVector(past.sig_b[i]);
        
        static constexpr double cc = -mises_rate;
        C.addTensorProduct(Px, Pn,  dphi*theta*cc);
        C.addTensorProduct(Pn, Pn, -dphi*theta*n.dot(past.sig_b[i])*phi_nnmm*cc);
      }
    }
  
    static inline void update(const NonlinearJ2 &m,
                              const NonlinearJ2::State &past,
                              double lamda,
                              const VectorND<9>& n,
                              NonlinearJ2::State &pres
                            ) noexcept {
      const size_t nc = past.sig_b.size();
      for (size_t i=0;i<nc;i++) {
        if (m.Ck_[i] == 0.0 || m.gammak_[i] == 0.0)
          continue;

        double phi;
        switch (bs_integration) {
          case BackStressIntegration::BackwardEuler:
            phi = 1.0/(1.0 + m.gammak_[i]*lamda);
            break;
          case BackStressIntegration::Exponential:
            phi = std::exp(-m.gammak_[i]*lamda);
            break;
        }
        pres.sig_b[i] = past.sig_b[i]*phi;
        const double chi = (1.0 - phi)*(m.Ck_[i]/m.gammak_[i]);
        pres.sig_b[i].addVector(1.0, n, chi*flow_rate/(mises_scale*norm_unit));
        // for backward Euler, this is equivalent to:
        // pres.sig_b[i].addVector(1.0, n, m.Ck_[i]*lamda*phi);
      }
    }
  };

  struct Isotropic {
    static inline double H(const NonlinearJ2& m,
                    const NonlinearJ2::State &past,
                    double dep
    ) noexcept {
      double mises_strain = past.mises_strain + dep;
      double H = m.Hiso_;
      for (int i=0;i<2;i++) {
        if (m.Q_[i] != 0.0 && m.b_[i] != 0.0)
          H += m.Q_[i]*m.b_[i]*std::exp(-m.b_[i]*mises_strain);
      }
      return H;
    }

    static inline double Y(const NonlinearJ2& m,
                    const NonlinearJ2::State &past,
                    double dep) noexcept {
      double mises_strain = past.mises_strain + dep;
      double Fy = m.Fy + m.Hiso_*mises_strain;
      for (int i=0;i<2;i++) {
        if (m.Q_[i] != 0.0 && m.b_[i] != 0.0)
          Fy += m.Q_[i]*(1.0 - std::exp(-m.b_[i] * mises_strain));
      }
      return Fy;
    }

    static inline void update(const NonlinearJ2& m,
                              const NonlinearJ2::State &past,
                              double d_mises_strain,
                              // outputs:
                              NonlinearJ2::State &pres) noexcept {
      pres.mises_strain = past.mises_strain + d_mises_strain;
    }
  };
};

} // namespace OpenSees

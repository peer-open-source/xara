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
#pragma once
#include <NDMaterial.h>
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

  // Copies
  NDMaterial *getCopy() override;
  NDMaterial *getCopy(const char *type) override;

  // Comm/Print
  int sendSelf(int ctag, Channel &) override;
  int recvSelf(int ctag, Channel &, FEM_ObjectBroker &) override;
  void Print(OPS_Stream &s, int flag) override;

private:
  // Parameters
  double E, nu, fy, G;
  // Isotropic hardening
  double a_, DInf_, b_, QInf_, Hiso_;
  // Kinematic hardening
  std::vector<double> Ck_, gammak_;
  // Solver controls
  double newton_tolerance;
  int    MaxIter_;
  // Mass density
  double density_;

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
    OpenSees::VectorND<6> eps_p;   // plastic strain
    std::vector<OpenSees::VectorND<9>> sig_b; // back-stress components (deviatoric, 9x1 each)
    double gamma;                  // iso hardening scalar
    double Estr;                   // strain energy
    double Ehst;                   // hysteretic energy
  } past, pres;
  std::vector<double> wk_;

private:
  // Core update
  int updateState();
  bool isLinearHardening() const {return (Ck_.size() == 1 && gammak_[0] == 0.0);}


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


  // constants for different conventions in scaling the yield check
  static constexpr double SQRT23 = 0.81649658092772603; // sqrt(2/3)
  static constexpr double mises_scale = SQRT23;
  static constexpr double norm_unit = 1.0; // /SQRT23; // sqrt(3/2) for J2 norm
  static constexpr double flow_rate = 1.0; // SQRT23;
  // derived
  static constexpr double phi_n = norm_unit;
  // d(mises_strain)/d(lamda) = flow_rate*norm_unit*mises_scale
  static constexpr double mises_rate = flow_rate*phi_n*mises_scale;
  static constexpr double phi_m = flow_rate*phi_n*phi_n;
  static constexpr double phi_nnmm = (phi_n*phi_n)/(phi_m*phi_m);

  // method for integrating backstress; currently set at compile time
  enum class BackStressIntegration {
    Exponential, // Hartloper version
    BackwardEuler
  } bs_integration =  BackStressIntegration::Exponential; // BackStressIntegration::BackwardEuler;


  static inline double metric(const VectorND<9>& v) noexcept {
    return v.norm()*phi_n;
  }

  inline void newton_update(const VectorND<9>& s_tr, 
                            double lambda,
                            // outputs:
                            VectorND<9>& m,
                            double& g, 
                            double& Dg,
                            double& phi_a) const noexcept {
    double mises_strain = (flow_rate*lambda)*phi_n*mises_scale;

    // Partial stress, Z
    m = s_tr;
    m -= Kinematic::A(*this, past, mises_strain);

    phi_a  = metric(m); // phi(a)

    // Newton residual function
    g = phi_a - 2.0*G*mises_strain*(phi_n/mises_scale)
      - Kinematic::H(*this, past, mises_strain)*(phi_n/mises_scale)
      - Isotropic::Y(*this, past, mises_strain)*phi_n*mises_scale;

    if (phi_a < 1e-16) {
      m.zero();
      Dg = 0.0;
      return;
    }

    m    *= phi_m/phi_a; // m = a * phi(m)/phi(a)

    // Derivative of newton residual
    // Note: d(norm(Z))/d(lamda) == n . X' == dX(n)
    // const double dphi_a = 
    Dg = Kinematic::dX(*this, past, mises_strain, m)*mises_rate/flow_rate
       - Kinematic::dH(*this, past, mises_strain)*mises_rate*(phi_n/mises_scale)
       - 2.0*G*mises_rate*(phi_n/mises_scale)
       - Isotropic::H(*this, past, mises_strain)*mises_scale*phi_n*mises_rate;
  }

  inline void scale_tolerance(double sig_nrm, double& gtol, double& ftol) const noexcept {
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

  // Nonlinear Chaboche kinematic hardening
  struct Kinematic {
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
        switch (m.bs_integration) {
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
        switch (m.bs_integration) {
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
                            double lamda) noexcept {
      double dH = 0.0;
      const size_t nc = m.Ck_.size();
      for (size_t i=0;i<nc;i++) {
        if (m.Ck_[i] == 0.0 || m.gammak_[i] == 0.0)
          continue;
        
        double phi, dphi;
        switch (m.bs_integration) {
          case BackStressIntegration::BackwardEuler:
            phi = 1.0/(1.0 + m.gammak_[i]*lamda);
            dphi = -m.gammak_[i]*phi*phi;
            break;
          case BackStressIntegration::Exponential:
            phi = std::exp(-m.gammak_[i]*lamda);
            dphi = -m.gammak_[i]*phi;
            break;
        }
        // const double phi = 1.0/(1.0 + m.gammak_[i]*lamda);
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
      // const double gamma = past.gamma + SQRT23*lamda;
      const size_t nc = m.gammak_.size();
      for (size_t i=0;i<nc;i++) {
        double phi, dphi;
        switch (m.bs_integration) {
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
                                  MatrixND<6,6,double> &C) noexcept {
      const size_t nc = past.sig_b.size();
      const VectorND<6> Pn = P * n;
      for (size_t i=0; i<nc; i++) {
        double phi, dphi;
        switch (m.bs_integration) {
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
        VectorND<6> Px = P * past.sig_b[i];
        
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
        switch (m.bs_integration) {
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
                    const NonlinearJ2::State &s,
                    double lamda) noexcept {
      double gamma = s.gamma + lamda;
      double H = m.Hiso_;
      if (m.DInf_ != 0.0 && m.a_ != 0.0)
        H += m.DInf_*m.a_*std::exp(-m.a_*gamma);
      if (m.QInf_ != 0.0 && m.b_ != 0.0)
        H += m.QInf_*m.b_*std::exp(-m.b_*gamma);
      return H;
    }

    static inline double Y(const NonlinearJ2& m,
                    const NonlinearJ2::State &s,
                    double lamda) noexcept {
      double gamma = s.gamma + lamda;
      double Fy = m.fy + m.Hiso_*gamma;
      if (m.QInf_ != 0.0 && m.b_ != 0.0)
        Fy += m.QInf_*(1.0 - std::exp(-m.b_ * gamma));
      if (m.DInf_ != 0.0 && m.a_ != 0.0)
        Fy += m.DInf_*(1.0 - std::exp(-m.a_ * gamma));
      return Fy;
    }

    static inline void update(const NonlinearJ2& m,
                              const NonlinearJ2::State &past,
                              double lamda,
                              NonlinearJ2::State &pres) noexcept {
      pres.gamma = past.gamma + lamda;
    }
  };
};

} // namespace OpenSees

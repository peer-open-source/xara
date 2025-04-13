//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#pragma once
#include <State.h>
#include <Field.h>
#include <material/section/SectionForceDeformation.h>

struct FrameSectionConstants {
  // n-n
  double A;
  double Ay, 
         Az;
  // m-m
  double Iy, 
         Iz, 
         Iyz;
  // w-w
  double Cw, 
         Ca;
  // n-m
  double Qy, 
         Qz;
  // n-w
  double Rw, 
         Ry, 
         Rz;
  // m-w
  double Sa, 
         Sy, 
         Sz;
};

enum FrameStress : int {
  End         =     0,
  N           =     2, //
  Vy          =     3, // 0b00000010
  Vz          =     5, // 0b00000100
  T           =     6, // 0b00000000
  My          =     4, // 0b00000000
  Mz          =     1, // 0b00000000
  R           =     7, // (Obselete, also bishear)
  Q           =     8, // (Obselete, also bimoment)
  Bimoment    =     9, // 
  Wagner      =    10, // (Obselete, this is redundant)
  Bishear     =    11,
  By, Bz,
  Qy, Qz,
  Max,
};

struct FrameLayout {
  int n[3], m[3], w[3], v[3];
};

typedef int FrameStressLayout[FrameStress::Max];

static inline constexpr FrameLayout
WarpIndex(const int n, const FrameStressLayout& layout) {
    
    FrameLayout L {{-1, -1, -1}, {-1, -1, -1},
                   {-1, -1, -1}, {-1, -1, -1}};
    // Save layout locations
    for (int i=0; i<n; i++) {
      switch (layout[i]) {
        case FrameStress::N:        L.n[0] = i;  break;
        case FrameStress::Vy:       L.n[1] = i;  break;
        case FrameStress::Vz:       L.n[2] = i;  break;
        case FrameStress::T:        L.m[0] = i;  break;
        case FrameStress::My:       L.m[1] = i;  break;
        case FrameStress::Mz:       L.m[2] = i;  break;
        case FrameStress::Bimoment: L.w[0] = i;  break;
        case FrameStress::By:       L.w[1] = i;  break;
        case FrameStress::Bz:       L.w[2] = i;  break;
        case FrameStress::Bishear:  L.v[0] = i;  break;
        case FrameStress::Qy:       L.v[1] = i;  break;
        case FrameStress::Qz:       L.v[2] = i;  break;
        default:
        ;
      }
    }
  return L;
}

class FrameSection : public SectionForceDeformation {

public:
  FrameSection(int tag, int clstag, double mass=0, bool use_mass=false)
    : SectionForceDeformation(tag, clstag),
      density(mass), has_mass(use_mass)
  {}

  virtual FrameSection* getFrameCopy() =0;
  virtual FrameSection* getFrameCopy(const FrameStressLayout& layout) {
    return getFrameCopy();
  }

  virtual SectionForceDeformation* getCopy() {
    return this->getFrameCopy();
  }

  virtual int getIntegral(Field field, State state, double& value) const {
    if ((field == Field::Density) && has_mass) {
        value = density;
        return 0;
    }
    return -1;
  }

  virtual int setFiberValue(int tag, int field, double value) {
    return -1;
  }

  template <int n, const FrameStressLayout& scheme>
  int setTrialState(OpenSees::VectorND<n, double> e);


  template <int n, const FrameStressLayout& scheme>
  OpenSees::VectorND<n> 
  getResultant() {

    OpenSees::VectorND<n> sout;

    const ID& layout = this->getType();

    int m = this->getOrder();

    const Vector& s = this->getStressResultant();
    for (int i=0; i<n; i++) {
      sout[i] = 0.0;
      for (int j=0; j<m; j++)
        if (layout(j) == scheme[i])
          sout[i] = s(j);
    }

    return sout;
  }

  template <int n, const FrameStressLayout& scheme>
  OpenSees::MatrixND<n,n> getTangent(State state) {

    OpenSees::MatrixND<n,n> kout;

    const ID& layout = this->getType();

    int m = this->getOrder();

    const Matrix& ks = (state == State::Init)
                     ? this->getInitialTangent()
                     : this->getSectionTangent();


    constexpr FrameLayout e = WarpIndex(n, scheme);

    int sect_bishear[3] = {-1,-1,-1};

    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        kout(i,j) = 0.0;

        for (int k=0; k<m; k++) {
          switch (layout(k)) {
            case FrameStress::Bishear:
              sect_bishear[0]  = k;
              break;
            case FrameStress::Qy:
              sect_bishear[1]  = k;
              break;
            case FrameStress::Qz:
              sect_bishear[2]  = k;
              break;
          }

          if (layout(k) == scheme[i]) {
            for (int l=0; l<m; l++)
              if (layout(l) == scheme[j])
                kout(i,j) = ks(k,l);
          }
        }
      }
    }

    // If element has a twisting DOF and no Bishear
    // DOF, then twist == alpha, where alpha is the
    // bishear DOF.
    if constexpr (e.m[0] != -1 && e.v[0] == -1)
      if (sect_bishear[0] != -1) {
        kout(e.m[0],e.m[0]) += ks(sect_bishear[0],sect_bishear[0]);
        for (int i=0; i<n; i++)
          for (int j=0; j<m; j++)
            if (layout(j) == scheme[i]) {
              // if (std::fabs(ks(j,sect_bishear[0])) > 1e-12)
                kout(i,e.m[0]) += ks(j,sect_bishear[0]);
              // if (std::fabs(ks(sect_bishear[0],j)) > 1e-12)
                kout(e.m[0],i) += ks(sect_bishear[0],j);
            }
      }

    // If element has a shear (Vy) DOF and no Qy
    if constexpr (e.n[1] != -1 && e.v[1] == -1)
      if (sect_bishear[1] != -1) {
        kout(e.n[1],e.n[1]) += ks(sect_bishear[1],sect_bishear[1]);
        for (int i=0; i<n; i++)
          for (int j=0; j<m; j++)
            if (layout(j) == scheme[i]) {
              kout(i,e.n[1]) += ks(j,sect_bishear[1]);
              kout(e.n[1],i) += ks(sect_bishear[1],j);
            }
      }

    // If element has a shear (Vz) DOF and no Qz
    if constexpr (e.n[2] != -1 && e.v[2] == -1)
      if (sect_bishear[2] != -1) {
        kout(e.n[2],e.n[2]) += ks(sect_bishear[2],sect_bishear[2]);
        for (int i=0; i<n; i++)
          for (int j=0; j<m; j++)
            if (layout(j) == scheme[i]) {
              kout(i,e.n[2]) += ks(j,sect_bishear[2]);
              kout(e.n[2],i) += ks(sect_bishear[2],j);
            }
      }
    return kout;
  }

  template <int n, const FrameStressLayout& scheme>
  OpenSees::MatrixND<n,n, double> getFlexibility(State state=State::Pres);

private:
  double density;
  bool has_mass;
};

template <int n, const FrameStressLayout& scheme>
int 
FrameSection::setTrialState(OpenSees::VectorND<n, double> e) {
  double strain_data[FrameStress::Max]{};

  const int m = this->getOrder();
  Vector trial(strain_data, m);
  trial.Zero();

  const ID& layout = this->getType();


  constexpr FrameLayout l = WarpIndex(n, scheme);

  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++)
      if (layout(j) == scheme[i])
        trial[j] = e[i];
  }

  // Case 2 and 3
  // If element has a twisting DOF and no Bishear
  // DOF, then twist == alpha, where alpha is the
  // bishear DOF.
  // Note that elem_twist and elem_bishear are computable
  // at compile time, so this branch can theoretically be 
  // optimized out by the compiler, however this might be 
  // optimistic
  //
  if (l.v[0] == -1) {
    for (int j=0; j<m; j++)
      switch (layout(j)) {
        case FrameStress::Bishear:
          // Set alpha = tau
          if (l.m[0] != -1)
            trial[j] = e[l.m[0]];
          break;
        case FrameStress::Qy:
          // Set alpha_y = gamma_y
          if (l.n[1] != -1)
            trial[j] = e[l.n[1]];
          break;
        case FrameStress::Qz:
          // Set alpha_z = gamma_z
          if (l.n[2] != -1)
            trial[j] = e[l.n[2]];
          break;
        default:
          ;
      }
  }
  return this->setTrialSectionDeformation(trial);
}

template <int n, const FrameStressLayout& scheme>
OpenSees::MatrixND<n,n, double> 
FrameSection::getFlexibility(State state)
{
  OpenSees::MatrixND<n,n,double> K = getTangent<n,scheme>(state);
  OpenSees::MatrixND<n,n,double> F;
  K.invert(F);
  return F;

  OpenSees::MatrixND<n,n,double> Fout;


  const ID& layout = this->getType();

  int m = this->getOrder();

  const Matrix& Fs = (state == State::Init)
                    ? this->getInitialFlexibility()
                    : this->getSectionFlexibility();

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {  
      Fout(i,j) = 0.0;
      for (int k=0; k<m; k++) {
        if (layout(k) == scheme[i]) {
          for (int l=0; l<m; l++)
            if (layout(l) == scheme[j]) {
              Fout(i,j) = Fs(k,l);
            }
        }
      }
    }
  }

  return Fout;
}

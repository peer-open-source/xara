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
#pragma once
#include "FrameSection.h"

namespace {


struct FrameLayout {
  int n[3], m[3], w[3], v[3];
};

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


static inline constexpr
OpenSees::MatrixND<12,12>
SectionCondensation(const FrameLayout& e)
{
  OpenSees::MatrixND<12,12> L{};
  L.addDiagonal(1.0);

  if ((e.m[0] != -1) && (e.v[0] == -1)) {
    L(9, 3) = 1.0;
    L(9, 9) = 0.0;
  }
  if (e.n[1] != -1) {
    if (e.v[1] == -1) {
      L(10, 1) = 1.0;
      L(10,10) = 0.0;
    }
  }
  if (e.n[2] != -1) {
    if (e.v[2] == -1) {
      L(11, 2) = 1.0;
      L(11,11) = 0.0;
    }
  }
  return L;
}
} // namespace

template <int n, const FrameStressLayout& scheme>
OpenSees::VectorND<n>
FrameSection::getResultant() {

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

    constexpr FrameLayout e = WarpIndex(n, scheme);
    if constexpr ((e.n[1] != -1) && (e.v[1] == -1))
      if (m == 12)
        sout[e.n[1]] += s[10];

    if constexpr ((e.n[2] != -1) && (e.v[2] == -1))
      if (m == 12)
        sout[e.n[2]] += s[11];

    return sout;
  }

template <int n, const FrameStressLayout& scheme>
OpenSees::MatrixND<n,n> 
FrameSection::getTangent(State state) {

  OpenSees::MatrixND<n,n> kout;

  int m = 12;

  OpenSees::MatrixND<12,12> Ks = this->getFullTangent(state);

  constexpr FrameLayout e = WarpIndex(n, scheme);

  constexpr OpenSees::MatrixND<12,12> L = SectionCondensation(e);

  OpenSees::MatrixND<12,12> LKL = L^(Ks*L);
  for (int i=0; i<12; i++)
    for (int j=0; j<12; j++)
      Ks(i,j) = LKL(i,j);

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      kout(i,j) = 0.0;
      for (int k=0; k<m; k++) {
        if (FullLayout[k] == scheme[i]) {
          for (int l=0; l<m; l++)
            if ((FullLayout[l] == scheme[j]))
              kout(i,j) = Ks(k,l);
        }
      }
    }
  }

  return kout;
}


template <int n, const FrameStressLayout& scheme>
int 
FrameSection::setTrialState(const OpenSees::VectorND<n>& e) {
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
  if constexpr (l.v[0] == -1) 
  {
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
  // TODO: clean this up, validate
  OpenSees::MatrixND<n,n,double> F;
  Cholesky<n>(K).invert(F);
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

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

using OpenSees::MatrixND;
using OpenSees::VectorND;

struct FrameLayout {
  int n[3], m[3], w[3], v[3];
};

static inline constexpr FrameLayout
WarpIndex(const int n, const FrameStressLayout& layout) noexcept {
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
MatrixND<12,12>
ConstraintMatrix(const FrameLayout& e) noexcept
{
  MatrixND<12,12> L{};

  L.addDiagonal(1.0);

  if ((e.m[0] != -1) && (e.v[0] == -1)) {
    L(9, 3) = 1.0;
    L(9, 9) = 0.0;
  }

  if ((e.n[1] != -1) && (e.v[1] == -1)) {
    if (e.v[1] == -1) {
      L(10, 1) = 1.0;
      L(10,10) = 0.0;
    }
  }
  if ((e.n[2] != -1) && (e.v[2] == -1)) {
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
FrameSection::getResultant() noexcept {
  using namespace OpenSees;

  VectorND<n> sout{};

  const VectorND<12> s_full = this->getFullStress();

  for (int i=0; i<n; i++)
    for (int j=0; j<12; j++)
      if (FullLayout[j] == scheme[i])
        sout[i] = s_full(j);

  return sout;
}


template <int n, const FrameStressLayout& scheme>
OpenSees::MatrixND<n,n> 
FrameSection::getTangent(State state) noexcept
{

  constexpr static int m = 12;

  OpenSees::MatrixND<n,n> kout{};

  OpenSees::MatrixND<12,12> Ks = this->getFullTangent(state);


  if (!getenv("XARA_FIBER_THREADS")) {
    constexpr FrameLayout iw = WarpIndex(n, scheme);
    constexpr
    OpenSees::MatrixND<12,12> L  = ConstraintMatrix(iw);
    OpenSees::MatrixND<12,12> LKL = Ks*L;
    for (int i=0; i<m; i++)
      for (int j=0; j<m; j++)
        Ks(i,j) = LKL(i,j);
  }


  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      kout(i,j) = 0.0;
      for (int k=0; k<m; k++) {
        if (FullLayout[k] == scheme[i]) {
          for (int l=0; l<m; l++)
            if (FullLayout[l] == scheme[j])
              kout(i,j) = Ks(k,l);
        }
      }
    }
  }
  return kout;
}


template <int n, const FrameStressLayout& scheme>
int 
FrameSection::setTrialState(const OpenSees::VectorND<n>& e) noexcept
{
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
  if ((l.v[0] == -1) && !getenv("XARA_FIBER_THREADS")) 
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
OpenSees::MatrixND<n,n> 
FrameSection::getFlexibility(State state) noexcept
{
  static OpenSees::MatrixND<n,n> K;
  K = getTangent<n,scheme>(state);
  OpenSees::MatrixND<n,n> F{};
  K.invert(F);
  // Cholesky<n>(K).invert(F);
  return F;
}


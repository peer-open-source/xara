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
#include <State.h>
#include <Field.h>
#include <Cholesky.tpp>
#include <material/section/SectionForceDeformation.h>
#include <cmath>

// TODO: Maybe make this public under ElasticFrameSection
struct FrameSectionConstants {
  // n-n
  double A;
  double Ay, Az;
  // m-m
  double Iy, Iz, Iyz;
  // w-w
  double Cw, Ca;
  // n-m
  double Qy, Qz;
  // n-w
  double Rw, Ry, Rz;
  // m-w
  double Sa, Sy, Sz;
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


typedef int FrameStressLayout[FrameStress::Max];



class FrameSection : public SectionForceDeformation {

public:
  FrameSection(int tag, int clstag, double mass=0, bool use_mass=false)
    : SectionForceDeformation(tag, clstag),
      density(mass), has_mass(use_mass)
  {}

  struct Tangent {
    OpenSees::MatrixND<3,3> nn,     nw, nv,
                            mn, mm, mw, mv, 
                                    ww,
                                        vv;
    void zero() {
      nn.zero();            nw.zero(); nv.zero();
      mn.zero(); mm.zero(); mw.zero(); mv.zero();
      ww.zero();
      vv.zero();
    }
  };

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

  //
  // New response API
  //
  template <int n, const FrameStressLayout& scheme>
  int setTrialState(const OpenSees::VectorND<n>& e);

  virtual OpenSees::MatrixND<12,12>
  getFullTangent(State state) {
    const Matrix& ks = (state == State::Init)
                      ? this->getInitialTangent()
                      : this->getSectionTangent();

    int m = this->getOrder();
    const ID& layout = this->getType();

    OpenSees::MatrixND<12,12> Ks;
    for (int i=0; i<12; i++) {
      for (int j=0; j<12; j++) {
        Ks(i,j) = 0.0;
        for (int k=0; k<m; k++)
          if (layout(k) == FullLayout[i]) {
            for (int l=0; l<m; l++)
              if (layout(l) == FullLayout[j])
                Ks(i,j) = ks(k,l);
          }
      }
    }
    return Ks;
  }

  template <int n, const FrameStressLayout& scheme>
  OpenSees::VectorND<n> getResultant();

  template <int n, const FrameStressLayout& scheme>
  OpenSees::MatrixND<n,n> getTangent(State state); 

  template <int n, const FrameStressLayout& scheme>
  OpenSees::MatrixND<n,n, double> getFlexibility(State state=State::Pres);

  template <int n, const FrameStressLayout& scheme>
  OpenSees::VectorND<n> 
  getResultantGradient(int grad, bool conditional) {

    OpenSees::VectorND<n> sout;

    const ID& layout = this->getType();

    int m = this->getOrder();

    const Vector& s = this->getStressResultantSensitivity(grad, conditional);
    for (int i=0; i<n; i++) {
      sout[i] = 0.0;
      for (int j=0; j<m; j++)
        if (layout(j) == scheme[i])
          sout[i] = s(j);
    }

    return sout;
  }


private:
  double density;
  bool   has_mass;

  static constexpr FrameStressLayout FullLayout = {
    FrameStress::N, FrameStress::Vy, FrameStress::Vz,
    FrameStress::T, FrameStress::My, FrameStress::Mz,
    FrameStress::Bimoment, FrameStress::By, FrameStress::Bz,
    FrameStress::Bishear,  FrameStress::Qy, FrameStress::Qz
  };
};

//
// Inlines
//
#include "FrameSection.tpp"

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
#include <MatrixND.h>
#include <VectorND.h>
#include <Cholesky.tpp>
#include <material/section/SectionForceDeformation.h>
#include <cmath>
#include "FrameSectionConstants.h"

using namespace OpenSees;

enum FrameStress : int {
  End         =     0,
  N           =     2,
  Vy          =     3,
  Vz          =     5,
  T           =     6,
  My          =     4,
  Mz          =     1,
  R           =     7, // (Obselete, also bishear)
  Q           =     8, // (Obselete, also bimoment)
  Bimoment    =     9, // 
  Wagner      =    10, // (Obselete, this is redundant)
  Bishear     =    11,
  By, Bz,
  Qy, Qz,
  Max=16,
};


typedef int FrameStressLayout[FrameStress::Max];



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

  virtual int getShape(Frame::Prism& shape) {
    
    // 1) Get exact reference properties; not all sections provide these

    double value;
    if (this->getIntegral(Field::Unit,   State::Init, value) == 0)
      shape.A = value;
    else
      shape.A = 1.0;

    if (this->getIntegral(Field::UnitZZ, State::Init, value) == 0)
      shape.Iy = value;

    if (this->getIntegral(Field::UnitYY, State::Init, value) == 0)
      shape.Iz = value;

    // 2) Get Young and Shear Modulus and determine if shear is supported
    // by the section. The shear areas we pull here may still be 
    // uncondensed.
    const ID& layout = this->getType();
    const Matrix& Ks = this->getInitialTangent();
    for (int i=0; i<layout.Size(); i++) {
      if (layout(i) == FrameStress::N && shape.A) {
        shape.E = Ks(i,i)/(shape.A);
      }
      else if (layout(i) == FrameStress::Vy && shape.A) {
        shape.G = Ks(i,i)/(shape.A);
        shape.Ay = shape.A;
      }
      else if (layout(i) == FrameStress::Vz && shape.A) {
        shape.G = Ks(i,i)/(shape.A);
        shape.Az = shape.A;
      }
    }
    //
    for (int i=0; i<layout.Size(); i++) {
      if (layout(i) == FrameStress::My && !shape.Iy && shape.E) {
        shape.Iy = Ks(i,i)/(*shape.E);
      }
      else if (layout(i) == FrameStress::Mz && !shape.Iz && shape.E) {
        shape.Iz = Ks(i,i)/(*shape.E);
      }
    }
    // In a 3D shear-free section, G wouldnt have been found yet
    for (int i=0; i<layout.Size(); i++) {
      if (layout(i) == FrameStress::T && shape.Iy && shape.Iz && !shape.G) {
        shape.G = Ks(i,i)/(*shape.Iy + *shape.Iz);
      }
    }

    // 3) Condense Warping
    static constexpr FrameStressLayout scheme = {
        FrameStress::N,
        FrameStress::Vy,
        FrameStress::Vz,
        FrameStress::T,
        FrameStress::My,
        FrameStress::Mz,
    };

    MatrixND<6,6> Kc = this->getTangent<6,scheme>(State::Init);

    if (shape.G) {
      shape.J  = Kc(3,3)/(*shape.G);
      if (shape.Ay)
        shape.Ay = Kc(1,1)/(*shape.G);
      if (shape.Az)
        shape.Az = Kc(2,2)/(*shape.G);
    }
    return 0;
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
  int setTrialState(const OpenSees::VectorND<n>& e) noexcept;

  virtual OpenSees::VectorND<12>
  getFullStress() {
    const Vector& s = this->getStressResultant();

    const int m = this->getOrder();
    const ID& layout = this->getType();

    OpenSees::VectorND<12> S{};
    for (int i=0; i<m; i++) {
      const int k = layout(i);
      for (int j=0; j<12; j++)
        if (k == FullLayout[j])
          S[j] = s(i);
    }
    return S;
  }

  virtual OpenSees::MatrixND<12,12>
  getFullTangent(State state) {
    const Matrix& ks = (state == State::Init)
                      ? this->getInitialTangent()
                      : this->getSectionTangent();

    int m = this->getOrder();
    const ID& layout = this->getType();

    OpenSees::MatrixND<12,12> Ks{};
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
  OpenSees::VectorND<n> getResultant() noexcept;

  template <int n, const FrameStressLayout& scheme>
  OpenSees::MatrixND<n,n> getTangent(State state) noexcept; 

  template <int n, const FrameStressLayout& scheme>
  OpenSees::MatrixND<n,n, double> getFlexibility(State state=State::Pres) noexcept;

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

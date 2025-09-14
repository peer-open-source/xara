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
// Claudio M. Perez
//
#pragma once
#include <array>
#include <vector>
#include <Domain.h>
#include <string.h>
#include <Versor.h>
#include <VectorND.h>
#include <Matrix3D.h>
#include <GroupSO3.h>
#include <FiniteElement.h>
#include <ElementalLoad.h>
#include <LoadPattern.h>
#include <FrameSection.h>

class Element;

namespace OpenSees {

#define LOAD_TAG_FrameLoad 141414

class FrameLoad: public ElementalLoad 
{
private:
  constexpr static int classTag = LOAD_TAG_FrameLoad;

public:
  enum Basis {
    Embedding,
    Reference,
    Director,
  };
  enum Shape {
    Dirac,
    Heaviside,
    Lagrange,
  };
  FrameLoad(int basis, 
            int shape, 
            std::vector<Vector3D>& p,
            std::vector<Vector3D>& m,
            std::vector<Vector3D>& r,
            LoadPattern& pattern)
  : ElementalLoad(classTag),
    basis(basis),
    shape(shape),
    pattern(pattern),
    p(p),
    m(m),
    r(r)
  {
    switch (shape) {
      case Dirac:
          gauss = {{r[0][0], 1.0}};
          break;
      case Lagrange:
      case Heaviside:
        // 4-point Gauss-Legendre quadrature on [0, 1]
        gauss = {{0.069431844, 0.173927423},
                 {0.330009478, 0.326072577},
                 {0.669990522, 0.326072577},
                 {0.930568156, 0.173927423}};
        break;
    }
  }

  ~FrameLoad() {
    // NOTE: This is abusing the load factor
    // argument; the element has to recognize that for
    // this load type, zero load factor means delete from
    // your list of loads.
    for (auto e: elements)
      e->addLoad(this, 0.0);
  }

  virtual void
  setDomain(Domain *theDomain) final
  {
    this->Load::setDomain(theDomain);
  
    if (theDomain == nullptr) {
      for (auto e: elements)
        e->addLoad(this, 0.0);
      return;
    }
  }

  int 
  addElement(Element& element) 
  {
    auto name = element.getClassType();
    if (strstr(name, "Frame") == nullptr) {
      opserr << "WARNING FrameLoad::addElement() - cannot add load to element of type " << name << '\n';
      return -1;
    }
    elements.push_back(&element);
    element.addLoad(this, 1.0);
    return 0;
  }

  int
  recvSelf(int commitTag, Channel& , FEM_ObjectBroker&) final
  {
    return -1;
  }

  int
  sendSelf(int commitTag, Channel& ) final
  {
    return -1;
  }

  const std::vector<std::array<double,2>>& 
  quadrature () {
    return gauss;
  }

  void applyLoad(double loadFactor) final {
    for (auto e: elements)
      e->update();
  }
  
  virtual const Vector&
  getData(int& type, double loadFactor) override final {
    type = classTag;
    static Vector v(0);
    return v;
  }

  int getBasis() const {
    return basis;
  }

  bool conservative() const {
    return false;
  }


public:
  template <int i, int nn, int ndf>
  void addLoadAtPoint(VectorND<nn*ndf>& pe, 
                      double x, double w, double jxs,
                      const Matrix3D& R0,
                      const Matrix3D& R) const
  {   
    if (w == 0.0)
      return;

    for (unsigned q = 0; q < r.size(); q++) {
      Vector3D px,mx,rx;
      rx = r[q];
      rx[0] = 0.0;
      rx = R*(R0*rx);
      switch (basis) {
        case Embedding:
            px = p[q];
            mx = m[q] + rx.cross(px);
            break;
        case Reference:
            px = R0 * p[q];
            mx = R0 * m[q] + rx.cross(px);
            break;
        case Director:
            px = R * p[q];
            mx = R * m[q] + rx.cross(px);
            break;
      }

      double scale = -w*pattern.getLoadFactor();
      switch (shape) {
        case Dirac:
          scale /= jxs;
          if (std::fabs(x - r[q][0]) > 1.0e-6)
            scale *= 0.0;
          break;
        case Heaviside:
          if (x < r[q][0])
            scale *= 0.0;
          break;
        case Lagrange:
          for (unsigned s=0; s<r.size(); s++)
            if (s != q)
              scale *= (x - r[s][0]) / (r[q][0] - r[s][0]);
          break;
      }
      pe.template assemble<  i*ndf>(px, scale);
      pe.template assemble<3+i*ndf>(mx, scale);
    }
  }

  template <int i, int j, int nn, int n>
  void addTangAtPoint(MatrixND<nn*n,nn*n>& K, 
                      double x, double w, double jxs,
                      const Matrix3D& R0,
                      const Matrix3D& R) const
  {
    if (w == 0.0)
        return;

    for (unsigned q = 0; q < r.size(); q++) {
      Vector3D px, mx, rx;
      rx = r[q];
      rx[0] = 0.0;

      rx = R * rx;
      rx = R*(R0*rx);
      if (rx.norm() == 0.0 && basis == Embedding)
        continue;
      switch (basis) {
        case Embedding:
          px = p[q];
          mx = m[q] + rx.cross(px);
          break;
        case Reference:
          px = R0 * p[q];
          mx = R0 * m[q] + rx.cross(px);
          break;
        case Director:
          px = R * p[q];
          mx = R * m[q] + rx.cross(px);
          break;
      }

      double scale = -w*pattern.getLoadFactor();
      switch (shape) {
        case Dirac:
          scale /= jxs;
          if (std::fabs(x - r[q][0]) > 1.0e-6)
            scale *= 0.0;
          break;
        case Heaviside:
          if (x < r[q][0])
            scale *= 0.0;
          break;
        case Lagrange:
          for (unsigned s=0; s<r.size(); s++)
            if (s != q)
            scale *= (x - r[s][0]) / (r[q][0] - r[s][0]);
          break;
      }

      Matrix3D Px = Hat(px);
      if (basis == Director) {
        K.assemble(        Px,   i*n, 3+j*n, -scale);
        K.assemble(Hat(rx)*Px, 3+i*n, 3+j*n, -scale);
      }
      
      if (rx.norm() != 0.0)
        K.assemble(Px*Hat(rx), 3+i*n, 3+j*n, scale);
    }
  }

  template <int nsr, const FrameStressLayout& scheme>
  void addBasicSolution(VectorND<nsr>&    s, double x, double L, 
                        const Matrix3D& R0, const Matrix3D& R) const
  // add particular solution for basic equations
  {

    VectorND<3> nm{}, M{};
    {
      Vector3D rx = r[0];
      rx[0] = 0.0;
      rx = R*(R0*rx);
      switch (basis) {
        case Embedding:
          nm = R^p[0]; // + rx.cross(m[0]);
          break;
        case Reference:
          nm = R^(R0^p[0]);
          M  = (R^m[0]) + rx.cross(nm);
          break;
        case Director:
          nm = p[0];
          M  = m[0] + rx.cross(nm);
          break;
      }
    }


    double scale = pattern.getLoadFactor();

    switch (shape) {
      case Heaviside: {
        double wa = nm[0]*scale; // Axial
        double wy = nm[1]*scale; // Transverse
        double wz = nm[2]*scale; // Transverse

        for (int i = 0; i < nsr; i++) {
          switch (scheme[i]) {
          case FrameStress::N:  s[i] +=  wa * (L - x); break;
          case FrameStress::Vy: s[i] += M[2]*scale + wy * (x - 0.5 * L); break;
          case FrameStress::Vz: s[i] +=-M[1]*scale - wz * (x - 0.5 * L); break;

          case FrameStress::My: s[i] += - wz * 0.5 * x * (x - L); break;
          case FrameStress::Mz: s[i] += + wy * 0.5 * x * (x - L); break;
          default:
            break;
          }
        }
        break;
      }
      case Dirac: {
        double N      = nm[0]*scale;
        double Py     = nm[1]*scale;
        double Pz     = nm[2]*scale;
        double aOverL = r[0][0];

        if (aOverL < 0.0 || aOverL > 1.0)
          break;

        double a = aOverL * L;

        double Vyi = Py * (1.0 - a/L);
        double Vyj = Py * aOverL;

        double Vzi = Pz * (1.0 - a/L);
        double Vzj = Pz * aOverL;

        for (int i = 0; i < nsr; i++) {
          if (x <= a) {
            switch (scheme[i]) {
            case FrameStress::N:  s[i] +=       N; break;
            case FrameStress::Vy: s[i] -=     Vyi; break;
            case FrameStress::Vz: s[i] -=     Vzi; break;
            case FrameStress::My: s[i] += x * Vzi; break;
            case FrameStress::Mz: s[i] -= x * Vyi; break;
            default:                  break;
            }
          } else {
            switch (scheme[i]) {
            case FrameStress::Vy: s[i] +=           Vyj; break;
            case FrameStress::Vz: s[i] +=           Vzj; break;
            case FrameStress::My: s[i] += (L - x) * Vzj; break;
            case FrameStress::Mz: s[i] -= (L - x) * Vyj; break;
            default:                  break;
            }
          }
        }
        break;
      }
    }
  }

  template <int nsr, const FrameStressLayout& scheme>
  void addBasicTangent(MatrixND<nsr,3>& Ks,
                       VectorND<nsr>& s) const
  {

    VectorND<3> nm = s.template extract<3>(0), 
                M  = s.template extract<3>(3);

    double scale = pattern.getLoadFactor();

    Matrix3D Px = Hat(nm);
    if (basis == Director) {
      MatrixND<6,3> K{};
      K.assemble(        Px,   0, 0, -scale);
      K.assemble(    Hat(M),   3, 0, -scale);

      for (int i = 0; i < nsr; i++) {
        MatrixND<1,3> ki = K.template extract<1,3>(i, 0);
        switch (scheme[i]) {
        case FrameStress::N:  Ks.assemble(ki, i, 0, 1.0); break;
        case FrameStress::Vy: Ks.assemble(ki, i, 0, 1.0); break;
        case FrameStress::Vz: Ks.assemble(ki, i, 0, 1.0); break;
        case FrameStress::My: Ks.assemble(ki, i, 0, 1.0); break;
        case FrameStress::Mz: Ks.assemble(ki, i, 0, 1.0); break;
        default:
          break;
        }
      }
    }
  }

  void
  addLinearSolution(double*p0, double L, const Matrix3D& R0, const Matrix3D& R)
  {
    double scale = pattern.getLoadFactor();
    Vector3D n{}, M{};
    {
      Vector3D rx = r[0];
      rx[0] = 0.0;
      rx = R*(R0*rx);
      switch (basis) {
        case Embedding:
          n = R^p[0];
          break;
        case Reference:
          n = R^(R0^p[0]);
          M = (R^m[0]) + rx.cross(n);
          break;
        case Director:
          n = p[0];
          break;
      }
    }

    if (shape == Heaviside) {
      double wa = n[0] * scale; // Axial
      double wy = n[1] * scale; // Transverse
      double wz = n[2] * scale; // Transverse

      double P  =     wa*L;
      double Vy = 0.5*wy*L;
      double Vz = 0.5*wz*L;
      double My = M[1]*scale;
      double Mz = M[2]*scale;

      // Reactions in basic system (projections on linear shape functions)
      p0[0] -=  P;
      p0[1] +=  Mz - Vy; // Vyi
      p0[2] += -Mz - Vy; // Vyj
      p0[3] +=  My - Vz;
      p0[4] += -My - Vz;
    }

    #if 0
    else if (shape == LOAD_TAG_Beam3dPartialUniformLoad) {
      double wy  = p[0](1) * scale;  // Transverse Y at start
      double wz  = p[0](2) * scale;  // Transverse Z at start
      double wa  = p[0](0) * scale;  // Axial at start
      double a   = data(3) * L;
      double b   = data(4) * L;
      double wyb = data(5) * scale;  // Transverse Y at end
      double wzb = data(6) * scale;  // Transverse Z at end
      double wab = data(7) * scale;  // Axial at end
      p0[0] -= wa * (b - a) + 0.5 * (wab - wa) * (b - a);
      double c = a + 0.5 * (b - a);
      double Fy = wy * (b - a); // resultant transverse load Y (uniform part)
      p0[1] -= Fy * (1 - c / L);
      p0[2] -= Fy * c / L;
      double Fz = wz * (b - a); // resultant transverse load Z (uniform part)
      p0[3] -= Fz * (1 - c / L);
      p0[4] -= Fz * c / L;
      c = a + 2.0 / 3.0 * (b - a);
      Fy = 0.5 * (wyb - wy) * (b - a); // resultant transverse load Y (triang. part)
      p0[1] -= Fy * (1 - c / L);
      p0[2] -= Fy * c / L;
      Fz = 0.5 * (wzb - wz) * (b - a); // resultant transverse load Z (triang. part)
      p0[3] -= Fz * (1 - c / L);
      p0[4] -= Fz * c / L;
    }
    #endif

    else if (shape == Dirac) {
      double N      = p[0](0) * scale;
      double Py     = p[0](1) * scale;
      double Pz     = p[0](2) * scale;
      double aOverL = r[0][0];

      if (aOverL < 0.0 || aOverL > 1.0)
        return;

      double V1 = Py * (1.0 - aOverL);
      double V2 = Py * aOverL;
      p0[0] -= N;
      p0[1] -= V1;
      p0[2] -= V2;
      V1 = Pz * (1.0 - aOverL);
      V2 = Pz * aOverL;
      p0[3] -= V1;
      p0[4] -= V2;
    }
  }
private:
  Vector3D getForce(double x,
                    const Matrix3D& R0,
                    const Matrix3D& R) const
  {
    Vector3D n{};
    return n;
  }

  Vector3D getCouple(double x,
                     const Matrix3D& R0,
                     const Matrix3D& R) const
  { 
    Vector3D m{};
    return m;
  }

private:
  const int basis;
  const int shape;
  LoadPattern& pattern;
  std::vector<Vector3D> p;
  std::vector<Vector3D> m;
  std::vector<Vector3D> r;
  std::vector<Element*> elements;
  std::vector<std::array<double,2>> gauss;
};
}

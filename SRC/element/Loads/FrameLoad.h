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
// Type I elements (ExactFrame, CosseratFrame, etc):
//
//   These elements are displacement-interpolated and do not assume a 
//   "basic" system. The interface consists of:
//   - addLoadAtPoint: evaluates the load at a point x along the element
//   - addTangAtPoint
//
// Type II elements (ForceFrame, MixedFrame, etc):
//   These elements are force-interpolated and assume a "basic" system. 
//   They require knowledge of a "particular solution" to the linear BVP.
//   The interface consists of:
//   - addParticularSolution: evaluates the particular solution at a point x along the element
//   - addParticularGradient: evaluates the derivative of the particular solution with respect to nodal variables
//   - addParticularBoundary: evaluates the particular solution at the boundaries of the element
//
//
// Claudio M. Perez
//
#pragma once
#include "Shape.h"
#include <array>
#include <vector>
#include <cassert>
#include <Domain.h>
#include <string.h>
#include <Versor.h>
#include <VectorND.h>
#include <Matrix3D.h>
#include <GroupSO3.h>
#include <FiniteElement.h>
#include <ElementalLoad.h>
#include <StaticPattern.h>
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
  FrameLoad(int tag,
            int basis, 
            int shape, 
            const std::vector<Vector3D>& p,
            const std::vector<Vector3D>& m,
            const std::vector<Vector3D>& r,
            StaticPattern& pattern)
  : ElementalLoad(tag,classTag),
    basis(basis),
    shape(shape),
    pattern(pattern),
    p(p),
    m(m),
    r(r)
  {
    assert(r.size() == 1);
    assert(p.size() == 1);
    assert(m.size() == 1);

    switch (shape) {
      case Dirac:
          gauss = {{r[0][0], 1.0}};
          break;
      case Lagrange:
        // Unimplemented
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

  int
  setDomain(Domain *theDomain) final
  {
    this->ElementalLoad::setDomain(theDomain);
  
    if (theDomain == nullptr) {
      for (auto e: elements)
        e->addLoad(this, 0.0);
      return 0;
    }
    return 0;
  }

  void Print(OPS_Stream &s, int flag) final {}

  int 
  addElement(Element& element) 
  {
    auto name = element.getClassType();
    if (strstr(name, "Frame") == nullptr) {
      opserr << "WARNING FrameLoad::addElement - cannot add load to element of type " << name << '\n';
      return -1;
    }
    elements.push_back(&element);
    element.addLoad(this, 1.0);
    return 0;
  }

  const std::vector<std::array<double,2>>& 
  quadrature () {
    return gauss;
  }

  void
  applyLoad(double loadFactor) final {
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
    return false; //(basis != Director) && (r[0][1] == 0.0) && (r[0][2] == 0.0);
  }

  bool proportional() const {
    return (basis != Director); // && (r[0][1] == 0.0) && (r[0][2] == 0.0);
  }


public:
  template <int i, int nn, int ndf>
  void addLoadAtPoint(VectorND<nn*ndf>& pe, 
                      double x, double w, double jxs,
                      const Matrix3D& R0,
                      const Matrix3D& R) const
  {
    // NOTE: Here R is the pure global rotation field taken directly from the nodes;
    // it does not include the element orientation.
    if (w == 0.0)
      return;

    // const Vector3D theta = LogSO3(R);
    // const Matrix3D T = dLogSO3(theta);// TanSO3(theta);//
    for (unsigned q = 0; q < r.size(); q++) {
      Vector3D px,mx,rx;
      rx = r[q];
      rx[0] = 0.0;
      rx = R*(R0*rx);
      switch (basis) {
        case Embedding:
            px = p[q];
            mx = (m[q] + rx.cross(px));
            break;
        case Reference:
            px = R0 * p[q];
            mx = R0 * m[q] + rx.cross(px);
            break;
        case Director:
            px = R * R0*p[q];
            mx = R * R0*m[q] + rx.cross(px);
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

    // const Vector3D theta = LogSO3(R);
    // const Matrix3D T = TanSO3(theta);

    for (unsigned q = 0; q < r.size(); q++) {
      Vector3D px, mx, rx;
      rx = r[q];
      rx[0] = 0.0;

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
          px = R * R0*p[q];
          mx = R * R0*m[q] + rx.cross(px);
          break;
      }
      // const Matrix3D dT = dTanSO3(theta, mx);

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
        K.assemble(Px*Hat(rx), 3+i*n, 3+j*n,  scale);
    }
  }

  //
  // Interface for corotational elements
  //
  template <int nsr, const FrameStressLayout& scheme>
  void addParticularSolution(VectorND<nsr>& s, double x, double L,
                             const Matrix3D& R0, const Matrix3D& R) const
  {
    Vector3D wn{}, wm{};
    this->localWrench(R0, R, wn, wm);
    this->template particularSolution<nsr,scheme>(s, x, L, wn, wm);
  }


  template <int NDF>
  void
  addParticularBoundary(VectorND<NDF*2>& p0, double L, const Matrix3D& R0, const Matrix3D& R)
  const
  {
    // assemble 
    VectorND<6> sx{};
    static constexpr FrameStressLayout scheme = {
      FrameStress::N,
      FrameStress::Vy,
      FrameStress::Vz,
      FrameStress::T,
      FrameStress::My,
      FrameStress::Mz,
    };
    for (int i=0; i<2; i++) {
      sx.zero();
      double sn = i==0 ? 1.0 : -1.0;
      addParticularSolution<6, scheme>(sx, double(i)*L, L, R0, R);
      p0[i*NDF + 0] -= sx[0]*sn; // N
      p0[i*NDF + 1] -= sx[1]*sn; // Vy %%
      p0[i*NDF + 2] -= sx[2]*sn; // Vz %%
      p0[i*NDF + 3] -= sx[3]*sn; // T
    }
    return;
  }


  int 
  addBasicIntegral(VectorND<6>& q0, double L,
                   Frame::Release release,
                   const Matrix3D& R0, 
                   const Matrix3D& R) const
  {
    switch (shape) {
      case Heaviside: {
        double scale = pattern.getLoadFactor();
        double wx = p[0][0] * scale; // Axial
        double wy = p[0][1] * scale; // Transverse
        double wz = p[0][2] * scale; // Transverse

        double P  =     wx*L; // +/- 
        double Vy = 0.5*wy*L;
        double Vz = 0.5*wz*L;
        // Fixed end forces in basic system
        double Mz = Vy/6.0*L; // wy*L*L/12
        double My = Vz/6.0*L; // wz*L*L/12
        q0[0] -= 0.5*P;
        if (!(release.i & Frame::Release::Mz) && 
            !(release.j & Frame::Release::Mz)) {
          q0[1] -= Mz;
          q0[2] += Mz;
        }
        if (release.i & Frame::Release::Mz)
          q0[2] += wy/8*L*L;
          
        if (release.j & Frame::Release::Mz)
          q0[1] -= wy/8*L*L;
        
        if (!(release.i & Frame::Release::My) && 
            !(release.j & Frame::Release::My)) {
          q0[3] += My;
          q0[4] -= My;
        }
        if (release.i & Frame::Release::My)
          q0[4] -= wz/8*L*L;

        if (release.j & Frame::Release::My)
          q0[3] += wz/8*L*L;
      }
    }

    return 0;
  }


  template <int nsr, const FrameStressLayout& scheme>
  void addParticularGradient(MatrixND<nsr,3>& ds, 
                             double x, double L,
                             const Matrix3D& R0, const Matrix3D& R) const
  {

  }

  template <int NDF>
  void addBoundaryGradient(MatrixND<NDF*2,3>& dpf, double L,
                           const Matrix3D& R0, const Matrix3D& R) const
  {

  }


private:

  template <int NDF>
  void addBoundarySpin(MatrixND<NDF*2,3>& dpf, double L,
                       const Matrix3D& R0, const Matrix3D& R) const
  {
  }

  void localWrench(const Matrix3D& R0, const Matrix3D& R,
                   Vector3D& wn, 
                   Vector3D& wm,
                   MatrixND<6,3>* Omega = nullptr) const
  {
    //
    // Local load parameters in the corotated basic frame.
    //
    // On return wn holds the force (per length for Heaviside, total for Dirac)
    // and wm the moment including the eccentricity couple rx x wn, both with
    // components in the basic frame. If Omega is given it receives the
    // derivative of [wn; wm] with respect to the spin dw of the basic frame,
    //
    //     dR = R*Hat(dw)   =>   d[wn; wm] = Omega*dw
    //
    // Only Embedding and Reference loads change when the frame spins; a
    // Director load is constant in the basic frame and Omega is zero.
    //
    const Vector3D rx {0.0, r[0][1], r[0][2]};
    Vector3D ml{};
    switch (basis) {
      case Embedding:
        wn = R^p[0];
        ml = R^m[0];
        break;
      case Reference:
        wn = R^(R0*p[0]);
        ml = R^(R0*m[0]);
        break;
      case Director:
        wn = p[0];
        ml = m[0];
        break;
    }
    wm = ml + rx.cross(wn);
  }

  template <int nsr, const FrameStressLayout& scheme>
  void particularSolution(VectorND<nsr>& s, double x, double L,
                          const Vector3D& wn, 
                          const Vector3D& wm) const
  {
    //
    // Particular solution for load in the basic frame.
    // Linear in (wn, wm); the load factor is applied here.
    //
    const double scale = pattern.getLoadFactor();

    switch (shape) {
      case Heaviside: {
        // wm is moment/length, wn is force/length
        double wa = wn[0]*scale; // Axial
        double wy = wn[1]*scale; // Transverse
        double wz = wn[2]*scale; // Transverse

        for (int i = 0; i < nsr; i++) {
          switch (scheme[i]) {
          case FrameStress::N:  s[i] +=  wa * (L - x); break;
          case FrameStress::Vy: s[i] -=  wm[2]*scale + wy*(x - 0.5*L); break;
          case FrameStress::Vz: s[i] -= -wm[1]*scale + wz*(x - 0.5*L); break;
          case FrameStress::T : s[i] +=  wm[0]*(L-x)*scale; break;
          case FrameStress::My: s[i] += -wz*0.5*x*(x - L); break;
          case FrameStress::Mz: s[i] +=  wy*0.5*x*(x - L); break;
          default:
            break;
          }
        }
        break;
      }

      case Dirac: {
        double N      = wn[0]*scale;
        double Py     = wn[1]*scale;
        double Pz     = wn[2]*scale;
        double T      = wm[0]*scale;
        double My     = wm[1]*scale;
        double Mz     = wm[2]*scale;
        double aOverL = r[0][0];

        if (aOverL < 0.0 || aOverL > 1.0)
          break;

        double a = aOverL * L;

        if (x <= a && x < L) {
          double Vyi = -Py*(1.0 - a/L) + Mz/L;
          double Vzi = -Pz*(1.0 - a/L) - My/L;
          for (int i = 0; i < nsr; i++) {
            switch (scheme[i]) {
            case FrameStress::N:  s[i] +=       N; break;
            case FrameStress::Vy: s[i] -=     Vyi; break;
            case FrameStress::Vz: s[i] -=     Vzi; break;
            case FrameStress::T : s[i] +=       T; break;
            case FrameStress::My: s[i] -= x * Vzi; break;
            case FrameStress::Mz: s[i] += x * Vyi; break;
            default:                  break;
            }
          }
        } else {
          // x > a
          double Vyj = Py * aOverL + Mz/L;
          double Vzj = Pz * aOverL - My/L;
          for (int i = 0; i < nsr; i++) {
            switch (scheme[i]) {
            case FrameStress::Vy: s[i] -=           Vyj; break;
            case FrameStress::Vz: s[i] -=           Vzj; break;
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


private:
  const int basis;
  const int shape;
  StaticPattern& pattern;

  std::vector<Vector3D> p;
  std::vector<Vector3D> m;
  std::vector<Vector3D> r;
  std::vector<Element*> elements;
  std::vector<std::array<double,2>> gauss;
};
}

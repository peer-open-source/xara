//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#pragma once
#include <array>
#include <vector>
#include <Domain.h>
#include <string.h>
#include <Versor.h>
#include <Vector.h> // TODO: this is just for debugging
#include <VectorND.h>
#include <Matrix3D.h>
#include <Rotations.hpp>
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
            gauss = {{0.0, 0.2777777777777778},
                     {0.2777777777777778, 0.4444444444444444},
                     {0.7222222222222222, 0.4444444444444444},
                     {1.0, 0.2777777777777778}};
            if (shape == Heaviside)
                for (unsigned i=0; i<gauss.size(); i++) {
                    gauss[i][0] = r[0][0] + (1.0 - r[0][0]) * gauss[i][0];
                    gauss[i][1] *= (1 - r[0][0]);
                }
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
    addElement(Element& element) {
        auto name = element.getClassType();
        if (strstr(name, "Frame") == nullptr) {
            opserr << "WARNING FrameLoad::addElement() - cannot add load to element of type " << name << '\n';
            return -1;
        }
        elements.push_back(&element);
        element.addLoad(this, 1.0);
        return 0;
    }

    virtual int
    recvSelf(int commitTag, Channel& , FEM_ObjectBroker& ) final
    {
        return -1;
    }

    virtual int
    sendSelf(int commitTag, Channel& ) final
    {
        return -1;
    }

    const std::vector<std::array<double,2>>& 
    quadrature () {
        return gauss;
    }

    virtual void applyLoad(double loadFactor) final {
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


private:
    Vector3D getForce(double x,
                      const Matrix3D& R0,
                      const Matrix3D& R) const
    {
        Vector3D n;
        return n;
    }

    Vector3D getCouple(double x,
                       const Matrix3D& R0,
                       const Matrix3D& R) const
    { 
        Vector3D m;
        return m;
    }

public:
    template <int i, int nn, int n>
    void addLoadAtPoint(VectorND<nn*n>& pe, 
                        double x, double w,
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
                    if (fabs(x - r[q][0]) > 1.0e-6)
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
            pe.template assemble<  i*n>(px, scale);
            pe.template assemble<3+i*n>(mx, scale);
        }
        
    }

    template <int i, int j, int nn, int n>
    void addTangAtPoint(MatrixND<nn*n,nn*n>& K, 
                        double x, double w,
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
                    if (fabs(x - r[q][0]) > 1.0e-6)
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

    template <int n, const FrameStressLayout& scheme>
    void addBasicSolution(VectorND<n>&    s, double x, double L)
    // add particular solution for basic equations
    {
        double scale = pattern.getLoadFactor();

        switch (shape) {
        case LOAD_TAG_Beam3dUniformLoad: {
            double wa = p[0][2] * scale; // Axial
            double wy = p[0][0] * scale; // Transverse
            double wz = p[0][1] * scale; // Transverse

            for (int i = 0; i < n; i++) {
            switch (scheme[i]) {
            case SECTION_RESPONSE_P:  s[i] +=  wa * (L - x); break;
            case SECTION_RESPONSE_VY: s[i] +=  wy * (x - 0.5 * L); break;
            case SECTION_RESPONSE_VZ: s[i] += -wz * (x - 0.5 * L); break;
            case SECTION_RESPONSE_MY: s[i] +=  wz * 0.5 * x * (L - x); break;
            case SECTION_RESPONSE_MZ: s[i] +=  wy * 0.5 * x * (x - L); break;
            default:                  break;
            }
            }
            break;
        }
        case LOAD_TAG_Beam3dPointLoad: {
            double N      = p[0][2] * scale;
            double Py     = p[0][0] * scale;
            double Pz     = p[0][1] * scale;
            double aOverL = r[0][0];

            if (aOverL < 0.0 || aOverL > 1.0)
                break;
    
            double a = aOverL * L;
    
            double Vyi = Py * (1.0 - a/L);
            double Vyj = Py * aOverL;
    
            double Vzi = Pz * (1.0 - a/L);
            double Vzj = Pz * aOverL;
    
            for (int i = 0; i < n; i++) {
              if (x <= a) {
                switch (scheme[i]) {
                case SECTION_RESPONSE_P:  s[i] +=       N; break;
                case SECTION_RESPONSE_VY: s[i] -=     Vyi; break;
                case SECTION_RESPONSE_VZ: s[i] -=     Vzi; break;
                case SECTION_RESPONSE_MY: s[i] += x * Vzi; break;
                case SECTION_RESPONSE_MZ: s[i] -= x * Vyi; break;
                default:                  break;
                }
              } else {
                switch (scheme[i]) {
                case SECTION_RESPONSE_VY: s[i] +=           Vyj; break;
                case SECTION_RESPONSE_VZ: s[i] +=           Vzj; break;
                case SECTION_RESPONSE_MY: s[i] += (L - x) * Vzj; break;
                case SECTION_RESPONSE_MZ: s[i] -= (L - x) * Vyj; break;
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
  LoadPattern& pattern;
  std::vector<Vector3D> p;
  std::vector<Vector3D> m;
  std::vector<Vector3D> r;
  std::vector<Element*> elements;
  std::vector<std::array<double,2>> gauss;
};
}

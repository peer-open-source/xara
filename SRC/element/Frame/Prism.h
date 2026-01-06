#pragma once
#include <optional>
#include <FrameSection.h>

namespace OpenSees {
namespace Frame {
// TODO: Maybe make this public under ElasticFrameSection
struct Prism {
  // n-n
  std::optional<double> A;
  std::optional<double> Ay, Az;
  // m-m
  std::optional<double> Iy, Iz, Iyz;
  // w-w
  std::optional<double> Cw, Ca;
  // n-m
  std::optional<double> Qy, Qz;
  // n-w
  std::optional<double> Rw, Ry, Rz;
  // m-w
  std::optional<double> Sa, Sy, Sz;
  // Derived
  std::optional<double> J;    // Torsional constant


  std::optional<double> E, G;


  double thermal_coeff = 0.0,          // Thermal
          thermal_depth = 0.0;

  Prism() = default;
  Prism(FrameSection& section):
    A{}, Ay{}, Az{},
    Iy{}, Iz{}, Iyz{},
    Cw{}, Ca{},
    Qy{}, Qz{},
    Rw{}, Ry{}, Rz{},
    Sa{}, Sy{}, Sz{},
    J{},
    E{}, G{}
  {

    // 1) Get exact reference properties; not all sections provide these

    double value;
    if (section.getIntegral(Field::Unit,   State::Init, value) == 0)
      A = value;
    else
      A = 1.0;

    if (section.getIntegral(Field::UnitZZ, State::Init, value) == 0)
      Iy = value;

    if (section.getIntegral(Field::UnitYY, State::Init, value) == 0)
      Iz = value;

    // 2) Get Young and Shear Modulus and determine if shear is supported
    // by the section. The shear areas we pull here may still be 
    // uncondensed.
    const ID& layout = section.getType();
    const Matrix& Ks = section.getInitialTangent();
    for (int i=0; i<layout.Size(); i++) {
      if (layout(i) == FrameStress::N && A) {
        E = Ks(i,i)/(*A);
      }
      else if (layout(i) == FrameStress::Vy && A) {
        G = Ks(i,i)/(*A);
        Ay = A;
      }
      else if (layout(i) == FrameStress::Vz && A) {
        G = Ks(i,i)/(*A);
        Az = A;
      }
    }
    //
    for (int i=0; i<layout.Size(); i++) {
      if (layout(i) == FrameStress::My && !Iy && E) {
        Iy = Ks(i,i)/(*E);
      }
      else if (layout(i) == FrameStress::Mz && !Iz && E) {
        Iz = Ks(i,i)/(*E);
      }
    }
    // In a 3D shear-free section, G wouldnt have been found yet
    for (int i=0; i<layout.Size(); i++) {
      if (layout(i) == FrameStress::T && Iy && Iz && !G) {
        G = Ks(i,i)/(*Iy + *Iz);
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

    MatrixND<6,6> Kc = section.getTangent<6,scheme>(State::Init);

    if (G) {
      J  = Kc(3,3)/(*G);
      if (Ay)
        Ay = Kc(1,1)/(*G);
      if (Az)
        Az = Kc(2,2)/(*G);
    }
  }
};
} // namespace Frame
} // namespace OpenSees

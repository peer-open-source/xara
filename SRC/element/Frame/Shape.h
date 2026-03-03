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
#include <optional>
#include <Logging.h>
#include <MatrixND.h>
#include <VectorND.h>

// integrate(Rigid, [AEGM], [012], [012])


namespace OpenSees {

namespace Frame {
struct Release {
  enum {
    None  = 0,
    My    = 1 << 3,
    Mz    = 1 << 4,
  } i,j;
};


struct Shape {
  int ndm, ndf;
  std::optional<double> E, G;
  std::optional<double> density; // per unit length

  double alpha=0.0; // ?
  // Rigid
  // n-n
  double A;// std::optional<double> A;
  // m-m
  std::optional<double> Iy, Iz;
  double Iyz = 0.0;
  // n-m
  double Qy=0.0, Qz=0.0;
  double centroid[2];

  // Warping
  // w-w
  double Cw=0.0, Ca=0.0;
  // n-w
  double Rw=0.0, Ry=0.0, Rz=0.0;
  // m-w
  double Sa=0.0, Sy=0.0, Sz=0.0;

  // Mixed
  std::optional<double> Ay, Az, ky, kz;
  std::optional<double> J;    // Torsional constant
  // double J = consts.Iy + consts.Iz - consts.Ca;

  struct Mixed {
    std::optional<MatrixND<2,2>>  shear_align;
    std::optional<VectorND<2>>    shift_twist;
    std::optional<VectorND<2>>    shift_axial;
    VectorND<2>    shear_origin, twist_origin;
  } mixed;

  enum class MixedType {
    None,
    UT,
    Constant,
    Energetic,
    Equilibrium
  };
  // using MixedType = Shape::MixedType;
  MixedType mixed_form = MixedType::Energetic;

  //
  // Thermal
  double thermal_coeff = 0.0,
         thermal_depth = 0.0;

  // Prism() = default;

  Shape(int ndm, int ndf)
   : ndm(ndm), ndf(ndf)
   {
   }
  
  void print() const {
    opserr << "Shape: ndm=" << ndm << " ndf=" << ndf << "\n"
           << "  A=" << A 
           << " Ay=" << (Ay ? std::to_string(*Ay) : "n/a")
           << " Az=" << (Az ? std::to_string(*Az) : "n/a")
           << " Iy=" << (Iy ? std::to_string(*Iy) : "n/a") 
           << " Iz=" << (Iz ? std::to_string(*Iz) : "n/a") 
           << " Iyz=" << Iyz
           << " Cw=" << Cw << " Ca=" << Ca
           << " Qy=" << Qy << " Qz=" << Qz
           << " Rw=" << Rw << " Ry=" << Ry << " Rz=" << Rz
           << " Sa=" << Sa << " Sy=" << Sy << " Sz=" << Sz
           << " J="  << (J ? std::to_string(*J) : "n/a")
           << " MT = ";
    switch (mixed_form) {
      case MixedType::None:        opserr << "None"; break;
      case MixedType::UT:          opserr << "UT"; break;
      case MixedType::Constant:    opserr << "UG"; break;
      case MixedType::Energetic:   opserr << "UE"; break;
      case MixedType::Equilibrium: opserr << "NR"; break;
    }
    opserr << "\n";
  }

  int validate()
  {

    if (Ay && ky)
      return -1;
    if (Az && kz)
      return -1;

    if (E && *E <= 0.0) {
      opserr << OpenSees::PromptValueError 
             << "Young's modulus must be positive"
             << OpenSees::SignalMessageEnd;
      return -1;
    }

    if (G && *G <= 0.0) {
      opserr << OpenSees::PromptValueError 
             << "Shear modulus must be positive"
             << OpenSees::SignalMessageEnd;
      return -1;
    }

    if (density && *density < 0.0) {
      opserr << OpenSees::PromptValueError 
             << "Density must be non-negative"
             << OpenSees::SignalMessageEnd;
      return -1;
    }

    if (A <= 0.0) {
      opserr << OpenSees::PromptValueError 
             << "Area must be positive"
             << OpenSees::SignalMessageEnd;
      return -1;
    }
    if (Ay && *Ay == 0.0)
      Ay = std::nullopt;
    else if (Ay && *Ay < 0.0) {
      opserr << OpenSees::PromptValueError 
             << "Shear area Ay must be non-negative"
             << OpenSees::SignalMessageEnd;
      return -1;
    }

    if (Az && *Az == 0.0)
      Az = std::nullopt;
    else if (Az && *Az < 0.0) {
      opserr << OpenSees::PromptValueError 
             << "Shear area Az must be non-negative"
             << OpenSees::SignalMessageEnd;
      return -1;
    }

    // We allow zero in Iy and Iz because these may be set to zero
    // when a 2D simulation is run by restraining out-of-plane DOFs
    // in a 3D model; this is the case in a alot of models imported from
    // SAP2000, for example.
    if (Iy && *Iy < 0.0) {
      opserr << OpenSees::PromptValueError 
             << "Moment of inertia Iy must be positive"
             << OpenSees::SignalMessageEnd;
      return -1;
    }
    if (Iz && *Iz < 0.0) {
      opserr << OpenSees::PromptValueError 
             << "Moment of inertia Iz must be positive"
             << OpenSees::SignalMessageEnd;
      return -1;
    }
    return 0;
  }
};
} // namespace Frame
} // namespace OpenSees

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
// Description: This file implements a unified parser for plasticity materials.
//
// Written: cmp
// April 2025
//
#include <set>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <ModelRegistry.h>

#include "isotropy.h"

#include <HardeningMaterial.h>
#include <SimplifiedJ2.h>
#if defined(XARA_HAVE_GENERALIZEDJ2)
#include <GeneralizedJ2.h>
#endif
#include <NonlinearJ2.h>
#include <PlaneStressSimplifiedJ2.h>
#include <J2Plasticity.h>
#include <J2PlasticityThermal.h>
#include <MultiaxialCyclicPlasticity.h>
#include <DruckerPrager.h>
#include <J2BeamFiber2d.h>
#include <J2BeamFiber3d.h>
#include <J2BeamThread3d.h>
#include <J2PlateFiber.h>
#include <J2PlateFibre.h>


using namespace OpenSees;
template <typename Position>
static inline int
ParsePlasticity(ClientData clientData, Tcl_Interp *interp,
                Tcl_Size argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);

  ArgumentTracker<Position> tracker;
  std::set<int> positional;


  // Count the number of required isotropic parameters.
  // This is needed to accommodate uniaxial materials, which generally
  // only require one isotropic parameter (ie, E).
  int niso = (
    (Position::E      < Position::EndRequired) +
    (Position::G      < Position::EndRequired) +
    (Position::Nu     < Position::EndRequired) +
    (Position::K      < Position::EndRequired) +
    (Position::Lambda < Position::EndRequired)
  );

  // 
  // Values we're parsing for
  //
  int tag;
  double density = 0.0;
  // Isotropy
  IsotropicConstants consts {};
  // Plasticity
  double Fy, Fsat = 0, Fo = 0;
  // Hardening
  struct {
    std::vector<double> C{0.0}, gamma{0.0};
  } kinematic;
  struct {
    std::vector<double> Q{0.0}, b{0.0};
  } isotropic;
  struct {
    double speed = 0.0,
           limit = 0.0;
  } overstress;
  double Hiso=0.0;

  struct {
    double theta = 1.0;
    double Hmix  = 0;
  } hard{};
  bool mix = Position::Theta < Position::End;
  // Viscosity
  double eta=0;
  // Drucker-Prager
  double rho = 0, rho_bar = 0;
  double atm = 101.0;
  double delta2 = 0.0;
  double yield_tol = 1e-16;
  int max_iter = 15;

  //
  // 1. Keyword arguments
  //

  // Isotropy
  IsotropicParse iso {consts, niso};
  if (TclCommand_setIsotropicParameters((ClientData)&iso, interp, argc, argv) == TCL_OK) {
    tracker.consume(Position::E);
    tracker.consume(Position::G);
    tracker.consume(Position::Nu);
    tracker.consume(Position::K);
    tracker.consume(Position::Lambda);
  }

  // Other arguments
  for (int i=2; i<argc; i++) {
    if (iso.positions.find(i) != iso.positions.end()) {
      continue;
    }

    else if (strcmp(argv[i], "-rho") == 0 || 
             strcmp(argv[i], "-density") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &density) != TCL_OK) {
        opserr << "Invalid density value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
    }
    // Yielding
    else if (strcmp(argv[i], "-Fy") == 0 || 
             strcmp(argv[i], "-fy") == 0 ||
             strcmp(argv[i], "-yield-stress") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptParseError
               << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fy) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid yield stress value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::YieldStress);
    }
    else if ((strcmp(argv[i], "-Fo") == 0) || 
             (strcmp(argv[i], "-Ko") == 0)) {
      if (++i >= argc) {
        opserr << OpenSees::PromptParseError
               << "Missing value for option " << argv[i-1] 
               << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fo) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid initial saturation stress value " 
               << argv[i] 
               << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::SatStress0);
    }
    else if ((strcmp(argv[i], "-Fs") == 0) || 
             (strcmp(argv[i], "-Fsat") == 0) ||
             (strcmp(argv[i], "-fsat") == 0)) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fsat) != TCL_OK) {
        opserr << "Invalid saturation stress value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      isotropic.Q[0] = Fsat - Fy;
      tracker.consume(Position::SatStress);
    }
    else if (strcmp(argv[i], "-overstress") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &overstress.limit) != TCL_OK) {
        opserr << "Invalid overstress limit value " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
    }
    else if (strcmp(argv[i], "-transition") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &overstress.speed) != TCL_OK) {
        opserr << "Invalid overstress speed value " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
    }
    else if (strcmp(argv[i], "-tol") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &yield_tol) != TCL_OK) {
        opserr << "Invalid yield tolerance value " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
    }
    //
    // Hardening
    //
    else if (strcmp(argv[i], "-Q") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      int argc_iso;
      TCL_Char** argv_iso;
      if (Tcl_SplitList(interp, argv[i], &argc_iso, &argv_iso) != TCL_OK) {
        if (Tcl_GetDouble(interp, argv[i], &isotropic.Q[0]) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
        }
      }
      else {
        isotropic.Q.resize(argc_iso);
        for (int j=0; j<argc_iso; j++) {
          double val;
          if (Tcl_GetDouble(interp, argv_iso[j], &val) != TCL_OK) {
            opserr << "Invalid value for option " << argv[i-1] << "\n";
            return TCL_ERROR;
          }
          isotropic.Q[j] = val;
        }
        Tcl_Free((char*)argv_iso);
      }
    }
    else if ((strcmp(argv[i], "-b") == 0) || 
             (strcmp(argv[i], "-delta") == 0) || 
             (strcmp(argv[i], "-Hsat") == 0)) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      int argc_iso;
      TCL_Char** argv_iso;
      if (Tcl_SplitList(interp, argv[i], &argc_iso, &argv_iso) != TCL_OK) {
        if (Tcl_GetDouble(interp, argv[i], &isotropic.b[0]) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "Invalid value for option " << argv[i-1] 
                 << " " << argv[i]
                 << "\n";
          return TCL_ERROR;
        }
      }
      else {
        isotropic.b.resize(argc_iso);
        for (int j=0; j<argc_iso; j++) {
          double val;
          if (Tcl_GetDouble(interp, argv_iso[j], &val) != TCL_OK) {
            opserr << "Invalid value for option " << argv[i-1] << "\n";
            return TCL_ERROR;
          }
          isotropic.b[j] = val;
        }
        Tcl_Free((char*)argv_iso);
      }
      tracker.consume(Position::Hsat);
    }
    else if ((strcmp(argv[i], "-C") == 0) || 
             (strcmp(argv[i], "-Hkin") == 0) || 
             (strcmp(argv[i], "-kinematic-hardening") == 0)) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      
      int argc_kin;
      TCL_Char** argv_kin;
      if (Tcl_SplitList(interp, argv[i], &argc_kin, &argv_kin) != TCL_OK) {
        if (Tcl_GetDouble(interp, argv[i], &kinematic.C[0]) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
        }
      }
      else {
        kinematic.C.resize(argc_kin);
        for (int j=0; j<argc_kin; j++) {
          double val;
          if (Tcl_GetDouble(interp, argv_kin[j], &val) != TCL_OK) {
            opserr << "Invalid value for option " << argv[i-1] << "\n";
            return TCL_ERROR;
          }
          kinematic.C[j] = val;
        }
        Tcl_Free((char*)argv_kin);
      }
      mix = false;
      tracker.consume(Position::Hkin);
    }
    else if (strcmp(argv[i], "-gamma") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      int argc_kin;
      TCL_Char** argv_kin;
      if (Tcl_SplitList(interp, argv[i], &argc_kin, &argv_kin) != TCL_OK) {
        if (Tcl_GetDouble(interp, argv[i], &kinematic.gamma[0]) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
        }
      }
      else {
        kinematic.gamma.resize(argc_kin);
        for (int j=0; j<argc_kin; j++) {
          double val;
          if (Tcl_GetDouble(interp, argv_kin[j], &val) != TCL_OK) {
            opserr << "Invalid value for option " << argv[i-1] << "\n";
            return TCL_ERROR;
          }
          kinematic.gamma[j] = val;
        }
        Tcl_Free((char*)argv_kin);
      }
    }
    else if (strcmp(argv[i], "-Hiso") == 0 || 
             strcmp(argv[i], "-isotropic-hardening") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Hiso) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      mix = false;
      tracker.consume(Position::Hiso);
    }
    // else if (strcmp(argv[i], "-Hkin") == 0 || 
    //          strcmp(argv[i], "-kinematic-hardening") == 0) {
    //   if (++i >= argc) {
    //     opserr << "Missing value for option " << argv[i-1] << "\n";
    //     return TCL_ERROR;
    //   }
    //   if (Tcl_GetDouble(interp, argv[i], &kinematic.C[0]) != TCL_OK) {
    //     opserr << "Invalid value for option " << argv[i-1] << "\n";
    //     return TCL_ERROR;
    //   }
    //   mix = false;
    //   tracker.consume(Position::Hkin);
    // }
    else if (strcmp(argv[i], "-H") == 0 || strcmp(argv[i], "-Hmix") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &hard.Hmix) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      mix = true;
      tracker.consume(Position::Hmix);
      tracker.consume(Position::Hiso);
      tracker.consume(Position::Hkin);
    }

    else if (strcmp(argv[i], "-theta") == 0 || strcmp(argv[i], "-mix") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &hard.theta) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Theta);
    }
    // else if (strcmp(argv[i], "-Hsat") == 0 || 
    //          strcmp(argv[i], "-delta") == 0 || 
    //          strcmp(argv[i], "-delta1") == 0) {
    //   if (++i >= argc) {
    //     opserr << "Missing value for option " << argv[i-1] << "\n";
    //     return TCL_ERROR;
    //   }
    //   if (Tcl_GetDouble(interp, argv[i], &isotropic.b[0]) != TCL_OK) {
    //     opserr << "Invalid value for option " << argv[i-1] << "\n";
    //     return TCL_ERROR;
    //   }
    //   tracker.consume(Position::Hsat);
    // }
    // Drucker-Prager
    else if (strcmp(argv[i], "-delta2") == 0 ||
             strcmp(argv[i], "-Hten") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &delta2) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Delta2);
    }
    else if (strcmp(argv[i], "-atm") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &atm) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Atm);
    }
    else if (strcmp(argv[i], "-Rvol") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &rho) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Rho);
    }
    else if (strcmp(argv[i], "-Rbar") == 0 ||
             strcmp(argv[i], "-rhoBar") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &rho_bar) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::RhoBar);
    }
    // Viscosity
    else if (strcmp(argv[i], "-eta") == 0 || strcmp(argv[i], "-viscosity") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &eta) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Eta);
    }
    else
      positional.insert(i);
  }

  //
  // 2) Positional arguments
  //
  for (int i : positional) {
  
    if (tracker.current() == Position::EndRequired)
      tracker.increment();

    switch (tracker.current()) {
      // General
      case Position::Tag :
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid section Elastic tag.\n";
            return TCL_ERROR;           
        } else {
          tracker.increment();
          break;
        }
      case Position::Density:
        if (Tcl_GetDouble (interp, argv[i], &density) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid density.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      // Isotropy
      case Position::E:
        if (Tcl_GetDouble (interp, argv[i], &consts.E) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid E.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      case Position::G:
        if (Tcl_GetDouble (interp, argv[i], &consts.G) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid G.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      case Position::K:
        if (Tcl_GetDouble (interp, argv[i], &consts.K) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid K.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      case Position::Nu:
        if (Tcl_GetDouble (interp, argv[i], &consts.nu) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid nu.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      
      case Position::Lambda:
        if (Tcl_GetDouble (interp, argv[i], &consts.lambda) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Lame lambda.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      // Yielding
      case Position::YieldStress:
        if (Tcl_GetDouble(interp, argv[i], &Fy) != TCL_OK) {
            opserr << OpenSees::PromptParseError 
                   << "invalid yield stress."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;           
        } else {
          tracker.increment();
          break;
        }
      case Position::SatStress:
        if (Tcl_GetDouble (interp, argv[i], &Fsat) != TCL_OK) {
          opserr << OpenSees::PromptParseError 
                  << "invalid saturation stress " 
                  << argv[i]
                  << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      
      case Position::SatStress0:
        if (Tcl_GetDouble (interp, argv[i], &Fo) != TCL_OK) {
          opserr << OpenSees::PromptParseError 
                  << "invalid initial saturation stress.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      // Hardening
      case Position::Hiso:
        if (Tcl_GetDouble (interp, argv[i], &Hiso) != TCL_OK) {
          opserr << OpenSees::PromptParseError 
                  << "invalid Hiso.\n";
          return TCL_ERROR;
        } else {
          tracker.consume(Position::Hiso);
          tracker.consume(Position::Theta);
          break;
        }
      case Position::Hkin:
        if (Tcl_GetDouble (interp, argv[i], &kinematic.C[0]) != TCL_OK) {
          opserr << OpenSees::PromptParseError
                 << "invalid Hkin.\n";
          return TCL_ERROR;
        } else {
          tracker.consume(Position::Hkin);
          tracker.consume(Position::Theta);
          break;
        }
      case Position::Hsat:
        if (Tcl_GetDouble (interp, argv[i], &isotropic.b[0]) != TCL_OK) {
          opserr << OpenSees::PromptParseError
                 << "invalid Hsat.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Hmix:
        if (Tcl_GetDouble (interp, argv[i], &hard.Hmix) != TCL_OK) {
          opserr << OpenSees::PromptParseError
                 << "invalid Hmix.\n";
          return TCL_ERROR;
        } else {
          mix = true;
          tracker.consume(Position::Hmix);
          tracker.consume(Position::Hiso);
          tracker.consume(Position::Hkin);
          break;
        }
      case Position::Theta:
        if (Tcl_GetDouble (interp, argv[i], &hard.theta) != TCL_OK) {
          opserr << OpenSees::PromptParseError 
                 << "invalid hardening theta.\n";
          return TCL_ERROR;
        } else {
          tracker.consume(Position::Theta);
          break;
        }

      // Drucker 
      case Position::Delta2:
        if (Tcl_GetDouble (interp, argv[i], &delta2) != TCL_OK) {
          opserr << OpenSees::PromptParseError 
                 << "invalid delta2.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Rho:
        if (Tcl_GetDouble (interp, argv[i], &rho) != TCL_OK) {
          opserr << OpenSees::PromptParseError 
                  << "invalid Rvol.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Atm:
        if (Tcl_GetDouble (interp, argv[i], &atm) != TCL_OK) {
          opserr << OpenSees::PromptParseError 
                  << "invalid atm.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::RhoBar:
        if (Tcl_GetDouble (interp, argv[i], &rho_bar) != TCL_OK) {
          opserr << OpenSees::PromptParseError 
                 << "invalid Rbar.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      // Viscosity
      case Position::Eta:
        if (Tcl_GetDouble (interp, argv[i], &eta) != TCL_OK) {
          opserr << OpenSees::PromptParseError 
                 << "invalid eta.\n";
          return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }

      case Position::EndRequired:
        // This will not be reached
        break;

      case Position::End:
        opserr << OpenSees::PromptParseError 
               << "unexpected argument " << argv[i] << ".\n";
        return TCL_ERROR;
    }
  }

  if (mix) {
    Hiso = hard.theta  * hard.Hmix;
    kinematic.C[0] = (1.0 - hard.theta) * hard.Hmix;
  }

  //
  // 3. Check for required arguments
  //
  if (tracker.current() < Position::EndRequired) {
    opserr << OpenSees::PromptParseError
            << "missing required arguments: ";
    while (tracker.current() != Position::EndRequired) {
      switch (tracker.current()) {
        case Position::Tag :
          opserr << "tag ";
          break;
        // Isotropy
        case Position::E:
          opserr << "E ";
          break;
        case Position::G:
          opserr << "G ";
          break;
        case Position::K:
          opserr << "K ";
          break;
        case Position::Nu:
          opserr << "nu ";
          break;
        case Position::Lambda:
          opserr << "lambda ";
          break;
        // Yielding
        case Position::YieldStress:
          opserr << "Fy ";
          break;
        case Position::SatStress:
          opserr << "Fsat ";
          break;
        case Position::SatStress0:
          opserr << "Fo ";
          break;
        // Hardening
        case Position::Hiso:
          opserr << "Hiso ";
          break;
        case Position::Hkin:
          opserr << "Hkin ";
          break;
        case Position::Hmix:
          opserr << "Hmix ";
          break;
        case Position::Theta:
          opserr << "theta ";
          break;

        // Drucker-Prager
        case Position::Delta2:
          opserr << "Hten ";
          break;
        case Position::Rho:
          opserr << "Rvol ";
          break;
        case Position::RhoBar:
          opserr << "Rbar ";
          break;
        
        // Viscosity
        case Position::Eta:
          opserr << "eta ";
          break;

        case Position::EndRequired:
        case Position::End:
        default:
          break;
      }

      if (tracker.current() == Position::EndRequired)
        break;

      tracker.consume(tracker.current());
    }

    opserr << "\n";

    return TCL_ERROR;
  }

  //
  // Create the material (TODO)
  //
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);
  if ((strcmp(argv[1], "Hardening") == 0) ||
      (strcmp(argv[1], "Steel") == 0)) {

    UniaxialMaterial* theMaterial = new HardeningMaterial(tag, consts.E, Fy, 
                                                          Hiso, kinematic.C[0], eta);
    if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(argv[1], "J2BeamFiber") == 0) && (getenv("XARA_STATIC_MATERIALS") != nullptr)) {
    NDMaterial* theMaterial = nullptr;
    if (builder->getNDM() == 2)
      theMaterial = new J2BeamFiber2d(tag, consts.E, consts.nu, Fy, kinematic.C[0], Hiso);
    else 
      theMaterial = new J2BeamFiber3d(tag, consts.E, consts.nu, Fy, kinematic.C[0], Hiso);

    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }
  else if ((strcmp(argv[1], "J2BeamThread") == 0) ||
           (strcmp(argv[1], "J2BeamFiber") == 0)) {
    NDMaterial* theMaterial = nullptr;
    if (builder->getNDM() == 2)
      theMaterial = new J2BeamFiber2d(tag, consts.E, consts.nu, Fy, kinematic.C[0], Hiso);
    else 
      theMaterial = new J2BeamThread3d(tag, consts.E, consts.nu, 
                                       Fy, kinematic.C[0], Hiso, density);

    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(argv[1], "J2PlateFibre") == 0)) {
    NDMaterial* theMaterial = new J2PlateFibre(tag, consts.E, consts.nu, Fy, kinematic.C[0], Hiso);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(argv[1], "Simplified3DJ2") == 0) ||
           (strcmp(argv[1], "SimplifiedJ2") == 0) ||
           (strcmp(argv[1], "J2Simplified") == 0) ||
           (strcmp(argv[1], "J2L") == 0) ||
           (strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) ||
           (strcmp(argv[1], "3DJ2") == 0)) {

    NDMaterial* theMaterial = new SimplifiedJ2(tag, 3, consts.G, consts.K, 
                                               Fy, kinematic.C[0], 
                                               Hiso, density);
    if (strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) {
      theMaterial = new PlaneStressSimplifiedJ2(tag, 2, *theMaterial);
    }
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(argv[1], "J2") == 0) ||
           (strcmp(argv[1], "J2N") == 0) ||
           (strcmp(argv[1], "J2Plasticity")  == 0)) {

    NDMaterial* theMaterial = new J2Plasticity(tag, 0, consts.K, consts.G, 
                                               Fy, Fsat, 
                                               isotropic.b[0], 
                                               Hiso, eta, density);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if (strcmp(argv[1], "GeneralizedJ2") == 0) {
#if !defined(XARA_HAVE_GENERALIZEDJ2)
    opserr << OpenSees::PromptParseError
           << "GeneralizedJ2 material requires Xara to be built with GeneralizedJ2 support.\n";
    return TCL_ERROR;
#else
    NDMaterial* theMaterial = new GeneralizedJ2(tag,
                                               consts.E, consts.nu,
                                               Fy,
                                               overstress.limit,
                                               Hiso, 
                                               kinematic.C[0],
                                               overstress.speed,
                                               density, 
                                               GeneralizedJ2::HRule::GP);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
#endif
  }

  else if (strcmp(argv[1], "NonlinearJ2") == 0) {
    double b[2]{};
    double Q[2]{};

    if (isotropic.b.size() > 0)
      b[0] = isotropic.b[0];
    if (isotropic.b.size() > 1)
      b[1] = isotropic.b[1];
    if (isotropic.Q.size() > 0)
      Q[0] = isotropic.Q[0];
    if (isotropic.Q.size() > 1)
      Q[1] = isotropic.Q[1];

    NDMaterial* theMaterial = new NonlinearJ2(tag, consts.E, consts.nu,
                                            Fy, density,
                                            Hiso,
                                            b[0], Q[0],
                                            b[1], Q[1],
                                            kinematic.C,
                                            kinematic.gamma,
                                            yield_tol,
                                            max_iter);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(argv[1], "J2PlasticityThermal") == 0) ||
           (strcmp(argv[1], "J2Thermal") == 0)) {
    NDMaterial* theMaterial = new J2PlasticityThermal(tag, 0, consts.K, consts.G, 
              Fy, 
              Fsat,
              isotropic.b[0],
              Hiso, eta, density);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
    delete theMaterial;
    return TCL_ERROR;
    }
    return TCL_OK;
  }
  else if (strcmp(argv[1], "DruckerPrager") == 0 ||
           strcmp(argv[1], "DP") == 0) {

    NDMaterial* theMaterial = new DruckerPrager(
                                    tag, 0, consts.K, consts.G,
                                    Fy, rho, rho_bar, Fsat, Fo,
                                    isotropic.b[0], delta2, hard.Hmix, hard.theta,
                                    density, atm);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }
  return TCL_ERROR;
}


int
TclCommand_newPlasticMaterial(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  // 
  if (strcmp(argv[1], "Simplified3DJ2") == 0 ||
      strcmp(argv[1], "SimplifiedJ2") == 0 ||
      strcmp(argv[1], "J2Simplified") == 0 ||
      strcmp(argv[1], "J2L") == 0 ||
      strcmp(argv[1], "3DJ2") == 0 ||
      strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) {

    // "SimplifiedJ2"  tag?  G?  K?  Fy? Hkin?  Hiso?
    enum class Position : int {
      Tag, G, K, YieldStress, EndRequired,
      Hkin, Hiso,
      End,
      E, Nu, Lambda, Eta, Theta, Hmix, Hsat,
      SatStress, SatStress0,
      Delta2, Rho, RhoBar, Atm,
      Density
    };
    return ParsePlasticity<Position>(clientData, interp, argc, argv);
  }
  else if ((strcmp(argv[1], "J2BeamFiber") == 0) ||
           (strcmp(argv[1], "J2BeamThread") == 0)) {
    // J2BeamFiber $tag $E $v $sigmaY $Hiso $Hkin <$rho>
    enum class Position : int {
      Tag, E, G, YieldStress, EndRequired,
      Hkin, Hiso,
      Density,
      End,
      Nu, K, Eta, Lambda, Theta, Hmix, Hsat,
      SatStress, SatStress0,
      Delta2, Rho, RhoBar, Atm
    };
    return ParsePlasticity<Position>(clientData, interp, argc, argv);
  }
  else if ((strcmp(argv[1], "J2PlateFibre") == 0)) {
    // J2PlateFibre $tag $E $v $sigmaY $Hiso $Hkin <$rho>
    enum class Position : int {
      Tag, E, G, YieldStress, EndRequired,
      Hkin, Hiso,
      Density,
      End,
      Nu, K, Eta, Lambda, Theta, Hmix, Hsat,
      SatStress, SatStress0,
      Delta2, Rho, RhoBar, Atm
    };
    return ParsePlasticity<Position>(clientData, interp, argc, argv);
  }

  // "UniaxialJ2Plasticity" tag? E? sigmaY? Hkin? <Hiso?>
  else if (strcmp(argv[1], "UniaxialJ2Plasticity") == 0) {
  }

  else if (strcmp(argv[1], "HardeningMaterial") == 0 ||
           strcmp(argv[1], "Hardening")  == 0 ||
           strcmp(argv[1], "Hardening2") == 0 ||
           strcmp(argv[1], "Steel") == 0) {

    // "Hardening"  tag?  E?  Y?  Hiso?  Hkin?
    enum class Position : int {
      Tag, E, YieldStress, Hiso, EndRequired, 
      Hkin,
      End,
      // Keyword-only arguments
      Density,
      // Unused
      Eta, G, K, Nu, Lambda, Theta, Hmix, Hsat,
      SatStress, SatStress0, 
      Delta2, Rho, RhoBar, Atm,
    };
    return ParsePlasticity<Position>(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "J2") == 0 ||
           strcmp(argv[1], "J2N") == 0 ||
           strcmp(argv[1], "J2Plasticity")  == 0) {

    // "J2Plasticity" tag? K? G? sig0? sigInf? delta? Hiso? <eta?>
    enum class Position : int {
      Tag, K, G, YieldStress, SatStress, Hsat, Hiso, EndRequired, 
      Eta,                                           End,
      E, Nu, Lambda, Hkin, Theta, Hmix, SatStress0, 
      Delta2, Rho, RhoBar, Atm,
      Density
    };
    return ParsePlasticity<Position>(clientData, interp, argc, argv);
  }
  else if (strcmp(argv[1], "GeneralizedJ2") == 0) {

    // "GeneralizedJ2" tag? E? nu? sig0?  Hiso? Hkin? sigInf? <density?>
    enum class Position : int {
      Tag, E, Nu, YieldStress, EndRequired, 
      Hiso, Hkin,  SatStress, Density,
      End,
      G, K, Lambda, Eta, Theta, Hmix, Hsat, SatStress0,
      Delta2, Rho, RhoBar, Atm
    };
    return ParsePlasticity<Position>(clientData, interp, argc, argv);
  }
  else if (strcmp(argv[1], "NonlinearJ2") == 0) {

    // "NonlinearJ2" tag? E? nu?
    enum class Position : int {
      Tag, E, Nu, YieldStress, EndRequired, 
      End,
      G, K, Lambda, Eta, Theta, Hmix, Hsat, SatStress0,
      Delta2, Rho, RhoBar, Atm,
      Hiso, Hkin,  SatStress, Density,
    };
    return ParsePlasticity<Position>(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "J2PlasticityThermal") == 0 ||
           strcmp(argv[1], "J2Thermal") == 0) {

    // "J2Thermal" tag? K? G? sig0? sigInf? delta? Hiso? <eta?>
    enum class Position : int {
      Tag, K, G, YieldStress, SatStress, Hsat, Hiso, EndRequired,
      Eta,                                           End,
      E, Nu, Lambda, Hkin, Theta, Hmix, SatStress0,
      Delta2, Rho, RhoBar, Atm,
      Density
    };
    return ParsePlasticity<Position>(clientData, interp, argc, argv);
  }
  else if (strcmp(argv[1], "DP") == 0 ||
           strcmp(argv[1], "DruckerPrager")  == 0) {

    // DruckerPrager tag? K? G? sigma_y? rho? rho_bar? Kinf? Ko? delta1? delta2? H? theta? <massDensity? atm?>
    enum class Position : int {
      Tag, K, G, YieldStress, Rho, RhoBar, 
        SatStress, SatStress0, Hsat, Delta2, Hmix, Theta, EndRequired, 
      Density, Atm, End,
      Eta, E, Nu, Lambda, Hiso, Hkin
    };
    return ParsePlasticity<Position>(clientData, interp, argc, argv);
  }
  return TCL_ERROR;
}


#include <UniaxialJ2Plasticity.h>
int
TclCommand_newUniaxialJ2Plasticity(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char** const argv)
{
    // ----- 1D J2 Plasticity ----
    if (argc < 7) {
      opserr << "WARNING invalid number of arguments\n";
      opserr << "Want: uniaxialMaterial UniaxialJ2Plasticity tag? E? sigmaY? Hkin? <Hiso?>"
             << "\n";
      return TCL_ERROR;
    }

    int tag;
    double E, sigmaY, Hkin, Hiso;
    Hiso = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial UniaxialJ2Plasticity tag"
             << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &sigmaY) != TCL_OK) {
      opserr << "WARNING invalid sigmaY\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &Hkin) != TCL_OK) {
      opserr << "WARNING invalid Hkin\n";
      return TCL_ERROR;
    }

    if (argc >= 7)
      if (Tcl_GetDouble(interp, argv[6], &Hiso) != TCL_OK) {
        opserr << "WARNING invalid Hiso\n";
        return TCL_ERROR;
      }

    // Parsing was successful, allocate the material
    UniaxialMaterial* theMaterial = new UniaxialJ2Plasticity(tag, E, sigmaY, Hkin, Hiso);

   assert(clientData != nullptr);
   ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);
   builder->addTaggedObject<UniaxialMaterial>(*theMaterial);
   return TCL_OK;
}


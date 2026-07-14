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
// Written: cmp
// April 2025
//
#include <set>
#include <string.h>
#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#include <cmath>
#include <limits>
#include <algorithm>


#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <ModelRegistry.h>
#include "isotropy.h"

#include <damage/FariaPlasticDamage3d.h>

static int find_faria_peak(double Fc, double An, double Bn, double& Fo);

int
TclCommand_newConcreteMaterial(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);
  enum class Position {
    Tag, E, Nu, PeakTension, PeakCompression, EndRequired,
    Beta, Ap, An, Bn,
    Density,
    G, K, Lambda,
    End,
  };
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
  double Fc, Ft=0;
  // Damage
  double beta = 0.6;
  double Ap = 0.5,
         An = 2.0,
         Bn = 0.75;
  bool scale_peak = false;
  

  //
  // 1. Keyword arguments
  //

  // Isotropy
  IsotropicParse iso {consts, niso};
  if (XaraCmd_setIsotropicParameters((ClientData)&iso, interp, argc, argv) == TCL_OK) {
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

    else if (strcmp(argv[i], "-rho") == 0 || strcmp(argv[i], "-density") == 0) {
      if (++i >= argc) {
          opserr << OpenSees::PromptValueError
                 << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &density) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "Invalid density " << argv[i] << "\n";
          return TCL_ERROR;
      }
    }
  
    else if ((strcmp(argv[i], "-scale_peak") == 0) || (strcmp(argv[i], "-scale-peak") == 0)) {
      scale_peak = true;
    }
  
    // Compression
    else if (strcmp(argv[i], "-Fc") == 0 || 
             strcmp(argv[i], "-fc") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fc) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid " << &argv[i-1][1] << " value " << argv[i] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::PeakCompression);
    }
    // Tension
    else if ((strcmp(argv[i], "-Ft") == 0)) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Ft) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid " << &argv[i-1][1] << " value " << argv[i] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::PeakTension);
    }

    //
    // Hardening
    //
    else if (strcmp(argv[i], "-beta") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &beta) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Beta);
    }
    //
    //
    else if (strcmp(argv[i], "-Ap") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] 
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Ap) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "Invalid value for option " << argv[i-1] 
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      tracker.consume(Position::Ap);
    }
    else if (strcmp(argv[i], "-An") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
               << "Missing value for option " << argv[i-1] 
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &An) != TCL_OK) {
        opserr << OpenSees::PromptValueError
                << "Invalid value for option " << argv[i-1] 
                << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      tracker.consume(Position::An);
    }
    else if (strcmp(argv[i], "-Bn") == 0) {
      if (++i >= argc) {
        opserr << OpenSees::PromptValueError
                << "Missing value for option " << argv[i-1] 
                << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Bn) != TCL_OK) {
        opserr << OpenSees::PromptValueError
                << "Invalid value for option " << argv[i-1] 
                << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      tracker.consume(Position::Bn);
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
        }
        tracker.increment();
        break;

      case Position::Density:
        if (Tcl_GetDouble (interp, argv[i], &density) != TCL_OK) {
            opserr << OpenSees::PromptParseError 
                   << "invalid density" << argv[i]
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;

      // Isotropy
      case Position::E:
        if (Tcl_GetDouble (interp, argv[i], &consts.E) != TCL_OK) {
            opserr << OpenSees::PromptParseError 
                   << "invalid E" << argv[i]
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;

      case Position::G:
        if (Tcl_GetDouble (interp, argv[i], &consts.G) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid G."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;

      case Position::K:
        if (Tcl_GetDouble (interp, argv[i], &consts.K) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid K."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;

      case Position::Nu:
        if (Tcl_GetDouble (interp, argv[i], &consts.nu) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid nu."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;
      
      case Position::Lambda:
        if (Tcl_GetDouble (interp, argv[i], &consts.lambda) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Lame lambda."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;
    
      // Yielding
      case Position::PeakCompression:
        if (Tcl_GetDouble(interp, argv[i], &Fc) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Fc."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;           
        }
        tracker.increment();
        break;

      case Position::PeakTension:
        if (Tcl_GetDouble (interp, argv[i], &Ft) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Ft."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Position::Beta:
        if (Tcl_GetDouble (interp, argv[i], &beta) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid beta."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Position::Ap:
        if (Tcl_GetDouble (interp, argv[i], &Ap) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Ap."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Position::An:
        if (Tcl_GetDouble (interp, argv[i], &An) != TCL_OK) {
            opserr << OpenSees::PromptParseError 
                   << "invalid An."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Position::Bn:
        if (Tcl_GetDouble (interp, argv[i], &Bn) != TCL_OK) {
            opserr << OpenSees::PromptParseError 
                   << "invalid Bn."
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
        }
        tracker.increment();
        break;
    

      case Position::EndRequired:
        // This will not be reached
        break;

      case Position::End:
        opserr << OpenSees::PromptParseError 
               << "unexpected argument " << argv[i] << "."
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
    }
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
        case Position::PeakCompression:
          opserr << "Fc ";
          break;
        case Position::PeakTension:
          opserr << "Ft ";
          break;
        case Position::Ap:
          opserr << "Ap ";
          break;
        // Hardening
        case Position::An:
          opserr << "An ";
          break;
        case Position::Bn:
          opserr << "Bn ";
          break;
        case Position::Beta:
          opserr << "beta ";
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

    opserr << OpenSees::SignalMessageEnd;

    return TCL_ERROR;
  }

  //
  // Validate and normalize parameters
  //
  Fc = std::abs(Fc);


  if (scale_peak) {
    double Fcy;
    if (find_faria_peak(Fc, An, Bn, Fcy) != 0) {
      opserr << OpenSees::PromptParseError
             << "Failed to find yield from peak stress.\n";
      return TCL_ERROR;
    }
    Fc = Fcy;
  }

  //
  // Create the material
  //
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);
  if ((strcmp(argv[1], "PlasticDamageConcrete3d") == 0) ||
      (strcasecmp(argv[1], "PlasticDamageConcrete") == 0) ||
      (strcmp(argv[1], "FariaPlasticDamage") == 0)) {

    NDMaterial* theMaterial = new FariaPlasticDamage3d(tag,
                                                       consts.E, consts.nu, 
                                                       Ft, Fc, 
                                                       beta, Ap, An, Bn, 
                                                       density);

    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  return TCL_ERROR;
}

#include <cmath>
#include <limits>
#include <algorithm>

static int 
find_faria_peak(double Fc, double An, double Bn, double& Fo)
{
  Fo = std::numeric_limits<double>::quiet_NaN();

  if (!std::isfinite(Fc) || !std::isfinite(An) || !std::isfinite(Bn))
    return -1;

  if (Fc <= 0.0 || An <= 0.0 || Bn <= 0.0)
    return -1;

  const double A = An;
  const double B = Bn;
  const double sqrt2 = std::sqrt(2.0);
  const double tolA = 1.0e-12;
  const int max_iter = 100;

  auto m = [A, B](double x) -> double {
    return (1.0 - A)*x + A*x*x*std::exp(B*(1.0 - x));
  };

  auto dm = [A, B](double x) -> double {
    return (1.0 - A)  + A*std::exp(B*(1.0 - x))*(2.0*x - B*x*x);
  };

  auto bisect_pos_to_neg = [&](double lo, double hi) -> double {
    // Assumes dm(lo) >= 0 and dm(hi) <= 0.
    for (int i = 0; i < max_iter; ++i) {
      const double mid = 0.5*(lo + hi);
      if (dm(mid) > 0.0)
        lo = mid;
      else
        hi = mid;
    }
    return 0.5*(lo + hi);
  };

  double M = std::numeric_limits<double>::quiet_NaN();

  // Case A = 1.
  if (std::abs(A - 1.0) <= tolA*std::max(1.0, std::abs(A))) {
    if (B < 2.0)
      M = 4.0/(B*B)*std::exp(B - 2.0);
    else
      M = 1.0;
  }

  // Case 0 < A < 1.
  // The true global envelope is unbounded as x -> infinity.
  // This returns the first local peak, which is usually the calibration target.
  else if (A < 1.0) {
    if (dm(1.0) <= 0.0) {
        M = 1.0;
    }
    else {
      const double xlo = std::max(1.0, 2.0/B);
      const double xcrit = (2.0 + sqrt2)/B;

      if (!(xcrit > xlo))
        return -1;

      if (!(dm(xlo) >= 0.0 && dm(xcrit) < 0.0))
        return -1;

      const double xpeak = bisect_pos_to_neg(xlo, xcrit);
      M = m(xpeak);
    }
  }

  // Case A > 1.
  // This has a finite global peak.
  else {
    if (B >= 2.0) {
      M = 1.0;
    }
    else {
      const double xright = 2.0/B;
      const double xcrit = std::max(1.0, (2.0 - sqrt2)/B);

      if (dm(xcrit) <= 0.0)
        M = 1.0;
      else {
        const double xpeak = bisect_pos_to_neg(xcrit, xright);
        M = std::max(1.0, m(xpeak));
      }
    }
  }

  if (!std::isfinite(M) || M <= 0.0)
    return -1;

  Fo = Fc / M;
  return 0;
}

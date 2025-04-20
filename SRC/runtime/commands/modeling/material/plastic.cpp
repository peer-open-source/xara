//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// cmp
//
#include <tcl.h>
#include <set>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include "isotropy.h"
#include <BasicModelBuilder.h>
#include <HardeningMaterial.h>

// #include <PlasticMaterial.h>

#include <SimplifiedJ2.h>
#include <PlaneStressSimplifiedJ2.h>
#include <J2Plasticity.h>
#include <MultiaxialCyclicPlasticity.h>



template <typename Position>
int
TclCommand_newPlasticParser(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);

  ArgumentTracker<Position> tracker;
  std::set<int> positional;

  int niso = (
    (Position::E      < Position::EndRequired) +
    (Position::G      < Position::EndRequired) +
    (Position::Nu     < Position::EndRequired) +
    (Position::K      < Position::EndRequired) +
    (Position::Lambda < Position::EndRequired)
  );

  int tag;
  double density = 0.0;
  // Isotropy
  IsotropicConstants consts {};
  // Plasticity
  double Fy, Fsat = 0;
  // Hardening
  double Hiso=0, Hkin=0;
  double delta = 0;
  // Viscosity
  double eta=0;


  // Isotropy
  IsotropicParse iso {consts, niso};
  if (TclCommand_setIsotropicParameters((ClientData)&iso, interp, argc, argv) == TCL_OK) {
    tracker.consume(Position::E);
    tracker.consume(Position::G);
    tracker.consume(Position::Nu);
    tracker.consume(Position::K);
    tracker.consume(Position::Lambda);
  }

  //
  // Plasticity
  //
  for (int i=2; i<argc; i++) {
    if (iso.positions.find(i) != iso.positions.end()) {
      continue;
    }

    else if (strcmp(argv[i], "-rho") == 0 || strcmp(argv[i], "-density") == 0) {
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
             strcmp(argv[i], "-Y") == 0  || 
             strcmp(argv[i], "-yield-stress") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Fy) != TCL_OK) {
          opserr << "Invalid yield stress value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::YieldStress);
    }
    //
    // Hardening
    //
    else if (strcmp(argv[i], "-Hiso") == 0 || strcmp(argv[i], "-isotropic-hardening") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Hiso) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Hiso);
    }
    else if (strcmp(argv[i], "-Hkin") == 0 || strcmp(argv[i], "-kinematic-hardening") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Hkin) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      tracker.consume(Position::Hkin);
    }
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
  // Positional arguments
  //
  for (int i : positional) {
  
    if (tracker.current() == Position::EndRequired)
      tracker.increment();

    switch (tracker.current()) {
      case Position::Tag :
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid section Elastic tag.\n";
            return TCL_ERROR;           
        } else {
          tracker.increment();
          break;
        }
      case Position::YieldStress:
        if (Tcl_GetDouble(interp, argv[i], &Fy) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid yield stress.\n";
            return TCL_ERROR;           
        } else {
          tracker.increment();
          break;
        }
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
      
      // case Position::SatStress:
      //   if (Tcl_GetDouble (interp, argv[i], &Fsat) != TCL_OK) {
      //       opserr << OpenSees::PromptParseError << "invalid saturation stress.\n";
      //       return TCL_ERROR;
      //   } else {
      //     tracker.increment();
      //     break;
      //   }
      
      case Position::Hiso:
        if (Tcl_GetDouble (interp, argv[i], &Hiso) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Hiso.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Hkin:
        if (Tcl_GetDouble (interp, argv[i], &Hkin) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid Hkin.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      case Position::Eta:
        if (Tcl_GetDouble (interp, argv[i], &eta) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid eta.\n";
            return TCL_ERROR;
        } else {
          tracker.increment();
          break;
        }
      // case Position::Density:
      //   if (Tcl_GetDouble (interp, argv[i], &density) != TCL_OK) {
      //       opserr << OpenSees::PromptParseError << "invalid density.\n";
      //       return TCL_ERROR;
      //   } else {
      //     tracker.increment();
      //     break;
      //   }

      case Position::EndRequired:
        // This will not be reached
        break;

      case Position::End:
        opserr << OpenSees::PromptParseError << "unexpected argument " << argv[i] << ".\n";
        return TCL_ERROR;
    }
  }

  if (tracker.current() < Position::EndRequired) {
    opserr << OpenSees::PromptParseError
            << "missing required arguments: ";
    while (tracker.current() != Position::EndRequired) {
      switch (tracker.current()) {
        case Position::Tag :
          opserr << "tag ";
          break;
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
        case Position::YieldStress:
          opserr << "Fy ";
          break;
        case Position::Hiso:
          opserr << "Hiso ";
          break;
        case Position::Hkin:
          opserr << "Hkin ";
          break;
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
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  if ((strcmp(argv[1], "Hardening") == 0) ||
      (strcmp(argv[1], "Steel") == 0)) {

    UniaxialMaterial* theMaterial = new HardeningMaterial(tag, consts.E, Fy, Hiso, Hkin, eta);
    if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(argv[1], "Simplified3DJ2") == 0) ||
           (strcmp(argv[1], "SimplifiedJ2") == 0) ||
           (strcmp(argv[1], "J2Simplified") == 0) ||
           (strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) ||
           (strcmp(argv[1], "3DJ2") == 0)) {

    NDMaterial* theMaterial = new SimplifiedJ2(tag, 3, consts.G, consts.K, Fy, Hkin, Hiso);
    if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
      delete theMaterial;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(argv[1], "J2") == 0) ||
           (strcmp(argv[1], "J2Plasticity")  == 0)) {

    NDMaterial* theMaterial = new J2Plasticity(tag, 0, consts.K, consts.G, Fy, Fsat, delta, Hiso, eta);
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
      strcmp(argv[1], "3DJ2") == 0 ||
      strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) {

    // "SimplifiedJ2"  tag?  G?  K?  Y? $Hkin  $Hiso
    enum class Position : int {
      Tag, G, K, YieldStress, Hkin, Hiso, EndRequired, 
      E, Nu, Lambda, Eta, End
    };
    return TclCommand_newPlasticParser<Position>(clientData, interp, argc, argv);
  }

  // "UniaxialJ2Plasticity" tag? E? sigmaY? Hkin? <Hiso?>
  else if (strcmp(argv[1], "UniaxialJ2Plasticity") == 0) {
  }

  // "Hardening"  tag?  E?  Y?  Hiso?  Hkin?
  else if (strcmp(argv[1], "HardeningMaterial") == 0 ||
           strcmp(argv[1], "Hardening")  == 0 ||
           strcmp(argv[1], "Hardening2") == 0 ||
           strcmp(argv[1], "Steel") == 0) {

    enum class Position : int {
      Tag, E, YieldStress, Hiso, EndRequired, 
      Hkin,
      End,
      Eta, G, K, Nu, Lambda
    };
    return TclCommand_newPlasticParser<Position>(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "J2") == 0 ||
           strcmp(argv[1], "J2Plasticity")  == 0) {

    // "J2Plasticity" tag? K? G? sig0? sigInf? delta? Hiso? <eta?>
    enum class Position : int {
      Tag, K, G, YieldStress, SatStress, Delta, Hiso, EndRequired, 
      Eta,
      E, Nu, Lambda, Hkin, End
    };
    return TclCommand_newPlasticParser<Position>(clientData, interp, argc, argv);
  }

  return TCL_ERROR;
}


int
TclCommand_newJ2Simplified(ClientData clientData, Tcl_Interp* interp, int argc, const char** const argv)
{

  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  if (argc < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: nDmaterial Simplified3DJ2  $matTag  $G  $K  $sig0 $H_kin  $H_iso"
            << "\n";
    return TCL_ERROR;
  }

  int tag;
  double K, G, sig0, H_kin, H_iso;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid tag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[3], &G) != TCL_OK) {
    opserr << "WARNING invalid G\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &K) != TCL_OK) {
    opserr << "WARNING invalid K\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
    opserr << "WARNING invalid sig0\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[6], &H_kin) != TCL_OK) {
    opserr << "WARNING invalid H_kin\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7], &H_iso) != TCL_OK) {
    opserr << "WARNING invalid H_iso\n";
    return TCL_ERROR;
  }

  NDMaterial* theMaterial = nullptr;
  
  if ((strcmp(argv[1], "Simplified3DJ2") == 0) ||
      (strcmp(argv[1], "3DJ2") == 0) ||
      (strcmp(argv[1], "SimplifiedJ2") == 0) ||
      (strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0)) {

    theMaterial = new SimplifiedJ2(tag, 3, G, K, sig0, H_kin, H_iso);

    if (strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) {
      theMaterial = new PlaneStressSimplifiedJ2(tag, 2, *theMaterial);
    }
  }

  //
  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;
}

int
TclCommand_newJ2Material(ClientData clientData,
                         Tcl_Interp* interp,
                         int argc,
                         const char** const argv)
{
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  if (argc < 9) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: nDMaterial J2Plasticity tag? K? G? sig0? sigInf? delta? H? <eta?>" << "\n";
    return TCL_ERROR;
  }

  int tag;
  double K, G, sig0, sigInf, delta, H;
  double eta = 0.0;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid J2Plasticity tag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
    opserr << "WARNING invalid K\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
    opserr << "WARNING invalid G\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
    opserr << "WARNING invalid sig0\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[6], &sigInf) != TCL_OK) {
    opserr << "WARNING invalid sigInf\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7], &delta) != TCL_OK) {
    opserr << "WARNING invalid delta\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[8], &H) != TCL_OK) {
    opserr << "WARNING invalid H\n";
    return TCL_ERROR;
  }
  if (argc > 9 && Tcl_GetDouble(interp, argv[9], &eta) != TCL_OK) {
    opserr << "WARNING invalid eta\n";
    return TCL_ERROR;
  }


  //
  NDMaterial* theMaterial = nullptr;

  if ((strcmp(argv[1], "J2Plasticity") == 0) || 
      (strcmp(argv[1], "J2") == 0)) {
    theMaterial = new J2Plasticity(tag, 0, K, G, sig0, sigInf, delta, H, eta);
  }
  
  if (theMaterial == nullptr)
    return TCL_ERROR;

  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;
}


#if 0
int
TclCommand_newPlasticMaterial(ClientData clientData, Tcl_Interp* interp, int argc, const char**const argv)
{
  if ((strcmp(argv[1], "J2Plasticity") == 0) || (strcmp(argv[1], "J2") == 0))
  {
    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial J2Plasticity tag? K? G? sig0? sigInf? delta? H? <eta?>" << "\n";
      return TCL_ERROR;
    }

    int tag;
    double K, G, sig0, sigInf, delta, H;
    double eta = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid J2Plasticity tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
      opserr << "WARNING invalid K\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
      opserr << "WARNING invalid sig0\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &sigInf) != TCL_OK) {
      opserr << "WARNING invalid sigInf\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &delta) != TCL_OK) {
      opserr << "WARNING invalid delta\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &H) != TCL_OK) {
      opserr << "WARNING invalid H\n";
      return TCL_ERROR;
    }
    if (argc > 9 && Tcl_GetDouble(interp, argv[9], &eta) != TCL_OK) {
      opserr << "WARNING invalid eta\n";
      return TCL_ERROR;
    }

    theMaterial = new J2Plasticity(tag, 0, K, G, sig0, sigInf, delta, H, eta);
  }
}
#endif

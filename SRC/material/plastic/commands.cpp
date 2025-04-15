//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <set>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <isotropy.h>
#include <BasicModelBuilder.h>
// #include <PlasticMaterial.h>

#include <J2Plasticity.h>
#include <MultiaxialCyclicPlasticity.h>

struct IsotropicParse {
  IsotropicConstants &constants;
  std::set<int>       positions;
};

int
TclCommand_setIsotropicParameters(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // We expect the syntax:
    //    material Isotropic $tag <options>
    // The tag is still the first non-option argument (argv[2]).
    // Among the options, exactly two independent elastic parameters must be provided.
    // Valid elastic options: -E, -G, -K, -nu, -lambda.

    bool gotParam1 = false,
         gotParam2 = false;

    int flag1 = 0, 
        flag2 = 0;
    double val1 = 0.0, 
           val2 = 0.0;

    IsotropicParse* data = static_cast<IsotropicParse*>(clientData);
    IsotropicConstants* iso  = &data->constants;
    std::set<int>& positions =  data->positions;

    assert(iso != nullptr);

    // Process the remaining arguments.
    for (int i = 0; i < argc; i++) {
      if (strcmp(argv[i], "-E") == 0 || strcmp(argv[i], "-youngs-modulus") == 0) {
        if (++i >= argc) {
            opserr << "Missing value for option " << argv[i-1] << "\n";
            return TCL_ERROR;
        }
        double val;
        if (Tcl_GetDouble(interp, argv[i], &val) != TCL_OK) {
            opserr << "Invalid value for option " << argv[i-1] << "\n";
            return TCL_ERROR;
        }
        if (!gotParam1) { 
            gotParam1 = true; 
            val1 = val; 
            flag1 = static_cast<int>(Isotropy::Parameter::YoungModulus);
        }
        else if (!gotParam2) { 
            gotParam2 = true; 
            val2 = val; 
            flag2 = static_cast<int>(Isotropy::Parameter::YoungModulus); 
        }
        else {
            opserr << "Too many elastic parameter options provided.\n";
            return TCL_ERROR;
        }
        positions.insert(i-1);
        positions.insert(i);
      }
      else if (strcmp(argv[i], "-G") == 0 || strcmp(argv[i], "-shear-modulus") == 0) {
        if (++i >= argc) {
            opserr << "Missing value for option " << argv[i-1] << "\n";
            return TCL_ERROR;
        }
        double val;
        if (Tcl_GetDouble(interp, argv[i], &val) != TCL_OK) {
            opserr << "Invalid value for option " << argv[i-1] << "\n";
            return TCL_ERROR;
        }
        if (!gotParam1) {
          gotParam1 = true;
          val1 = val;
          flag1 = static_cast<int>(Isotropy::Parameter::ShearModulus);
        }
        else if (!gotParam2) {
          gotParam2 = true;
          val2 = val;
          flag2 = static_cast<int>(Isotropy::Parameter::ShearModulus);
        }
        else {
            opserr << "Too many elastic parameter options provided.\n";
            return TCL_ERROR;
        }
        positions.insert(i-1);
        positions.insert(i);
      }
      else if (strcmp(argv[i], "-K") == 0 || strcmp(argv[i], "-bulk-modulus") == 0) {
        if (++i >= argc) {
            opserr << "Missing value for option " << argv[i-1] << "\n";
            return TCL_ERROR;
        }
        double val;
        if (Tcl_GetDouble(interp, argv[i], &val) != TCL_OK) {
            opserr << "Invalid value for option " << argv[i-1] << "\n";
            return TCL_ERROR;
        }
        if (!gotParam1) {
          gotParam1 = true;
          val1 = val;
          flag1 = static_cast<int>(Isotropy::Parameter::BulkModulus);
        }
        else if (!gotParam2) {
          gotParam2 = true;
          val2 = val;
          flag2 = static_cast<int>(Isotropy::Parameter::BulkModulus);
        }
        else {
            opserr << "Too many elastic parameter options provided.\n";
            return TCL_ERROR;
        }
        positions.insert(i-1);
        positions.insert(i);
      }
      else if (strcmp(argv[i], "-nu") == 0 || strcmp(argv[i], "-poissons-ratio") == 0) {
        if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
        }
        double val;
        if (Tcl_GetDouble(interp, argv[i], &val) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
        }
        if (!gotParam1) { 
          gotParam1 = true; 
          val1 = val; 
          flag1 = static_cast<int>(Isotropy::Parameter::PoissonsRatio);
        }
        else if (!gotParam2) {
          gotParam2 = true; 
          val2 = val; 
          flag2 = static_cast<int>(Isotropy::Parameter::PoissonsRatio);
        }
        else {
          opserr << "Too many elastic parameter options provided.\n";
          return TCL_ERROR;
        }
        positions.insert(i-1);
        positions.insert(i);
      }
      else if (strcmp(argv[i], "-lambda") == 0 || strcmp(argv[i], "-lame-lambda") == 0) {
        if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
        }
        double val;
        if (Tcl_GetDouble(interp, argv[i], &val) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
        }
        if (!gotParam1) {
          gotParam1 = true;
          val1 = val;
          flag1 = static_cast<int>(Isotropy::Parameter::LameLambda);
        }
        else if (!gotParam2) {
          gotParam2 = true;
          val2 = val;
          flag2 = static_cast<int>(Isotropy::Parameter::LameLambda);
        }
        else {
          opserr << "Too many elastic parameter options provided.\n";
          return TCL_ERROR;
        }
        positions.insert(i-1);
        positions.insert(i);
      }
    }
    
    if (!gotParam1 || !gotParam2) {
      opserr << "Must specify exactly two independent elastic parameters.\n";
      return TCL_ERROR;
    }
    
    // Compute canonical Young's modulus and Poisson's ratio.
    int ret = isotropic_constants(flag1, val1, flag2, val2, *iso);
    if (ret != 0) {
        return TCL_ERROR;
    }
    return TCL_OK;
}


template <typename Position, int NDM>
int
TclCommand_newPlasticParser(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);
  // BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  ArgumentTracker<Position> tracker;
  std::set<int> positional;


  int tag;
  double density = 0.0;
  // Isotropy
  IsotropicConstants consts {};
  // Plasticity
  double Fy;
  // Hardening
  double Hiso=0, Hkin=0;


  // Isotropy
  IsotropicParse iso {consts};
  if (TclCommand_setIsotropicParameters((ClientData)&iso, interp, argc, argv) == TCL_OK) {
    tracker.consume(Position::E);
    tracker.consume(Position::G);
    tracker.consume(Position::Nu);
    tracker.consume(Position::K);
    tracker.consume(Position::Lambda);
  }

  // Plasticity
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
    else if (strcmp(argv[i], "-Hiso") == 0 || strcmp(argv[i], "-isotropic-hardening") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &Hiso) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
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
            opserr << OpenSees::PromptParseError << "invalid K.\n";
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
  // Create the material
  //

  return TCL_OK;
}

int
TclCommand_newPlasticMaterial(ClientData clientData, Tcl_Interp *interp,
                            int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  int ndm = builder->getNDM();

  if (ndm==2) {
    enum class Position : int {
      Tag, K, G, YieldStress, EndRequired, 
      E, Nu, Lambda, End
    };
    return TclCommand_newPlasticParser<Position, 2>(clientData, interp, argc, argv);

  } else {
    enum class Position : int {
      Tag, E, Nu, YieldStress, EndRequired, 
      K, G, Lambda, End
    };
    return TclCommand_newPlasticParser<Position, 3>(clientData, interp, argc, argv);
  }
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

  if ((strcmp(argv[1], "J2Plasticity") == 0) || 
      (strcmp(argv[1], "J2") == 0)) {
    builder->addTaggedObject<NDMaterial>(*new J2Plasticity(tag, 0, K, G, sig0, sigInf, delta, H, eta));

    return TCL_OK;
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
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
      opserr << "WARNING invalid sig0\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &sigInf) != TCL_OK) {
      opserr << "WARNING invalid sigInf\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &delta) != TCL_OK) {
      opserr << "WARNING invalid delta\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &H) != TCL_OK) {
      opserr << "WARNING invalid H\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }
    if (argc > 9 && Tcl_GetDouble(interp, argv[9], &eta) != TCL_OK) {
      opserr << "WARNING invalid eta\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial = new J2Plasticity(tag, 0, K, G, sig0, sigInf, delta, H, eta);
  }
}
#endif

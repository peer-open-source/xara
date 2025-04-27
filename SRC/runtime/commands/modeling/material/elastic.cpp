//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// cmp
//
#include <tcl.h>
#include <string.h>
#include <BasicModelBuilder.h>

#include <set>
#include <cstdlib>
#include "BasicModelBuilder.h"
#include "Logging.h"
#include "Parsing.h"
#include "ArgumentTracker.h"
#include "isotropy.h"

#include <ElasticMaterial.h>
#include <ElasticIsotropic.h>
#include <ElasticIsotropicMaterial.h>
#include <ElasticOrthotropicMaterial.h>
#include <material/Solid/ElasticIsotropicThreeDimensional.h>
#include <material/nD/ElasticIsotropicAxiSymm.h>
#include <material/nD/ElasticIsotropicMaterial.h>
#include <material/nD/ElasticIsotropicMaterialThermal.h>
#include <material/nD/ElasticIsotropicPlateFiber.h>

#include <material/nD/ElasticIsotropic3DThermal.h>
#include "ElasticIsotropicPlaneStress2D.h"
#include "ElasticIsotropicPlaneStrain2D.h"
#include <material/elastic/ElasticIsotropicBeamFiber2d.h>
#include <material/elastic/ElasticIsotropicBeamFiber.h>

#include <ElasticIsotropicMaterial.h>
// #include <ElasticCrossAnisotropic.h>
#include <PlaneStressMaterial.h>
#include <PlateFiberMaterial.h>
#include <BeamFiberMaterial.h>



template <typename Position>
int
TclCommand_newElasticParser(ClientData clientData, Tcl_Interp *interp,
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
  // Keywords
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
      
      case Position::Eta:
        if (Tcl_GetDouble (interp, argv[i], &eta) != TCL_OK) {
            opserr << OpenSees::PromptParseError << "invalid eta.\n";
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
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  if ((strcmp(argv[1], "ElasticIsotropic") == 0) ||
      (strcmp(argv[1], "Elastic") == 0)) {
    double E = consts.E;
    double nu = consts.nu;
    if (builder->addTaggedObject<NDMaterial>(*new ElasticIsotropicMaterial(tag, E, nu, density)) != TCL_OK ) {
      return TCL_ERROR;
    }
    if (strcmp(argv[0], "material") == 0) {
      if (builder->addTaggedObject<UniaxialMaterial>(*new ElasticMaterial(tag, E, eta, E)) != TCL_OK ) {
        return TCL_ERROR;
      }
      if (builder->addTaggedObject<Mate<3>>(*new ElasticIsotropic<3>(tag, E, nu, density)) != TCL_OK ) {
        return TCL_ERROR;
      }
    }
    return TCL_OK;
  }
  else if (strcmp(argv[1], "ElasticIsotropic3D") == 0) {
    double E = consts.E;
    double nu = consts.nu;
    if (builder->addTaggedObject<NDMaterial>(*new ElasticIsotropicThreeDimensional(tag, E, nu, density)) != TCL_OK ) {
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  return TCL_ERROR;
}

int
TclCommand_newElasticMaterial(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  // 
  if (strcmp(argv[1], "ElasticIsotropic") == 0 ||
      strcmp(argv[1], "PlaneStressSimplifiedJ2") == 0) {

    // "ElasticIsotropic" tag?  E?  nu?  rho?
    enum class Position : int {
      Tag, E, Nu,  EndRequired, 
      Density, End,
      G, K, Lambda, Eta
    };
    return TclCommand_newElasticParser<Position>(clientData, interp, argc, argv);
  }

  return TCL_ERROR;
}



enum class MaterialSymmetry {
  Triclinic     , // 21
  Monoclinic    , // 13
  Orthorhombic  , //  9
  Tetragonal    , //  7
//Tetragonal    , //  6
  Rhombohedral  , //  7
//Rhombohedral  , //  6
  Hexagonal     , //  5
  Cubic         , //  3
  Isotropic       //  2
};


int
TclCommand_newElasticAnisotropic(ClientData clientData, Tcl_Interp* interp, int argc, const char**const argv)
{
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  MaterialSymmetry symm = MaterialSymmetry::Isotropic;
  PlaneType type = PlaneType::None;

  if ((strcmp(argv[1], "ElasticIsotropic") == 0) || 
      (strcmp(argv[1], "Elastic") == 0) ||
      (strcmp(argv[1], "ElasticBeamFiber") == 0) ||
      (strcmp(argv[1], "ElasticIsotropic3D") == 0))
  {

    int tag;
    double E, v;
    double rho = 0.0;

    int loc = 2;
    if (Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
      opserr << "WARNING invalid ElasticIsotropic tag" << "\n";
      return TCL_ERROR;
    }
    loc++;

    if (Tcl_GetDouble(interp, argv[loc], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      return TCL_ERROR;
    }
    loc++;

    if (Tcl_GetDouble(interp, argv[loc], &v) != TCL_OK) {
      opserr << "WARNING invalid v\n";
      return TCL_ERROR;
    }
    loc++;

    if (argc > loc && Tcl_GetDouble(interp, argv[loc], &rho) != TCL_OK) {
      opserr << "WARNING invalid rho\n";
      return TCL_ERROR;
    }
    loc++;


    builder->addTaggedObject<UniaxialMaterial>(*new ElasticMaterial(tag, E, 0.0, E));
    builder->addTaggedObject<NDMaterial>(*new ElasticIsotropicMaterial(tag, E, v, rho));
    builder->addTaggedObject<Mate<3>>   (*new ElasticIsotropic<3>(tag, E, v, rho));

    // builder->addTaggedObject<NDMaterial,"PlaneFrame" >(*new ElasticIsotropicBeamFiber(tag, E, v, rho));
    // builder->addTaggedObject<NDMaterial,"PlaneStrain">(*new ElasticIsotropicPlaneStrain2D(tag, E, v, rho));
    // builder->addTaggedObject<Mate<2>,   "PlaneStrain">(*new ElasticIsotropic<2,PlaneType::Strain>(tag, E, v, rho));
    // builder->addTaggedObject<NDMaterial,"PlaneStress">(*new ElasticIsotropicPlaneStress2D(tag, E, v, rho));
    // builder->addTaggedObject<Mate<2>,   "PlaneStress">(*new ElasticIsotropic<2,PlaneType::Stress>(tag, E, v, rho));

    return TCL_OK;
  }

  else if (strcmp(argv[1], "ElasticCrossAnisotropic") == 0)
  {
    if (argc < 8) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial ElasticCrossAnisotropic tag? Ehh? Ehv? nuhv? nuvv? Ghv? <rho?>"
             << "\n";
      return TCL_ERROR;
    }

    int tag;
    double Eh, Ev, nuhv, nuhh, Ghv;
    double rho = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid ElasticCrossAnisotropic tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &Eh) != TCL_OK) {
      opserr << "WARNING invalid Eh\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &Ev) != TCL_OK) {
      opserr << "WARNING invalid Ev\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &nuhv) != TCL_OK) {
      opserr << "WARNING invalid nuhv\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &nuhh) != TCL_OK) {
      opserr << "WARNING invalid nuhh\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &Ghv) != TCL_OK) {
      opserr << "WARNING invalid Ghv\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    if (argc > 8 && Tcl_GetDouble(interp, argv[8], &rho) != TCL_OK) {
      opserr << "WARNING invalid rho\n";
      opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << "\n";
      return TCL_ERROR;
    }

    // theMaterial = new ElasticCrossAnisotropic(tag, Eh, Ev, nuhv, nuhh, Ghv, rho);
  }
}


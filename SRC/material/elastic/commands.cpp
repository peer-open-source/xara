//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <string.h>
#include <BasicModelBuilder.h>

#include <set>
#include <cctype>      // for isdigit()
#include <cstdlib>
#include "BasicModelBuilder.h"
#include "Logging.h"
#include "Parsing.h"
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

#if 1
// TclCommand_newIsotropicMaterial
//   Default (positional) syntax:
//      material Isotropic $tag $young $poisson <$rho>
//   Option syntax:
//      material Isotropic $tag <options>
//   Valid elastic parameter options are:
//      -E (or -youngs-modulus), -G (or -shear-modulus),
//      -K (or -bulk-modulus), -nu (or -poissons-ratio),
//      -lambda (or -lame-lambda)
//   Optionally, -rho (or -density) specifies density.
int
TclCommand_newIsotropicMaterial(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
    // Retrieve the model builder.
    BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

    // Determine which form is being used.
    // If none of the arguments (starting at argv[2]) begins with '-' (except when followed by a digit)
    // we assume the default positional form.
    bool defaultForm = true;
    for (int i = 2; i < argc; i++) {
        if (argv[i][0] == '-' && argv[i][1] != '\0' && !isdigit(argv[i][1])) {
            defaultForm = false;
            break;
        }
    }

    int tag;
    double E=0.0, nu=0.0, rho = 0.0;
    bool useRho = false;

    if (defaultForm) {
        // Default positional form: material Isotropic $tag $young $poisson <$rho>
        if (argc < 5 || argc > 6) {
            opserr << "Invalid number of arguments for default isotropic material.\n";
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
            opserr << "Invalid tag: " << argv[2] << "\n";
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
            opserr << "Invalid Young's modulus: " << argv[3] << "\n";
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[4], &nu) != TCL_OK) {
            opserr << "Invalid Poisson's ratio: " << argv[4] << "\n";
            return TCL_ERROR;
        }
        useRho = (argc == 6);
        if (useRho) {
            if (Tcl_GetDouble(interp, argv[5], &rho) != TCL_OK) {
                opserr << "Invalid density: " << argv[5] << "\n";
                return TCL_ERROR;
            }
        }
    }
    else {
        // Option-based form.
        // We expect the syntax:
        //    material Isotropic $tag <options>
        // The tag is still the first non-option argument (argv[2]).
        // Among the options, exactly two independent elastic parameters must be provided.
        // Valid elastic options: -E, -G, -K, -nu, -lambda.
        // An optional option: -rho (or -density) for density.
        bool gotParam1 = false, gotParam2 = false;
        int flag1 = 0, flag2 = 0;
        double val1 = 0.0, val2 = 0.0;

        // Parse tag from argv[2].
        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
            opserr << "Invalid tag: " << argv[2] << "\n";
            return TCL_ERROR;
        }
        // Process the remaining arguments.
        for (int i = 3; i < argc; i++) {
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
                    gotParam2 = true; val2 = val; flag2 = static_cast<int>(Isotropy::Parameter::YoungModulus); }
                else {
                    opserr << "Too many elastic parameter options provided.\n";
                    return TCL_ERROR;
                }
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
                if (!gotParam1) { gotParam1 = true; val1 = val; flag1 = static_cast<int>(Isotropy::Parameter::ShearModulus); }
                else if (!gotParam2) { gotParam2 = true; val2 = val; flag2 = static_cast<int>(Isotropy::Parameter::ShearModulus); }
                else {
                    opserr << "Too many elastic parameter options provided.\n";
                    return TCL_ERROR;
                }
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
                if (!gotParam1) { gotParam1 = true; val1 = val; flag1 = static_cast<int>(Isotropy::Parameter::BulkModulus); }
                else if (!gotParam2) { gotParam2 = true; val2 = val; flag2 = static_cast<int>(Isotropy::Parameter::BulkModulus); }
                else {
                    opserr << "Too many elastic parameter options provided.\n";
                    return TCL_ERROR;
                }
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
                if (!gotParam1) { gotParam1 = true; val1 = val; flag1 = static_cast<int>(Isotropy::Parameter::LameLambda); }
                else if (!gotParam2) { gotParam2 = true; val2 = val; flag2 = static_cast<int>(Isotropy::Parameter::LameLambda); }
                else {
                    opserr << "Too many elastic parameter options provided.\n";
                    return TCL_ERROR;
                }
            }
            else if (strcmp(argv[i], "-rho") == 0 || strcmp(argv[i], "-density") == 0) {
                if (++i >= argc) {
                    opserr << "Missing value for option " << argv[i-1] << "\n";
                    return TCL_ERROR;
                }
                if (Tcl_GetDouble(interp, argv[i], &rho) != TCL_OK) {
                    opserr << "Invalid density value for option " << argv[i-1] << "\n";
                    return TCL_ERROR;
                }
                useRho = true;
            }
            else {
                opserr << "Unknown option: " << argv[i] << "\n";
                return TCL_ERROR;
            }
        }
        
        if (!gotParam1 || !gotParam2) {
            opserr << "Must specify exactly two independent elastic parameters.\n";
            return TCL_ERROR;
        }
        
        // Use the conversion utility to compute canonical Young's modulus and Poisson's ratio.
        int ret = isotropic_convert(flag1, val1, flag2, val2,
                                       static_cast<int>(Isotropy::Parameter::YoungModulus), E);
        if (ret != 0) {
            opserr << "Error converting elastic parameters to Young's modulus.\n";
            return TCL_ERROR;
        }
        ret = isotropic_convert(flag1, val1, flag2, val2,
                                   static_cast<int>(Isotropy::Parameter::PoissonsRatio), nu);
        if (ret != 0) {
            opserr << "Error converting elastic parameters to Poisson's ratio.\n";
            return TCL_ERROR;
        }
    }
    
    // Create the isotropic material.
    NDMaterial *theMaterial = nullptr;
    if (useRho)
        theMaterial = new ElasticIsotropicMaterial(tag, E, nu, rho);
    else
        theMaterial = new ElasticIsotropicMaterial(tag, E, nu);

    if (theMaterial == nullptr || builder->addTaggedObject<NDMaterial>(*theMaterial) < 0) {
        if (theMaterial != nullptr)
            delete theMaterial;
        return TCL_ERROR;
    }
    return TCL_OK;
}
#endif


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
TclCommand_newElasticMaterial(ClientData clientData, Tcl_Interp* interp, int argc, const char**const argv)
{
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  MaterialSymmetry symm = MaterialSymmetry::Isotropic;
  PlaneType type = PlaneType::None;

  if ((strcmp(argv[1], "ElasticIsotropic") == 0) || 
      (strcmp(argv[1], "Elastic") == 0) ||
      (strcmp(argv[1], "ElasticBeamFiber") == 0) ||
      (strcmp(argv[1], "ElasticIsotropic3D") == 0))

  {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial ElasticIsotropic tag? E? v? <rho?>" << "\n";
      return TCL_ERROR;
    }

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


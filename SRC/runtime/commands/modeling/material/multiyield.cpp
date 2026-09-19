//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function invoked when the user 
// invokes the nDMaterial command in the interpreter.
//
// Developed by:
//   Frank McKenna (fmckenna@ce.berkeley.edu)
//   Gregory L. Fenves (fenves@ce.berkeley.edu)
//   Filip C. Filippou (filippou@ce.berkeley.edu)
//
// With extensive contributions from
//   Boris Jeremic    (jeremic@ucdavis.edu)
//   Zaohui Yang      (zhyang@ucdavis.edu)
//   Zhao Cheng       (zcheng@ucdavis.edu)
//
//===----------------------------------------------------------------------===//
//
#include <Parsing.h>
#include <Logging.h>
#include <ModelRegistry.h>

#include <MultiYieldSurfaceClay.h>
#include <FluidSolidPorousMaterial.h>
#include <PressureDependentElastic3D.h>
#include <PressureDependMultiYield.h>
#include <PressureDependMultiYield02.h>
#include <PressureDependMultiYield03.h>
#include <PressureIndependMultiYield.h>


int
XaraMatCmd_MultiYieldSurface(ClientData clientData,
                             Tcl_Interp *interp,
                             ArgSize argc,
                             TCL_Char ** const argv)
{
  NDMaterial *theMaterial = nullptr;
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  // Pressure Independent Multi-yield, by ZHY
  if (strcmp(argv[1], "PressureIndependMultiYield") == 0) {
    const int numParam = 6;
    const int totParam = 10;
    int tag;
    double param[totParam];
    param[6] = 0.0;
    param[7] = 100.;
    param[8] = 0.0;
    param[9] = 20;

    const char *arg[] = {"nd",
                         "rho",
                         "refShearModul",
                         "refBulkModul",
                         "cohesi",
                         "peakShearStra",
                         "frictionAng (=0)",
                         "refPress (=100)",
                         "pressDependCoe (=0.0)",
                         "numberOfYieldSurf (=20)"};
    if (argc < (3 + numParam)) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PressureIndependMultiYield tag? " << arg[0];
      opserr << "? "
             << "\n";
      opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? "
             << "\n";
      opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? "
             << "\n";
      opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? " << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureIndependMultiYield tag" << "\n";
      return TCL_ERROR;
    }

    for (int i = 3; (i < argc && i < 13); ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        return TCL_ERROR;
      }

    static double *gredu = 0;
    // user defined yield surfaces
    if (param[9] < 0 && param[9] > -40) {
      param[9] = -int(param[9]);
      gredu = new double[int(2 * param[9])];
      for (int i = 0; i < 2 * param[9]; ++i)
        if (Tcl_GetDouble(interp, argv[i + 13], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          return TCL_ERROR;
        }
    }

    PressureIndependMultiYield *temp = new PressureIndependMultiYield(
        tag, param[0], param[1], param[2], param[3], param[4], param[5],
        param[6], param[7], param[8], param[9], gredu);
    theMaterial = temp;

    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }

  // Pressure Independent Multi-yield, by Quan Gu
  else if (strcmp(argv[1], "MultiYieldSurfaceClay") == 0) {
    const int numParam = 6;
    const int totParam = 10;
    int tag;
    double param[totParam];
    param[6] = 0.0;
    param[7] = 100.;
    param[8] = 0.0;
    param[9] = 20;

    const char *arg[] = {"nd",
                         "rho",
                         "refShearModul",
                         "refBulkModul",
                         "cohesi",
                         "peakShearStra",
                         "frictionAng (=0)",
                         "refPress (=100)",
                         "pressDependCoe (=0.0)",
                         "numberOfYieldSurf (=20)"};
    if (argc < (3 + numParam)) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial MultiYieldSurfaceClay tag? " << arg[0];
      opserr << "? "
             << "\n";
      opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? "
             << "\n";
      opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? "
             << "\n";
      opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? " << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid MultiYieldSurfaceClay tag" << "\n";
      return TCL_ERROR;
    }

    for (int i = 3; (i < argc && i < 13); ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        opserr << "nDMaterial MultiYieldSurfaceClay: " << tag << "\n";
        return TCL_ERROR;
      }

    static double *gredu = 0;
    // user defined yield surfaces
    if (param[9] < 0 && param[9] > -40) {
      param[9] = -int(param[9]);
      gredu = new double[int(2 * param[9])];
      for (int i = 0; i < 2 * param[9]; ++i)
        if (Tcl_GetDouble(interp, argv[i + 13], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          return TCL_ERROR;
        }
    }

    MultiYieldSurfaceClay *temp = new MultiYieldSurfaceClay(
        tag, param[0], param[1], param[2], param[3], param[4], param[5],
        param[6], param[7], param[8], param[9], gredu);
    theMaterial = temp;

    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }
  // ============

  else if (strcmp(argv[1], "PressureDependentElastic3D") == 0) {
    if (argc < 6) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PressureDependentElastic3D tag? E? v? rho?"
             << "\n";
      return TCL_ERROR;
    }

    int  tag = 0;
    double E = 0.0;
    double v = 0.0;
    double rho = 0.0;
    double expp = 0.0;
    double prp = 0.0;
    double pop = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureDependentElastic3D tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &v) != TCL_OK) {
      opserr << "WARNING invalid v\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &rho) != TCL_OK) {
      opserr << "WARNING invalid rho\n";
      return TCL_ERROR;
    }

    if (argc == 6) {
      theMaterial = new PressureDependentElastic3D(tag, E, v, rho);
    }

    else if (argc == 7) {
      // get the exponent of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
        opserr << "WARNING invalid expp\n";
        return TCL_ERROR;
      }
      theMaterial = new PressureDependentElastic3D(tag, E, v, rho, expp);
    }

    else if (argc == 8) {
      // get the exponent pressure of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
        opserr << "WARNING invalid expp\n";
        return TCL_ERROR;
      }
      // get the reference pressure of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[7], &prp) != TCL_OK) {
        opserr << "WARNING invalid prp\n";
        return TCL_ERROR;
      }
      theMaterial = new PressureDependentElastic3D(tag, E, v, rho, expp, prp);
    }

    else if (argc >= 9) {
      // get the exponent of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
        opserr << "WARNING invalid expp\n";
        return TCL_ERROR;
      }
      // get the reference pressure of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[7], &prp) != TCL_OK) {
        opserr << "WARNING invalid prp\n";
        return TCL_ERROR;
      }
      // get the cutoff pressure po of the pressure sensitive elastic material)
      if (Tcl_GetDouble(interp, argv[8], &pop) != TCL_OK) {
        opserr << "WARNING invalid pop\n";
        return TCL_ERROR;
      }
      theMaterial =
          new PressureDependentElastic3D(tag, E, v, rho, expp, prp, pop);
    }
  }

  // Pressure Dependent Multi-yield, by ZHY
  else if (strcmp(argv[1], "PressureDependMultiYield") == 0) {
    const int numParam = 15;
    const int totParam = 24;
    int tag;
    double param[totParam];
    param[15] = 20;
    param[16] = 0.6;
    param[17] = 0.9;
    param[18] = 0.02;
    param[19] = 0.7;
    param[20] = 101.;
    param[21] = .3;
    param[22] = 0.;
    param[23] = 1.;

    const char *arg[] = {"nd",
                         "rho",
                         "refShearModul",
                         "refBulkModul",
                         "frictionAng",
                         "peakShearStra",
                         "refPress",
                         "pressDependCoe",
                         "phaseTransformAngle",
                         "contractionParam1",
                         "dilationParam1",
                         "dilationParam2",
                         "liquefactionParam1",
                         "liquefactionParam2",
                         "liquefactionParam4",
                         "numberOfYieldSurf (=20)",
                         "e (=0.6)",
                         "volLimit1 (=0.9)",
                         "volLimit2 (=0.02)",
                         "volLimit3 (=0.7)",
                         "Atmospheric pressure (=101)",
                         "cohesi (=.5)",
                         "Hv (=0)",
                         "Pv (=1.)"};
    if (argc < (3 + numParam)) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PressureDependMultiYield tag? " << arg[0];
      opserr << "? "
             << "\n";
      opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? "
             << "\n";
      opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? "
             << "\n";
      opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? "
             << "\n";
      opserr << arg[10] << "? " << arg[11] << "? " << arg[12] << "? "
             << "\n";
      opserr << arg[13] << "? " << arg[14] << "? " << arg[15] << "? "
             << "\n";
      opserr << arg[16] << "? " << arg[17] << "? " << arg[18] << "? "
             << "\n";
      opserr << arg[19] << "? " << arg[20] << "? " << arg[21] << "? " << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureDependMultiYield tag" << "\n";
      return TCL_ERROR;
    }

    for (int i = 3; (i < argc && i < 19); ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        return TCL_ERROR;
      }

    static double *gredu = 0;
    // user defined yield surfaces
    if (param[15] < 0 && param[15] > -40) {
      param[15] = -int(param[15]);
      gredu = new double[int(2 * param[15])];

      for (int i = 0; i < 2 * param[15]; ++i)
        if (Tcl_GetDouble(interp, argv[i + 19], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          return TCL_ERROR;
        }
    }

    if (gredu != 0) {
      for (int i = 19 + int(2 * param[15]); i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i],
                          &param[i - 3 - int(2 * param[15])]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[15])]
                 << "\n";
          return TCL_ERROR;
        }
    } else {
      for (int i = 19; i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[15])]
                 << "\n";
          return TCL_ERROR;
        }
    }

    PressureDependMultiYield *temp = new PressureDependMultiYield(
        tag, param[0], param[1], param[2], param[3], param[4], param[5],
        param[6], param[7], param[8], param[9], param[10], param[11], param[12],
        param[13], param[14], param[15], gredu, param[16], param[17], param[18],
        param[19], param[20], param[21], param[22], param[23]);

    theMaterial = temp;
    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }

  // Pressure Dependent Multi-yield, by ZHY
  else if (strcmp(argv[1], "PressureDependMultiYield02") == 0) {
    static constexpr int numParam = 13;
    static constexpr int totParam = 26;

    int tag;
    double param[totParam];
    param[numParam] = 20;
    param[numParam + 1] = 5.0;
    param[numParam + 2] = 3.;
    param[numParam + 3] = 1.;
    param[numParam + 4] = 0.;
    param[numParam + 5] = 0.6;
    param[numParam + 6] = 0.9;
    param[numParam + 7] = 0.02;
    param[numParam + 8] = 0.7;
    param[numParam + 9] = 101.;
    param[numParam +10] = 0.1;
    param[numParam +11] = 0.;
    param[numParam +12] = 1.;

    const char *arg[] = {"nd",
                         "rho",
                         "refShearModul",
                         "refBulkModul",
                         "frictionAng",
                         "peakShearStra",
                         "refPress",
                         "pressDependCoe",
                         "phaseTransformAngle",
                         "contractionParam1",
                         "contractionParam3",
                         "dilationParam1",
                         "dilationParam3",
                         "numberOfYieldSurf (=20)",
                         "contractionParam2=5.0",
                         "dilationParam2=3.0",
                         "liquefactionParam1=1.0",
                         "liquefactionParam2=0.0",
                         "e (=0.6)",
                         "volLimit1 (=0.9)",
                         "volLimit2 (=0.02)",
                         "volLimit3 (=0.7)",
                         "Atmospheric pressure (=101)",
                         "cohesi (=.1)",
                         "Hv (=0)",
                         "Pv (=1.)"};
    if (argc < (3 + numParam)) {
      opserr << "WARNING insufficient arguments\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureDependMultiYield02 tag" << "\n";
      return TCL_ERROR;
    }

    int in = 17;
    for (int i = 3; (i < argc && i < in); ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        return TCL_ERROR;
      }

    static double *gredu = 0;

    // user defined yield surfaces
    if (param[numParam] < 0 && param[numParam] > -100) {
      param[numParam] = -int(param[numParam]);
      gredu = new double[int(2 * param[numParam])];

      for (int i = 0; i < 2 * param[numParam]; ++i)
        if (Tcl_GetDouble(interp, argv[i + in], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          return TCL_ERROR;
        }
    }

    if (gredu != 0) {
      for (int i = in + int(2 * param[numParam]); i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i],
                          &param[i - 3 - int(2 * param[numParam])]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])]
                 << "\n";
          return TCL_ERROR;
        }
    } else {
      for (int i = in; i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])]
                 << "\n";
          return TCL_ERROR;
        }
    }

    PressureDependMultiYield02 *temp = new PressureDependMultiYield02(
        tag, param[0], param[1], param[2], param[3], param[4], param[5],
        param[6], param[7], param[8], param[9], param[10], param[11], param[12],
        param[13], gredu, param[14], param[15], param[16], param[17], param[18],
        param[19], param[20], param[21], param[22], param[23], param[24],
        param[25]);

    theMaterial = temp;
    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }

  // nDMaterial PressureDependMultiYield03  $tag  $nd  $rho  $refShearModul
  // $refBulkModul $frictionAng  $peakShearStra  $refPress  $pressDependCoe
  // $PTAng $mType $ca  $cb $cc $cd $ce $da $db $dc <$noYieldSurf=20
  // <$r1 $Gs1 …>  $liquefac1=1. $liquefac2=0. $pa=101 <$c=1.73>>
  // PressureDependMultiYield03 (based on PressureDependMultiYield02).
  else if (strcmp(argv[1], "PressureDependMultiYield03") == 0) {
    static constexpr int numParam = 18;
    static constexpr int totParam = 23;
    int tag;
    double param[totParam];
    param[numParam] = 20;
    param[numParam + 1] = 1.;
    param[numParam + 2] = 0.;
    param[numParam + 3] = 101.;
    param[numParam + 4] = 1.73;

    const char *arg[] = {"nd",
                         "rho",
                         "refShearModul",
                         "refBulkModul",
                         "frictionAng",
                         "peakShearStra",
                         "refPress",
                         "pressDependCoe",
                         "phaseTransformAngle",
                         "mType",
                         "ca",
                         "cb",
                         "cc",
                         "cd",
                         "ce",
                         "da",
                         "db",
                         "dc",
                         "numberOfYieldSurf (=20)",
                         "liquefactionParam1=1.0",
                         "liquefactionParam2=0.0",
                         "Atmospheric pressure (=101)",
                         "cohesi (=1.73)"};

    if (argc < (3 + numParam)) { // 3 refers to "nDMaterial
                                 // PressureDependMultiYield03  $tag"
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PressureDependMultiYield03 tag? " << arg[0];
      opserr << "? "
             << "\n";
      opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? "
             << "\n";
      opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? "
             << "\n";
      opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? "
             << "\n";
      opserr << arg[10] << "? " << arg[11] << "? " << arg[12] << "? "
             << "\n";
      opserr << arg[13] << "? " << arg[14] << "? " << arg[15] << "? "
             << "\n";
      opserr << arg[16] << "? " << arg[17] << "? " << arg[18] << "? "
             << "\n";
      opserr << arg[19] << "? " << arg[20] << "? " << arg[21] << "? " << arg[22]
             << "? " << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid PressureDependMultiYield03 tag" << "\n";
      return TCL_ERROR;
    }

    int in = 22;
    for (int i = 3; (i < argc && i < in); ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        return TCL_ERROR;
      }

    static double *gredu = 0;

    // user defined yield surfaces
    if (param[numParam] < 0 && param[numParam] > -100) {
      param[numParam] = -int(param[numParam]);
      gredu = new double[int(2 * param[numParam])];

      for (int i = 0; i < 2 * param[numParam]; ++i)
        if (Tcl_GetDouble(interp, argv[i + in], &gredu[i]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3] << "\n";
          return TCL_ERROR;
        }
    }

    if (gredu != 0) {
      for (int i = in + int(2 * param[numParam]); i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i],
                          &param[i - 3 - int(2 * param[numParam])]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])]
                 << "\n";
          opserr << "nDMaterial PressureDependMultiYield03: " << tag << "\n";
          return TCL_ERROR;
        }
    } else {
      for (int i = in; i < argc; ++i)
        if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
          opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])]
                 << "\n";
          opserr << "nDMaterial PressureDependMultiYield03: " << tag << "\n";
          return TCL_ERROR;
        }
    }

    PressureDependMultiYield03 *temp = new PressureDependMultiYield03(
        tag, param[0], param[1], param[2], param[3], param[4], param[5],
        param[6], param[7], param[8], param[9], param[10], param[11], param[12],
        param[13], param[14], param[15], param[16], param[17], param[18], gredu,
        param[19], param[20], param[21], param[22]);

    theMaterial = temp;
    if (gredu != 0) {
      delete[] gredu;
      gredu = 0;
    }
  }

  // Fluid Solid Porous, by ZHY
  else if (strcmp(argv[1], "FluidSolidPorous") == 0) {

    int tag;
    double param[4];
    const char *arg[] = {"nd", "soilMatTag", "combinedBulkModul",
                   "Atmospheric pressure"};
    if (argc < 6) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial FluidSolidPorous tag? " << arg[0];
      opserr << "? "
             << "\n";
      opserr << arg[1] << "? " << arg[2] << "? " << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid FluidSolidPorous tag" << "\n";
      return TCL_ERROR;
    }

    for (int i = 3; i < 6; ++i)
      if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[i - 3] << "\n";
        return TCL_ERROR;
      }

    NDMaterial *soil = builder->getTypedObject<NDMaterial>(param[1]);
    if (soil == nullptr)
      return TCL_ERROR;


    param[3] = 101.;
    if (argc == 7) {
      if (Tcl_GetDouble(interp, argv[6], &param[3]) != TCL_OK) {
        opserr << "WARNING invalid " << arg[3] << "\n";
        return TCL_ERROR;
      }
    }

    theMaterial =
        new FluidSolidPorousMaterial(tag, param[0], *soil, param[2], param[3]);
  }


  //
  //
  //
  if (theMaterial == nullptr) {
    opserr << "WARNING could not create material\n";
    return TCL_ERROR;
  }

  if (theMaterial == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "could not create nDMaterial " << argv[1]
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
    opserr << OpenSees::PromptValueError << "could not add material to the domain\n";
    opserr << *theMaterial << "\n";
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;
}

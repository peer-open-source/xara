//
// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the parsing routines for the
// TCL hystereticBackbone command.
//
#include <assert.h>
#include <OPS_Globals.h>
#include <Parsing.h>
#include <UniaxialMaterial.h>
#include <BasicModelBuilder.h>

#include <runtimeAPI.h>

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       G3_Runtime *rt, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

#include <TrilinearBackbone.h>
#include <MultilinearBackbone.h>
#include <ArctangentBackbone.h>
#include <CappedBackbone.h>
#include <MaterialBackbone.h>
#include <LinearCappedBackbone.h>
#include <ReeseSoftClayBackbone.h>
#include <ReeseStiffClayBelowWS.h>
#include <ReeseSandBackbone.h>
#include <ManderBackbone.h>
#include <RaynorBackbone.h>
#include <string.h>

extern OPS_Routine OPS_ArctangentBackbone;
// extern OPS_Routine OPS_ManderBackbone;
extern OPS_Routine OPS_TrilinearBackbone;
HystereticBackbone* TclCommand_newManderBackbone(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char ** const argv);
extern OPS_Routine OPS_BilinearBackbone;
extern OPS_Routine OPS_MultilinearBackbone;

int
TclCommand_addHystereticBackbone(ClientData clientData,
                                 Tcl_Interp *interp, int argc,
                                 TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  G3_Runtime *rt = G3_getRuntime(interp);

  if (argc < 3) {
    opserr << "WARNING insufficient number of hystereticBackbone arguments\n";
    opserr << "Want: hystereticBackbone type? tag? <specific hystereticBackbone args>"
           << endln;
    return TCL_ERROR;
  }



  // Pointer to a hysteretic backbone that will be added to the model builder
  HystereticBackbone *theBackbone = nullptr;

  // Check argv[1] for backbone type
  if (strcmp(argv[1], "Bilinear") == 0) {
    void *theBB = OPS_BilinearBackbone(rt, argc, argv);
    if (theBB != nullptr)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Trilinear") == 0) {
    void *theBB = OPS_TrilinearBackbone(rt, argc, argv);
    if (theBB != nullptr)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Multilinear") == 0) {
    void *theBB = OPS_MultilinearBackbone(rt, argc, argv);
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Arctangent") == 0) {
    void *theBB = OPS_ArctangentBackbone(rt, argc, argv);
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "ReeseSoftClay") == 0) {
    if (argc < 6) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: hystereticBackbone ReeseSoftClay tag? pu? y50? n?"
             << endln;
      return TCL_ERROR;
    }

    int tag;
    double pu, y50, n;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSoftClay tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &pu) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSoftClay pu" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &y50) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSoftClay y50" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &n) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSoftClay n" << endln;
      return TCL_ERROR;
    }

    theBackbone = new ReeseSoftClayBackbone(tag, pu, y50, n);
  }

  else if (strcmp(argv[1], "ReeseSand") == 0) {
    if (argc < 8) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: hystereticBackbone ReeseSand tag? kx? ym? pm? yu? pu?"
             << endln;
      return TCL_ERROR;
    }

    int tag;
    double kx, ym, pm, yu, pu;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &kx) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand kx" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &ym) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand ym" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &pm) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand pm" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6], &yu) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand yu" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[7], &pu) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand pu" << endln;
      return TCL_ERROR;
    }

    theBackbone = new ReeseSandBackbone(tag, kx, ym, pm, yu, pu);
  }

  else if (strcmp(argv[1], "ReeseStiffClayBelowWS") == 0) {
    if (argc < 7) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: hystereticBackbone ReeseStiffClayBelowWS tag? Esi? y50? "
                "As?  Pc?"
             << endln;
      return TCL_ERROR;
    }

    int tag;
    double Esi, y50, As, Pc;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseStiffClayBelowWS tag"
             << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &Esi) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseStiffClayBelowWS Esi"
             << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &y50) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseStiffClayBelowWS y50"
             << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &As) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseStiffClayBelowWS As"
             << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6], &Pc) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseStiffClayBelowWS Pc"
             << endln;
      return TCL_ERROR;
    }

    theBackbone = new ReeseStiffClayBelowWS(tag, Esi, y50, As, Pc);

  }

  else if (strcmp(argv[1], "Mander") == 0) {
    // void *theBB = OPS_ManderBackbone(rt, argc, argv);
    void *theBB = TclCommand_newManderBackbone(clientData, interp, argc, argv);

    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Raynor") == 0) {
    if (argc < 10) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: hystereticBackbone Raynor tag? Es? fy? fsu? Epsilonsh? "
                "Epsilonsm? C1? Ey?"
             << endln;
      return TCL_ERROR;
    }

    int tag;
    double Es, fy, fsu, Epsilonsh, Epsilonsm, C1, Ey;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Raynor tag" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &Es) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Raynor Es" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &fy) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Raynor fy" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &fsu) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Raynor fsu" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6], &Epsilonsh) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Raynor Epsilonsh" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[7], &Epsilonsm) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Raynor Epsilonsm" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &C1) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Raynor fy" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[9], &Ey) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Raynor fsu" << endln;
      return TCL_ERROR;
    }

    theBackbone =
        new RaynorBackbone(tag, Es, fy, fsu, Epsilonsh, Epsilonsm, C1, Ey);

  }

  else if (strcmp(argv[1], "Capped") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: hystereticBackbone Capped tag? hystereticBackboneTag? "
                "capTag?"
             << endln;
      return TCL_ERROR;
    }

    int tag, bTag, cTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Capped tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &bTag) != TCL_OK) {
      opserr
          << "WARNING invalid hystereticBackbone Capped hystereticBackboneTag"
          << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &cTag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Capped capTag" << endln;
      return TCL_ERROR;
    }

    HystereticBackbone *hystereticBackbone = 0;

    // hystereticBackbone = theTclBuilder->getHystereticBackbone(bTag);

    if (hystereticBackbone == 0) {
      opserr << "WARNING hystereticBackbone does not exist\n";
      opserr << "hystereticBackbone: " << bTag;
      opserr << "\nhystereticBackbone Capped: " << tag << endln;
      return TCL_ERROR;
    }

    HystereticBackbone *cap = 0;

    // cap = theTclBuilder->getHystereticBackbone(cTag);

    if (cap == 0) {
      opserr << "WARNING cap does not exist\n";
      opserr << "cap: " << cTag;
      opserr << "\nhystereticBackbone Capped: " << tag << endln;
      return TCL_ERROR;
    }

    // theBackbone = new CappedBackbone (tag, *hystereticBackbone, *cap);
    theBackbone = 0;
  }

  else if (strcmp(argv[1], "LinearCapped") == 0) {
    if (argc < 7) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: hystereticBackbone LinearCapped tag? backboneTag? eCap? "
                "E? sRes?"
             << endln;
      return TCL_ERROR;
    }

    int tag, bTag;
    double e, E, s;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone LinearCapped tag" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &bTag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone LinearCapped backboneTag"
             << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &e) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone LinearCapped eCap" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &E) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone LinearCapped E" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &s) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone LinearCapped sRes" << endln;
      return TCL_ERROR;
    }

    HystereticBackbone *backbone = 0;
    // backbone = theTclBuilder->getHystereticBackbone(bTag);

    if (backbone == 0) {
      opserr << "WARNING hystereticBackbone does not exist\n";
      opserr << "hystereticBackbone: " << bTag;
      opserr << "\nhystereticBackbone LinearCapped: " << tag << endln;
      return TCL_ERROR;
    }

    // theBackbone = new LinearCappedBackbone (tag, *backbone, e, E, s);
    theBackbone = 0;
  }

  else if (strcmp(argv[1], "Material") == 0) {
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: hystereticBackbone Material tag? matTag?" << endln;
      return TCL_ERROR;
    }

    int tag, matTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid tag\n";
      opserr << "hystereticBackbone Material: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag\n";
      opserr << "hystereticBackbone Material: " << tag << endln;
      return TCL_ERROR;
    }

    UniaxialMaterial *material = builder->getTypedObject<UniaxialMaterial>(matTag);

    if (material == 0) {
      opserr << "WARNING material does not exist\n";
      opserr << "material: " << matTag;
      opserr << "\nhystereticBackbone Material: " << tag << endln;
      return TCL_ERROR;
    }

    // theBackbone = new MaterialBackbone(tag, *material);
    theBackbone = 0;
  }

  else {
    opserr << "WARNING unknown type of hystereticBackbone: " << argv[1];
    opserr << "\nValid types: Bilinear, Trilinear, Arctangent," << endln;
    opserr << "\tCapped, LinearCapped, Material" << endln;
    return TCL_ERROR;
  }

  // Ensure we have created the Material, out of memory if got here and no
  // material
  if (theBackbone == 0) {
    opserr << "WARNING ran out of memory creating hystereticBackbone\n";
    opserr << argv[1] << endln;
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (!builder->addTypedObject<HystereticBackbone>(std::atoi(argv[2]), theBackbone)) {
    delete theBackbone; // invoke the material objects destructor, otherwise mem leak
    return TCL_ERROR;
  }


  return TCL_OK;
}

HystereticBackbone*
TclCommand_newManderBackbone(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{
  if (argc < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: hystereticBackbone Mander tag? fc? epsc? Ec?" << endln;
    return nullptr;
  }
  
  int tag;
  double fc, epsc, Ec;
  
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid hystereticBackbone Mander tag" << endln;
    return nullptr;
  }
  
  if (Tcl_GetDouble(interp, argv[3], &fc) != TCL_OK) {
    opserr << "WARNING invalid hystereticBackbone Mander fc" << endln;
    return nullptr;
  }
  
  if (Tcl_GetDouble(interp, argv[4], &epsc) != TCL_OK) {
    opserr << "WARNING invalid hystereticBackbone Mander epsc" << endln;
    return nullptr;
  }
  
  if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
    opserr << "WARNING invalid hystereticBackbone Mander Ec" << endln;
    return nullptr;
  }
  
  return new ManderBackbone (tag, fc, epsc, Ec);
}


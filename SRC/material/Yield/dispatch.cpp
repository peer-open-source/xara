#include <tcl.h>
#include "YieldSurface_BC.h"
#include <BasicModelBuilder.h>
#include <string.h>
#include <Vector.h>

#include "NullYS2D.h"
#include "Attalla2D.h"
#include "Orbison2D.h"
#include "Hajjar2D.h"
#include "ElTawil2D.h"
#include "ElTawil2DUnSym.h"

static void
printCommand(int argc, TCL_Char ** const argv)
{
  opserr << "Input command: ";
  for (int i = 0; i < argc; i++)
    opserr << argv[i] << " ";
  opserr << "\n";
}

int
TclBasicBuilderYieldSurface_BCCommand(ClientData clientData, Tcl_Interp *interp,
                                      int argc, TCL_Char ** const argv)
{
  int tag;
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  // Make sure there is a minimum number of arguments
  if (argc < 3) {
    opserr << "WARNING insufficient number of uniaxial material arguments\n";
    opserr << "Want: yieldSurfaceBC type? tag? <specific material args>"
           << "\n";
    return TCL_ERROR;
  }

  // Pointer to a ys that will be added to the model builder
  YieldSurface_BC *theYS = nullptr;

  if (strcmp(argv[1], "null") == 0) {
    if (argc < 4) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: yieldSurfaceBC null tag? dimensions?" << "\n";
      return TCL_ERROR;
    }
    int dim;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC null tag" << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3], &dim) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC null dimensions" << "\n";
      return TCL_ERROR;
    }

    switch (dim) {
    // case 1: 1D YS
    case 2:
      theYS = new NullYS2D(tag);
      break;
      // case 3: 3D YS
    default:
      opserr << "incorrect dimension for null ys\n";
      return TCL_ERROR;
    }

  } else if (strcmp(argv[1], "Orbison2D") == 0) {
    if (argc < 6) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      // Orbison2D(int tag, double xmax, double ymax, YS_HardeningModel &model);
      opserr << "Want: yieldSurfaceBC Orbison2D tag? xCap? yCap? ys_model_tag?"
             << "\n";
      return TCL_ERROR;
    }

    double xCap, yCap;
    // int matID1, matID2;
    int modelID;
    //		double isoRatio;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC Orbison2D tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &xCap) != TCL_OK) {
      opserr << "WARNING invalid xCap\n";
      opserr << "yieldSurfaceBC Orbison2D tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &yCap) != TCL_OK) {
      opserr << "WARNING invalid yCap\n";
      opserr << "yieldSurfaceBC Orbison2D tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[5], &modelID) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC Orbison2D matID1" << modelID
             << "\n";
      return TCL_ERROR;
    }

    YS_Evolution *theModel = builder->getTypedObject<YS_Evolution>(modelID);
    if (theModel == 0) {
      opserr << "WARNING yieldSurfaceBC Orbison2D no ys_model exists with tag: "
             << modelID << "\n";
      return TCL_ERROR;
    }

    // Parsing was successful, allocate the material
    theYS = new Orbison2D(tag, xCap, yCap, *theModel);
  }

  else if (strcmp(argv[1], "ElTawil2D") == 0) {
    if (argc < 7) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      // Orbison2D(int tag, double xmax, double ymax, YS_HardeningModel &model);
      opserr << "Want: yieldSurfaceBC ElTawil2D tag? xCap? yCap? ys_model_tag?"
             << "\n";
      return TCL_ERROR;
    }

    double xBal, yBal;
    double yPos, yNeg;

    int modelID;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC ElTawil2D tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &xBal) != TCL_OK) {
      opserr << "WARNING invalid xBal\n";
      opserr << "yieldSurfaceBC ElTawil2D tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &yBal) != TCL_OK) {
      opserr << "WARNING invalid yBal\n";
      opserr << "yieldSurfaceBC ElTawil2D tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &yPos) != TCL_OK) {
      opserr << "WARNING invalid xPos\n";
      opserr << "yieldSurfaceBC ElTawil2D tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &yNeg) != TCL_OK) {
      opserr << "WARNING invalid yNeg\n";
      opserr << "yieldSurfaceBC ElTawil2D tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &modelID) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC ElTawil2D matID1" << modelID
             << "\n";
      return TCL_ERROR;
    }

    YS_Evolution *theModel = builder->getTypedObject<YS_Evolution>(modelID);
    if (theModel == 0) {
      opserr << "WARNING yieldSurfaceBC ElTawil2D no ys_model exists with tag: "
             << modelID << "\n";
      return TCL_ERROR;
    }

    // Parsing was successful, allocate the material
    theYS = new ElTawil2D(tag, xBal, yBal, yPos, yNeg, *theModel);
  }

  else if (strcmp(argv[1], "ElTawil2DUnSym") == 0) {
    if (argc < 9) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      // Orbison2D(int tag, double xmax, double ymax, YS_HardeningModel &model);
      opserr << "Want: yieldSurfaceBC ElTawil2DUnSym tag? xPosBal? yPosBal? "
             << "xNegBal? yPos? yNeg? ys_model_tag?" << "\n";
      return TCL_ERROR;
    }

    double xPosBal, yPosBal;
    double xNegBal, yNegBal;
    double yPos, yNeg;

    int modelID;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC ElTawil2DUnSym tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &xPosBal) != TCL_OK) {
      opserr << "WARNING invalid xPosBal\n";
      opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &yPosBal) != TCL_OK) {
      opserr << "WARNING invalid yPosBal\n";
      opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &xNegBal) != TCL_OK) {
      opserr << "WARNING invalid xNegBal\n";
      opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &yNegBal) != TCL_OK) {
      opserr << "WARNING invalid yNegBal\n";
      opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &yPos) != TCL_OK) {
      opserr << "WARNING invalid xPos\n";
      opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &yNeg) != TCL_OK) {
      opserr << "WARNING invalid yNeg\n";
      opserr << "yieldSurfaceBC ElTawil2DUnSym tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[9], &modelID) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC ElTawil2DUnSym matID1"
             << modelID << "\n";
      return TCL_ERROR;
    }

    YS_Evolution *theModel = builder->getTypedObject<YS_Evolution>(modelID);
    if (theModel == 0) {
      opserr << "WARNING yieldSurfaceBC ElTawil2D no ys_model exists with tag: "
             << modelID << "\n";
      return TCL_ERROR;
    }

    // Parsing was successful, allocate the material
    theYS = new ElTawil2DUnSym(tag, xPosBal, yPosBal, xNegBal, yNegBal, yPos,
                               yNeg, *theModel);
  }

  else if (strcmp(argv[1], "Attalla2D") == 0) {
    // 	     Attalla2D( int tag, double xmax, double ymax, YS_HardeningModel
    // &model, 					double x_offset=0, double y_offset=0, 					double a01=0.19,  double
    // a02=0.54, double a03=-1.4, 					double a04=-1.64, double a05=2.21, double
    // a06=2.10);

    if (argc < 6 || argc > 14) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: yieldSurfaceBC Attalla2D tag? xCap? yCap? matXTag? "
                "maxYTag? isoRatio? <..>"
             << "\n";
      return TCL_ERROR;
    }

    double xCap, yCap;
    // int matID1, matID2;
    int modelID;
    // double isoRatio;
    // double x_offset = 0, y_offset = 0;
    Vector param(6);

    param[0] = 0.19;
    param[1] = 0.54;
    param[2] = -1.40;
    param[3] = -1.64;
    param[4] = 2.21;
    param[5] = 2.10;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC Attalla2D tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &xCap) != TCL_OK) {
      opserr << "WARNING invalid xCap\n";
      opserr << "yieldSurfaceBC Attalla2D tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &yCap) != TCL_OK) {
      opserr << "WARNING invalid yCap\n";
      opserr << "yieldSurfaceBC Attalla2D tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[5], &modelID) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC Attalla2D modelID" << modelID
             << "\n";
      return TCL_ERROR;
    }

    YS_Evolution *theModel = builder->getTypedObject<YS_Evolution>(modelID);
    if (theModel == 0) {
      opserr << "WARNING yieldSurfaceBC Orbison2D no ys_model exists with tag: "
             << modelID << "\n";
      return TCL_ERROR;
    }

    if (argc > 6) {
      int count = 6;
      double temp;

      for (int i = 0; i < 6; i++) {
        if (Tcl_GetDouble(interp, argv[count], &temp) != TCL_OK) {
          opserr << "WARNING invalid parameter " << i + 1 << "\n";
          opserr << "yieldSurfaceBC Attalla2D tag: " << tag << "\n";
          return TCL_ERROR;
        }
        param(i) = temp;
        count++;
      }
    }

    // Parsing was successful, allocate the material
    theYS = new Attalla2D(tag, xCap, yCap, *theModel, param(0), param(1),
                          param(2), param(3), param(4), param(5));

  }

  else if (strcmp(argv[1], "Hajjar2D") == 0) {
    if (argc < 9) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      opserr
          << "Want: yieldSurfaceBC Hajjar2D tag? ysModelTag? D? b? t? fc? fy?"
          << "\n";
      return TCL_ERROR;
    }

    // int matID1, matID2;
    int modelID;
    //		double isoRatio;
    double D, b, t, fc, fy;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC Hajjar2D  tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &modelID) != TCL_OK) {
      opserr << "WARNING invalid yieldSurfaceBC Hajjar2D  matID1" << modelID
             << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &D) != TCL_OK) {
      opserr << "WARNING invalid D \n";
      opserr << "yieldSurfaceBC Hajjar2D  tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
      opserr << "WARNING invalid b \n";
      opserr << "yieldSurfaceBC Hajjar2D  tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &t) != TCL_OK) {
      opserr << "WARNING invalid t \n";
      opserr << "yieldSurfaceBC Hajjar2D  tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &fc) != TCL_OK) {
      opserr << "WARNING invalid fc \n";
      opserr << "yieldSurfaceBC Hajjar2D  tag: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &fy) != TCL_OK) {
      opserr << "WARNING invalid fy \n";
      opserr << "yieldSurfaceBC Hajjar2D  tag: " << tag << "\n";
      return TCL_ERROR;
    }

    YS_Evolution *theModel = builder->getTypedObject<YS_Evolution>(modelID);
    if (theModel == 0) {
      opserr << "WARNING yieldSurfaceBC Orbison2D no ys_model exists with tag: "
             << modelID << "\n";
      return TCL_ERROR;
    }

    // Parsing was successful, allocate the material
    theYS = new Hajjar2D(tag, *theModel, D, b, t, fc, fy);
  }

  else {
    opserr << "Warning - unknown yield surface type \n";
    printCommand(argc, argv);
  }

  ///////////////////////////////////////////////////////////////
  // Now add the ys to the modelBuilder
  ///////////////////////////////////////////////////////////////

  if (builder->addTaggedObject<YieldSurface_BC>(*theYS) < 0) {
    opserr << "WARNING could not add YieldSurfaceBC to the domain\n";
    opserr << *theYS << "\n";
    delete theYS; // invoke the material objects destructor, otherwise mem leak
    return TCL_ERROR;
  }

  return TCL_OK;
}
#include "YieldSurface_BC.h"
#include <BasicModelBuilder.h>
#include <string.h>
#include <Vector.h>

#include "NullEvolution.h"
#include "Kinematic2D01.h"
#include "PeakOriented2D01.h"
#include "Isotropic2D01.h"
#include "CombinedIsoKin2D01.h"

#include "Kinematic2D02.h"
#include "PeakOriented2D02.h"
#include "CombinedIsoKin2D02.h"

#include "PlasticHardeningMaterial.h"
#include "YieldSurface_BC.h"


static int
addTclYS_Evolution(BasicModelBuilder *theBuilder, YS_Evolution *theModel)
{
  if (theModel == nullptr)
    return TCL_ERROR;

  if (theBuilder->addTaggedObject<YS_Evolution>(*theModel) < 0) {
    opserr << "WARNING could not add hardening model to the domain\n";
    opserr << *theModel << "\n";
    delete theModel; // invoke the material objects destructor, otherwise mem
                     // leak
    return TCL_ERROR;
  }

  return TCL_OK;
}

static PlasticHardeningMaterial *
getTclPlasticMaterial(Tcl_Interp *interp, TCL_Char *arg,
                      BasicModelBuilder *theBuilder)
{
  int id;
  if (Tcl_GetInt(interp, arg, &id) != TCL_OK) {
    opserr << "WARNING: TclModelYS_EvolutionCommand - Invalid plastic material "
              "tag \n";
    return 0;
  }

  PlasticHardeningMaterial *theMat = theBuilder->getTypedObject<PlasticHardeningMaterial>(id);
  if (theMat == 0) {
    opserr << "WARNING: TclModelYS_EvolutionCommand - no "
              "PlasticHardeningMaterial with id = "
           << id << " exists\n";
    return 0;
  } else
    return theMat;
}

static YieldSurface_BC *
getTclYieldSurface_BC(Tcl_Interp *interp, TCL_Char *arg,
                      BasicModelBuilder *builder)
{
  int id;
  if (Tcl_GetInt(interp, arg, &id) != TCL_OK) {
    opserr << "WARNING: TclModelYS_EvolutionCommand - Invalid YieldSurface_BC "
              "tag \n";
    return 0;
  }

  YieldSurface_BC *theYS = builder->getTypedObject<YieldSurface_BC>(id);
  if (theYS == 0) {
    opserr << "WARNING: TclModelYS_EvolutionCommand - no YieldSurface_BC with "
              "id = "
           << id << " exists\n";
    return 0;
  } else
    return theYS;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////////////

static int
TclNullEvolutionCommand(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  YS_Evolution *theModel = 0;
  int tag;
  double isox;
  double isoy;
  double isoz;
  int dim = 0;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
    return TCL_ERROR;

  if (argc > 3) {
    if (Tcl_GetDouble(interp, argv[3], &isox) != TCL_OK)
      return TCL_ERROR;
    dim++;
  }
  if (argc > 4) {
    if (Tcl_GetDouble(interp, argv[4], &isoy) != TCL_OK)
      return TCL_ERROR;
    dim++;
  }
  if (argc > 5) {
    if (Tcl_GetDouble(interp, argv[5], &isoz) != TCL_OK)
      return TCL_ERROR;
    dim++;
  }

  //		opserr << "Dim = " << dim << "\n";
  //		opserr << "\a";

  // Parsing was successful, allocate the material
  if (dim == 1)
    theModel = new NullEvolution(tag, isox);
  else if (dim == 2)
    theModel = new NullEvolution(tag, isox, isoy);
  else if (dim == 3)
    theModel = new NullEvolution(tag, isox, isoy, isoz);
  else
    theModel = 0;

  return addTclYS_Evolution(builder, theModel);
}

static int
TclKinematic2D01Command(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  YS_Evolution *theModel = 0;
  int tag;
  double minIsoFactor;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
    return TCL_ERROR;

  PlasticHardeningMaterial *theMatX =
      getTclPlasticMaterial(interp, argv[4], builder);
  if (theMatX == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *theMatY =
      getTclPlasticMaterial(interp, argv[5], builder);
  if (theMatY == 0)
    return TCL_ERROR;

  double dir;
  if (Tcl_GetDouble(interp, argv[6], &dir) != TCL_OK)
    return TCL_ERROR;

  // Parsing was successful, allocate the material
  theModel = new Kinematic2D01(tag, minIsoFactor, *theMatX, *theMatY, dir);

  return addTclYS_Evolution(builder, theModel);
}

static int
TclIsotropic2D01Command(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  YS_Evolution *theModel = 0;

  int tag;
  double minIsoFactor;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
    return TCL_ERROR;

  PlasticHardeningMaterial *theMatX =
      getTclPlasticMaterial(interp, argv[4], builder);
  if (theMatX == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *theMatY =
      getTclPlasticMaterial(interp, argv[5], builder);
  if (theMatY == 0)
    return TCL_ERROR;

  // Parsing was successful, allocate the material
  theModel = new Isotropic2D01(tag, minIsoFactor, *theMatX, *theMatY);

  return addTclYS_Evolution(builder, theModel);
}

static int
TclPeakOriented2D01Command(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  YS_Evolution *theModel = 0;
  int tag;
  double minIsoFactor;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
    return TCL_ERROR;

  PlasticHardeningMaterial *theMatX =
      getTclPlasticMaterial(interp, argv[4], builder);
  if (theMatX == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *theMatY =
      getTclPlasticMaterial(interp, argv[5], builder);
  if (theMatY == 0)
    return TCL_ERROR;

  // Parsing was successful, allocate the material
  theModel = new PeakOriented2D01(tag, minIsoFactor, *theMatX, *theMatY);

  return addTclYS_Evolution(builder, theModel);
}

int
TclCombinedIsoKin2D01Command(ClientData clientData, Tcl_Interp *interp, int argc,
                             TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  YS_Evolution *theModel = 0;

  int tag;
  double minIsoFactor, iso_ratio, kin_ratio, shr_iso_ratio, shr_kin_ratio;
  int deformable;
  bool deform = false;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
    return TCL_ERROR;
  if (Tcl_GetDouble(interp, argv[3], &iso_ratio) != TCL_OK)
    return TCL_ERROR;
  if (Tcl_GetDouble(interp, argv[4], &kin_ratio) != TCL_OK)
    return TCL_ERROR;
  if (Tcl_GetDouble(interp, argv[5], &shr_iso_ratio) != TCL_OK)
    return TCL_ERROR;
  if (Tcl_GetDouble(interp, argv[6], &shr_kin_ratio) != TCL_OK)
    return TCL_ERROR;
  if (Tcl_GetDouble(interp, argv[7], &minIsoFactor) != TCL_OK)
    return TCL_ERROR;

  PlasticHardeningMaterial *kpx_pos =
      getTclPlasticMaterial(interp, argv[8], builder);
  if (kpx_pos == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *kpx_neg =
      getTclPlasticMaterial(interp, argv[9], builder);
  if (kpx_neg == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *kpy_pos =
      getTclPlasticMaterial(interp, argv[10], builder);
  if (kpx_pos == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *kpy_neg =
      getTclPlasticMaterial(interp, argv[11], builder);
  if (kpx_neg == 0)
    return TCL_ERROR;

  if (Tcl_GetInt(interp, argv[12], &deformable) != TCL_OK)
    return TCL_ERROR;

  double dir;
  if (Tcl_GetDouble(interp, argv[13], &dir) != TCL_OK)
    return TCL_ERROR;

  if (deformable == 1)
    deform = true;

  // Parsing was successful, allocate the material
  theModel = new CombinedIsoKin2D01(tag, iso_ratio, kin_ratio, shr_iso_ratio,
                                    shr_kin_ratio, minIsoFactor, *kpx_pos,
                                    *kpx_neg, *kpy_pos, *kpy_neg, deform, dir);

  return addTclYS_Evolution(builder, theModel);
}

////////////////////////////////////////////////////////////////////////////////////////
static int
TclKinematic2D02Command(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  YS_Evolution *theModel = 0;
  int tag;
  double minIsoFactor;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
    return TCL_ERROR;

  YieldSurface_BC *ys = getTclYieldSurface_BC(interp, argv[4], builder);
  if (ys == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *theMatX =
      getTclPlasticMaterial(interp, argv[5], builder);
  if (theMatX == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *theMatY =
      getTclPlasticMaterial(interp, argv[6], builder);
  if (theMatY == 0)
    return TCL_ERROR;

  int algo;
  double resfact, appfact, dir;

  if (Tcl_GetInt(interp, argv[7], &algo) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[8], &resfact) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[9], &appfact) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[10], &dir) != TCL_OK)
    return TCL_ERROR;

  // Parsing was successful, allocate the material
  theModel = new Kinematic2D02(tag, minIsoFactor, *ys, *theMatX, *theMatY, algo,
                               resfact, appfact, dir);

  return addTclYS_Evolution(builder, theModel);
}

int
TclPeakOriented2D02Command(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  YS_Evolution *theModel = 0;

  int tag;
  double minIsoFactor;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
    return TCL_ERROR;

  YieldSurface_BC *ys = getTclYieldSurface_BC(interp, argv[4], builder);
  if (ys == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *kinX =
      getTclPlasticMaterial(interp, argv[5], builder);
  if (kinX == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *kinY =
      getTclPlasticMaterial(interp, argv[6], builder);
  if (kinY == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *isoX =
      getTclPlasticMaterial(interp, argv[7], builder);
  if (isoX == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *isoY =
      getTclPlasticMaterial(interp, argv[8], builder);
  if (isoY == 0)
    return TCL_ERROR;
  int algo;
  if (Tcl_GetInt(interp, argv[9], &algo) != TCL_OK)
    return TCL_ERROR;

  // Parsing was successful, allocate the material
  theModel = new PeakOriented2D02(tag, minIsoFactor, *ys, *kinX, *kinY, *isoX,
                                  *isoY, algo);

  return addTclYS_Evolution(builder, theModel);
}

static int
TclCombinedIsoKin2D02Command(ClientData clientData, Tcl_Interp *interp, int argc,
                             TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  YS_Evolution *theModel = 0;
  int tag, deformable;
  bool deform = false;
  double minIsoFactor, isoRatio, kinRatio;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[3], &minIsoFactor) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[4], &isoRatio) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[5], &kinRatio) != TCL_OK)
    return TCL_ERROR;

  YieldSurface_BC *ys = getTclYieldSurface_BC(interp, argv[6], builder);
  if (ys == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *kinX =
      getTclPlasticMaterial(interp, argv[7], builder);
  if (kinX == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *kinY =
      getTclPlasticMaterial(interp, argv[8], builder);
  if (kinY == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *isoXPos =
      getTclPlasticMaterial(interp, argv[9], builder);
  if (isoXPos == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *isoXNeg =
      getTclPlasticMaterial(interp, argv[10], builder);
  if (isoXNeg == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *isoYPos =
      getTclPlasticMaterial(interp, argv[11], builder);
  if (isoYPos == 0)
    return TCL_ERROR;

  PlasticHardeningMaterial *isoYNeg =
      getTclPlasticMaterial(interp, argv[12], builder);
  if (isoYNeg == 0)
    return TCL_ERROR;

  if (Tcl_GetInt(interp, argv[13], &deformable) != TCL_OK)
    return TCL_ERROR;

  if (deformable == 1)
    deform = true;
  int algo;
  if (Tcl_GetInt(interp, argv[14], &algo) != TCL_OK)
    return TCL_ERROR;

  double resfact, appfact, dir;

  if (Tcl_GetDouble(interp, argv[15], &resfact) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[16], &appfact) != TCL_OK)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[17], &dir) != TCL_OK)
    return TCL_ERROR;

  // Parsing was successful, allocate the material
  theModel = new CombinedIsoKin2D02(
      tag, minIsoFactor, isoRatio, kinRatio, *ys, *kinX, *kinY, *isoXPos,
      *isoXNeg, *isoYPos, *isoYNeg, deform, algo, resfact, appfact, dir);

  return addTclYS_Evolution(builder, theModel);
}

int
TclBasicBuilderYS_EvolutionModelCommand(ClientData clientData,
                                        Tcl_Interp *interp, int argc,
                                        TCL_Char ** const argv)
{
  // TODO: BasicModelBuilder* theBuilder = (BasicModelBuilder*)clientData;

  if (strcmp(argv[1], "null") == 0) {
    return TclNullEvolutionCommand(clientData, interp, argc, argv);
  } else if (strcmp(argv[1], "kinematic2D01") == 0) {
    return TclKinematic2D01Command(clientData, interp, argc, argv);
  } else if (strcmp(argv[1], "isotropic2D01") == 0) {
    return TclIsotropic2D01Command(clientData, interp, argc, argv);
  } else if (strcmp(argv[1], "peakOriented2D01") == 0) {
    return TclPeakOriented2D01Command(clientData, interp, argc, argv);
  } else if (strcmp(argv[1], "combinedIsoKin2D01") == 0) {
    return TclCombinedIsoKin2D01Command(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "kinematic2D02") == 0) {
    return TclKinematic2D02Command(clientData, interp, argc, argv);
  } else if (strcmp(argv[1], "peakOriented2D02") == 0) {
    return TclPeakOriented2D02Command(clientData, interp, argc, argv);
  } else if (strcmp(argv[1], "combinedIsoKin2D02") == 0) {
    return TclCombinedIsoKin2D02Command(clientData, interp, argc, argv);
  } else {
    opserr << "Unknown YS_Evolution type: " << argv[1] << "\n";
    return TCL_ERROR;
  }
}

///

#include "MultiLinearKp.h"
#include "ExponReducing.h"
#include "NullPlasticMaterial.h"
#include <BasicModelBuilder.h>


static int
TclMultiLinearCommand(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  // Pointer to a uniaxial material that will be added to the model builder
  PlasticHardeningMaterial *theMaterial = 0;

  int tag;
  if (strcmp(argv[1], "multiLinearKp") == 0) {
    int numPoints = (argc - 3) / 2;

    if (numPoints < 2) {
      opserr << "WARNING invalid uniaxialMaterial MultilinearUniaxial tag"
             << "\n";
      opserr << "Minimum of 2 points are required\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial MultilinearUniaxial tag"
             << "\n";
      return TCL_ERROR;
    }

    Vector defo(numPoints);
    Vector kp(numPoints);

    double temp;
    int indx = 3, j1, j2;

    for (j1 = 0; j1 < numPoints; j1++) {
      if (Tcl_GetDouble(interp, argv[indx], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << temp << '\n';
        opserr << "MultilinearUniaxial material: " << tag << "\n";
        return TCL_ERROR;
      }

      defo(j1) = temp;
      indx++;
    }

    for (j2 = 0; j2 < numPoints; j2++) {
      if (Tcl_GetDouble(interp, argv[indx], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << temp << '\n';
        opserr << "MultilinearUniaxial material: " << tag << "\n";
        return TCL_ERROR;
      }

      kp(j2) = temp;
      indx++;
    }

    // Parsing was successful, allocate the material
    theMaterial = new MultiLinearKp(tag, defo, kp);
  }

  // Now add the material to the modelBuilder
  if (builder->addTaggedObject<PlasticHardeningMaterial>(*theMaterial) < 0) {
    opserr << "WARNING could not add uniaxialMaterial to the domain\n";
    opserr << *theMaterial << "\n";
    delete theMaterial; // invoke the material objects destructor, otherwise mem
                        // leak
    return TCL_ERROR;
  }

  return TCL_OK;
}

// QuadrReducing(int tag, double kp0, double kp_half);
/*int TclQuadrReducingCommand(ClientData clientData, Tcl_Interp *interp, int
argc, char **argv, TclBasicBuilder *theTclBuilder)
{
    // Pointer to a uniaxial material that will be added to the model builder
    PlasticHardeningMaterial *theMaterial = 0;

        int tag;
        double kp_0;
        double kp_half;

        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
        {
                opserr << "WARNING invalid  PlaticHardening quadrReducing tag"
<< "\n"; return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[3], &kp_0) != TCL_OK)
        {
                opserr << "WARNING invalid  PlaticHardening quadrReducing kp_0"
<< "\n"; return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[4], &kp_half) != TCL_OK)
        {
                opserr << "WARNING invalid  PlaticHardening quadrReducing
kp_half" << "\n"; return TCL_ERROR;
        }

        theMaterial = new QuadrReducing(tag, kp_0, kp_half);
    if (builder->addRegistryObject("YS_PlasticMaterial", tag, (void*)theMaterial) < 0)
        {
                opserr << "WARNING could not add uniaxialMaterial to the
domain\n"; opserr << *theMaterial << "\n"; delete theMaterial; // invoke the
material objects destructor, otherwise mem leak return TCL_ERROR;
    }

        return TCL_OK;
}
*/

int
TclExponReducingCommand(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (argc < 5) {
    opserr << "TclExponReducingCommand - argc != 5 \n";
    return TCL_ERROR;
  }

  PlasticHardeningMaterial *theMaterial = 0;

  int tag;
  double arg1, arg2, arg3;

  // plasticMaterial exponReducing (int tag, double kp0, double alfa); //5
  // plasticMaterial exponReducing (int tag, double kp0, double x0, double tol);
  // //6
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid  PlaticHardening exponReducing tag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[3], &arg1) != TCL_OK) {
    opserr << "WARNING invalid double PlaticHardening exponReducing" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &arg2) != TCL_OK) {
    opserr << "WARNING invalid double PlaticHardening exponReducing" << "\n";
    return TCL_ERROR;
  }

  if (argc == 6) {
    if (Tcl_GetDouble(interp, argv[5], &arg3) != TCL_OK) {
      opserr << "WARNING invalid double PlaticHardening exponReducing" << "\n";
      return TCL_ERROR;
    }

    theMaterial = new ExponReducing(tag, arg1, arg2, arg3);
    //		opserr << "factor = " << arg3 << "\n";
  } else
    theMaterial = new ExponReducing(tag, arg1, arg2);

  if (builder->addTaggedObject<PlasticHardeningMaterial>(*theMaterial) < 0) {
    opserr << "WARNING could not add uniaxialMaterial to the domain\n";
    opserr << *theMaterial << "\n";
    delete theMaterial; // invoke the material objects destructor, otherwise mem
                        // leak
    return TCL_ERROR;
  }

  return TCL_OK;
}

static int
TclNullPlasticMaterialCommand(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  PlasticHardeningMaterial *theMaterial = 0;

  int tag;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid  PlaticHardening quadrReducing tag" << "\n";
    return TCL_ERROR;
  }

  theMaterial = new NullPlasticMaterial(tag);
  if (builder->addTaggedObject<PlasticHardeningMaterial>(*theMaterial) < 0) {
    opserr << "WARNING could not add uniaxialMaterial to the domain\n";
    opserr << *theMaterial << "\n";
    delete theMaterial; // invoke the material objects destructor, otherwise mem
                        // leak
    return TCL_ERROR;
  }

  return TCL_OK;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int
TclBasicBuilderPlasticMaterialCommand(ClientData clientData, Tcl_Interp *interp,
                                      int argc, TCL_Char ** const argv)
{

  if (strcmp(argv[1], "multiLinearKp") == 0) {
    return TclMultiLinearCommand(clientData, interp, argc, argv);
  }
#if 0
  else if (strcmp(argv[1],"quadrReducing") == 0) {
          return TclQuadrReducingCommand(clientData, interp, argc, argv, theTclBuilder);
  }
#endif
  else if (strcmp(argv[1], "exponReducing") == 0) {
    return TclExponReducingCommand(clientData, interp, argc, argv);

  } else if (strcmp(argv[1], "null") == 0) {
    return TclNullPlasticMaterialCommand(clientData, interp, argc, argv);
  } else {
    opserr << "Unknown PlasticMaterial: \nValid types: null, multiLinearKp, "
           << "quadrReducing, exponReducing \n";

    return TCL_ERROR;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// @ rkaul@stanford.edu
// @ ggd@stanford.edu

#include <BasicModelBuilder.h>
#include <string.h>
#include <Vector.h>

#include <YieldSurface_BC.h>
#include <YS_Section2D01.h>
#include <YS_Section2D02.h>

#include <SoilFootingSection2d.h>

// Added by S.Gajan <sgajan@ucdavis.edu>

SectionForceDeformation *
TclBasicBuilderYS_SectionCommand(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*) clientData;

  if (argc < 3) {
    opserr << "WARNING insufficient number of arguments\n";
    printCommand(argc, argv);
    return 0;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid section tag\n";
    printCommand(argc, argv);
    return 0;
  }

  SectionForceDeformation *theModel = nullptr;

  if (strcmp(argv[1], "YS_Section2D01") == 0 ||
      strcmp(argv[1], "YS_Section2d01") == 0) {

    if (argc < 7) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: section YS_Section2D01 tag? E? A? Iz? ysTag? <algo?>"
             << "\n";
      return 0;
    }

    int algo, ysTag;
    double E, A, Iz;
    int indx = 3;

    if (Tcl_GetDouble(interp, argv[indx++], &E) != TCL_OK) {
      opserr << "WARNING invalid E" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &A) != TCL_OK) {
      opserr << "WARNING invalid A" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &Iz) != TCL_OK) {
      opserr << "WARNING invalid Iz" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetInt(interp, argv[indx++], &ysTag) != TCL_OK) {
      opserr << "WARNING invalid ysTag" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    YieldSurface_BC *ys = builder->getTypedObject<YieldSurface_BC>(ysTag);

    if (ys == 0) {
      opserr << "WARNING yield surface does not exist\n";
      opserr << "yieldSurface: " << ysTag;
      opserr << "\nsection YieldSurface: " << tag << "\n";
      return 0;
    }

    bool useKr = true;
    if (argc > indx) {
      if (Tcl_GetInt(interp, argv[indx++], &algo) != TCL_OK) {
        opserr << "WARNING invalid algo" << "\n";
        opserr << " section: " << tag << "\n";
        return 0;
      }
      if (algo == 0)
        useKr = false;
    }

    theModel = new YS_Section2D01(tag, E, A, Iz, ys, useKr);
  }

  else if (strcmp(argv[1], "YS_Section2D02") == 0 ||
           strcmp(argv[1], "YS_Section2d02") == 0) {

    if (argc < 8) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: section YS_Section2D01 tag? E? A? Iz? maxPlastRot? "
                "ysTag? <algo?>"
             << "\n";
      return 0;
    }

    int algo, ysTag;
    double E, A, Iz, maxPlstkRot;
    int indx = 3;

    if (Tcl_GetDouble(interp, argv[indx++], &E) != TCL_OK) {
      opserr << "WARNING invalid E" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &A) != TCL_OK) {
      opserr << "WARNING invalid A" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &Iz) != TCL_OK) {
      opserr << "WARNING invalid Iz" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &maxPlstkRot) != TCL_OK) {
      opserr << "WARNING maxPlstkRot " << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetInt(interp, argv[indx++], &ysTag) != TCL_OK) {
      opserr << "WARNING invalid ysTag" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    YieldSurface_BC *ys = builder->getTypedObject<YieldSurface_BC>(ysTag);

    if (ys == 0) {
      opserr << "WARNING yield surface does not exist\n";
      opserr << "yieldSurface: " << ysTag;
      opserr << "\nsection YieldSurface: " << tag << "\n";
      return 0;
    }

    bool useKr = true;
    if (argc > indx) {
      if (Tcl_GetInt(interp, argv[indx++], &algo) != TCL_OK) {
        opserr << "WARNING invalid algo" << "\n";
        opserr << " section: " << tag << "\n";
        return 0;
      }
      if (algo == 0)
        useKr = false;
    }

    theModel = new YS_Section2D02(tag, E, A, Iz, maxPlstkRot, ys, useKr);
  }

  // Added by S.Gajan <sgajan@ucdavis.edu>

  else if ((strcmp(argv[1], "soilFootingSection2d") == 0) ||
           (strcmp(argv[1], "SoilFootingSection2d") == 0)) {

    if (argc < 10) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc, argv);
      opserr << "Want: section soilFootingSection2d tag? FS? Vult? L? Kv? dL?"
             << "\n";
      return 0;
    }

    double FS, Vult, L, Kv, Kh, Rv, deltaL;
    int indx = 3;

    if (Tcl_GetDouble(interp, argv[indx++], &FS) != TCL_OK) {
      opserr << "WARNING invalid FS" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &Vult) != TCL_OK) {
      opserr << "WARNING invalid Vult" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &L) != TCL_OK) {
      opserr << "WARNING invalid L" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &Kv) != TCL_OK) {
      opserr << "WARNING invalid Kv" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &Kh) != TCL_OK) {
      opserr << "WARNING invalid Kh" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &Rv) != TCL_OK) {
      opserr << "WARNING invalid Rv" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    if (Tcl_GetDouble(interp, argv[indx++], &deltaL) != TCL_OK) {
      opserr << "WARNING invalid Kv" << "\n";
      opserr << " section: " << tag << "\n";
      return 0;
    }

    theModel = new SoilFootingSection2d(tag, FS, Vult, L, Kv, Kh, Rv, deltaL);
  }

  return theModel;
}

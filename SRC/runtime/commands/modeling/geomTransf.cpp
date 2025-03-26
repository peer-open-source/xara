//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: Geometric transformation command
//
// cmp
//
#include <tcl.h>
#include <string.h>
#include <assert.h>
#include <BasicModelBuilder.h>
#include <G3_Logging.h>

#include <LinearCrdTransf2d.h>
#include <LinearCrdTransf3d.h>
#include <PDeltaCrdTransf2d.h>
#include <PDeltaCrdTransf3d.h>
#include <CorotCrdTransf2d.h>
#include <CorotCrdTransf3d.h>
#include <LinearCrdTransf2dInt.h>
#include <CorotCrdTransfWarping2d.h>

#include <LinearFrameTransf3d.h>
#include <PDeltaFrameTransf3d.h>
#include <CorotFrameTransf3d.h>
#include <CorotFrameTransf3d03.h>

//
// Create a coordinate transformation
//
int
TclCommand_addGeomTransf(ClientData clientData, Tcl_Interp *interp, int argc,
                         const char ** const argv)

{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "insufficient number of arguments\n";
    return TCL_ERROR;
  }

  int ndm = builder->getNDM();
  int ndf = builder->getNDF(); // number of degrees of freedom per node
  //
  // 2D Case
  //
  if ((ndm == 2 && ndf == 3) || (ndm == 2 && ndf == 4)) {

    int tag;
    Vector jntOffsetI(2),
           jntOffsetJ(2);

    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      return TCL_ERROR;
    }

    int argi = 2;
    if (Tcl_GetInt(interp, argv[argi++], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid tag\n";
      return TCL_ERROR;
    }

    // Additional options at end of command

    while (argi != argc) {
      if (strcmp(argv[argi], "-jntOffset") == 0) {
        argi++;
        for (int i = 0; i < 2; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetI(i)) != TCL_OK) {
            opserr << G3_ERROR_PROMPT << "invalid jntOffset value\n";
            return TCL_ERROR;
          }
        }

        for (int i = 0; i < 2; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetJ(i)) != TCL_OK) {
            opserr << G3_ERROR_PROMPT << "invalid jntOffset value\n";
            return TCL_ERROR;
          }
        }
      }

      else {
        opserr << G3_ERROR_PROMPT << "unexpected argument " << argv[argi] << "\n";
        return TCL_ERROR;
      }
    }

    //
    // construct the transformation
    //

    FrameTransform2d *crdTransf2d = nullptr;

    if (strcmp(argv[1], "Linear") == 0)
        crdTransf2d = new LinearCrdTransf2d(tag, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "LinearInt") == 0)
      crdTransf2d =
          new LinearCrdTransf2dInt(tag, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "PDelta") == 0 ||
             strcmp(argv[1], "LinearWithPDelta") == 0)
      crdTransf2d = new PDeltaCrdTransf2d(tag, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "Corotational") == 0 && ndf == 3)
      crdTransf2d = new CorotCrdTransf2d(tag, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "Corotational") == 0 && ndf == 4)
      crdTransf2d =
          new CorotCrdTransfWarping2d(tag, jntOffsetI, jntOffsetJ);

    else {
      opserr << G3_ERROR_PROMPT << "invalid Type: " << argv[1] << "\n";
      return TCL_ERROR;
    }

    //
    if (builder->addTaggedObject<FrameTransform2d>(*crdTransf2d) != TCL_OK)
      return TCL_ERROR;
  }

  else if (ndm == 3 && ndf >= 6) {
    int tag;
    Vector vecxzPlane(3);                // vector that defines local xz plane
    Vector jntOffsetI(3), jntOffsetJ(3); // joint offsets in global coordinates

    if (argc < 6) {
      opserr << G3_ERROR_PROMPT 
             << "insufficient arguments\n";
      return TCL_ERROR;
    }

    int argi = 2;
    if (Tcl_GetInt(interp, argv[argi++], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid tag\n";
      return TCL_ERROR;
    }

    // parse orientation vector
    bool parsed_xz = false;
    if (!parsed_xz) {
      const char ** xzarg;
      int xznum;
      Tcl_SplitList(interp, argv[argi], &xznum, &xzarg);
      if (xznum == 3) {
        for (int i=0; i<3; ++i)
           if (Tcl_GetDouble(interp, xzarg[i], &vecxzPlane(i)) != TCL_OK) {
             opserr << G3_ERROR_PROMPT << "Failed  to parse vectxz\n";
             return TCL_ERROR;
           }

        argi++;
        parsed_xz = true;
      }
    } 

    if (!parsed_xz) {
      if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(0)) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid vecxzPlaneX\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(1)) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid vecxzPlaneY\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(2)) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid vecxzPlaneZ\n";
        return TCL_ERROR;
      }
    }

    //
    // Additional keyword options
    //

    while (argi != argc) {
      if (strcmp(argv[argi], "-jntOffset") == 0) {
        argi++;
        for (int i = 0; i < 3; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetI(i)) != TCL_OK) {
            opserr << G3_ERROR_PROMPT << "invalid jntOffset\n";
            return TCL_ERROR;
          }
        }

        for (int i = 0; i < 3; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetJ(i)) != TCL_OK) {
            opserr << G3_ERROR_PROMPT << "invalid jntOffset\n";
            return TCL_ERROR;
          }
        }
      } else {
        opserr << G3_ERROR_PROMPT << "unexpected argument: " << argv[argi] << "\n";
        return TCL_ERROR;
      }
    }

    //
    // construct the transformation
    //
    FrameTransform3d *crdTransf3d=nullptr;

    if (strcmp(argv[1], "Linear") == 0)
      if (!getenv("CRD"))
        crdTransf3d = new LinearFrameTransf3d(tag, vecxzPlane, jntOffsetI, jntOffsetJ);
      else
        crdTransf3d = new LinearCrdTransf3d(tag, vecxzPlane, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "PDelta") == 0 ||
             strcmp(argv[1], "LinearWithPDelta") == 0)
      if (!getenv("CRD"))
        crdTransf3d = new PDeltaFrameTransf3d(tag, vecxzPlane, jntOffsetI, jntOffsetJ);
      else
        crdTransf3d = new PDeltaCrdTransf3d(tag, vecxzPlane, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "Corotational") == 0)
      if (getenv("CRD03")) {
        crdTransf3d = new CorotFrameTransf3d03(tag, vecxzPlane, jntOffsetI, jntOffsetJ);
      }
      else
        crdTransf3d = new CorotFrameTransf3d(tag, vecxzPlane, jntOffsetI, jntOffsetJ);
//    else
//      crdTransf3d = new CorotCrdTransf3d(tag, vecxzPlane, jntOffsetI, jntOffsetJ);

    else {
      opserr << G3_ERROR_PROMPT << "invalid Type\n";
      return TCL_ERROR;
    }

    if (crdTransf3d == nullptr) {
      opserr << G3_ERROR_PROMPT << "Failed to create transform\n";
      return TCL_ERROR;
    }

    // add the transformation to the modelBuilder
    if (builder->addTaggedObject<FrameTransform3d>(*crdTransf3d) != TCL_OK) {
      opserr << G3_ERROR_PROMPT 
             << "Failed to add transformation to model\n";
      return TCL_ERROR;
    }

  } else {
    opserr << G3_ERROR_PROMPT 
           << "ndm = " << ndm << " and ndf = " << ndf
           << " is incompatible with available frame elements\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

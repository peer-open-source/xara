/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the flatSliderBearing element.
//
#include <assert.h>
#include <BasicModelBuilder.h>
#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>

#include <FlatSliderSimple2d.h>
#include <FlatSliderSimple3d.h>
#include <FrictionModel.h>
#include <UniaxialMaterial.h>


#include <api/runtimeAPI.h>


int
TclCommand_addFlatSliderBearing(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char **const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  Domain* theTclDomain = builder->getDomain();
  const int eleArgStart = 1;


  Element *theElement = nullptr;
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();
  int tag;

  if (ndm == 2) {
    // check plane frame problem has 3 dof per node
    if (ndf != 3) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for plane problem need 3 - flatSliderBearing\n";
      return TCL_ERROR;
    }

    // check the number of arguments is correct
    if ((argc - eleArgStart) < 10) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: flatSliderBearing eleTag iNode jNode frnMdlTag kInit -P "
                "matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-shearDist "
                "sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>\n";
      return TCL_ERROR;
    }

    // get the id and end nodes
    int iNode, jNode, frnMdlTag, matTag, argi, i, j;
    int recvMat = 0;
    double kInit;
    double shearDistI = 0.0;
    int doRayleigh = 0;
    double mass = 0.0;
    int maxIter = 25;
    double tol = 1e-12;
    // double kFactUplift = 1e-12;

    if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
      opserr << "WARNING invalid flatSliderBearing eleTag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[4 + eleArgStart], &frnMdlTag) != TCL_OK) {
      opserr << "WARNING invalid frnMdlTag\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    FrictionModel *theFrnMdl = builder->getTypedObject<FrictionModel>(frnMdlTag);
    if (theFrnMdl == 0) {
      opserr << "WARNING friction model not found\n";
      opserr << "frictionModel: " << frnMdlTag << endln;
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5 + eleArgStart], &kInit) != TCL_OK) {
      opserr << "WARNING invalid kInit\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    UniaxialMaterial *theMaterials[2];
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-P") == 0) {
        theMaterials[0] = 0;
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid matTag\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[0] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[0] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-Mz") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid matTag\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[1] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[1] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    if (recvMat != 2) {
      opserr << "WARNING wrong number of materials\n";
      opserr << "got " << recvMat << " materials, but want 2 materials\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // check for optional arguments
    Vector x = 0;
    Vector y = 0;
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-orient") == 0) {
        j = i + 1;
        int numOrient = 0;
        while (j < argc && strcmp(argv[j], "-shearDist") != 0 &&
               strcmp(argv[j], "-doRayleigh") != 0 &&
               strcmp(argv[j], "-mass") != 0 && strcmp(argv[j], "-iter") != 0) {
          numOrient++;
          j++;
        }
        if (numOrient == 6) {
          argi = i + 1;
          x.resize(3);
          y.resize(3);
          double value;
          // read the x values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "flatSliderBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              x(j) = value;
            }
          }
          // read the y values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "flatSliderBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              y(j) = value;
            }
          }
        } else {
          opserr << "WARNING insufficient arguments after -orient flag\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-shearDist") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &shearDistI) != TCL_OK) {
          opserr << "WARNING invalid -shearDist value\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-doRayleigh") == 0)
        doRayleigh = 1;
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-mass") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK) {
          opserr << "WARNING invalid -mass value\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 2 < argc && strcmp(argv[i], "-iter") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &maxIter) != TCL_OK) {
          opserr << "WARNING invalid maxIter\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i + 2], &tol) != TCL_OK) {
          opserr << "WARNING invalid tol\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }

    // now create the flatSliderBearing
    theElement = new FlatSliderSimple2d(tag, iNode, jNode, *theFrnMdl, kInit,
                                        theMaterials, y, x, shearDistI,
                                        doRayleigh, mass, maxIter, tol);

    if (theElement == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // then add the flatSliderBearing to the domain
    if (theTclDomain->addElement(theElement) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      delete theElement;
      return TCL_ERROR;
    }
  }

  else if (ndm == 3) {
    // check space frame problem has 6 dof per node
    if (ndf != 6) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for space problem need 6 - flatSliderBearing \n";
      return TCL_ERROR;
    }

    // check the number of arguments is correct
    if ((argc - eleArgStart) < 14) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: flatSliderBearing eleTag iNode jNode frnMdlTag kInit -P "
                "matTag -T matTag -My matTag -Mz matTag <-orient <x1 x2 x3> y1 "
                "y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter "
                "maxIter tol>\n";
      return TCL_ERROR;
    }

    // get the id and end nodes
    int iNode, jNode, frnMdlTag, matTag, argi, i, j;
    int recvMat = 0;
    double kInit;
    double shearDistI = 0.0;
    int doRayleigh = 0;
    double mass = 0.0;
    int maxIter = 25;
    double tol = 1e-12;
    double kFactUplift = 1e-12;

    if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
      opserr << "WARNING invalid flatSliderBearing eleTag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[4 + eleArgStart], &frnMdlTag) != TCL_OK) {
      opserr << "WARNING invalid frnMdlTag\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    FrictionModel *theFrnMdl = builder->getTypedObject<FrictionModel>(frnMdlTag);
    if (theFrnMdl == 0) {
      opserr << "WARNING friction model not found\n";
      opserr << "frictionModel: " << frnMdlTag << endln;
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5 + eleArgStart], &kInit) != TCL_OK) {
      opserr << "WARNING invalid kInit\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    UniaxialMaterial *theMaterials[4];
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-P") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid axial matTag\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[0] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[0] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-T") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid torsional matTag\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[1] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[1] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-My") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid moment y matTag\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[2] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[2] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-Mz") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid moment z matTag\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[3] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[3] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    if (recvMat != 4) {
      opserr << "WARNING wrong number of materials\n";
      opserr << "got " << recvMat << " materials, but want 4 materials\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // check for optional arguments
    Vector x(0);
    Vector y(3);
    y(0) = 0.0;
    y(1) = 1.0;
    y(2) = 0.0;
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-orient") == 0) {
        j = i + 1;
        int numOrient = 0;
        while (j < argc && strcmp(argv[j], "-shearDist") != 0 &&
               strcmp(argv[j], "-doRayleigh") != 0 &&
               strcmp(argv[j], "-mass") != 0 && strcmp(argv[j], "-iter") != 0 &&
               strcmp(argv[j], "-kFactUplift") != 0) {
          numOrient++;
          j++;
        }
        if (numOrient == 3) {
          argi = i + 1;
          double value;
          // read the y values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "flatSliderBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              y(j) = value;
            }
          }
        } else if (numOrient == 6) {
          argi = i + 1;
          x.resize(3);
          double value;
          // read the x values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "flatSliderBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              x(j) = value;
            }
          }
          // read the y values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "flatSliderBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              y(j) = value;
            }
          }
        } else {
          opserr << "WARNING insufficient arguments after -orient flag\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-shearDist") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &shearDistI) != TCL_OK) {
          opserr << "WARNING invalid -shearDist value\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-doRayleigh") == 0)
        doRayleigh = 1;
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-mass") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK) {
          opserr << "WARNING invalid -mass value\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 2 < argc && strcmp(argv[i], "-iter") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &maxIter) != TCL_OK) {
          opserr << "WARNING invalid maxIter\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i + 2], &tol) != TCL_OK) {
          opserr << "WARNING invalid tol\n";
          opserr << "flatSliderBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-kFactUplift") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &kFactUplift) != TCL_OK) {
          opserr << "WARNING invalid kFactUplift\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }

    // now create the flatSliderBearing
    theElement = new FlatSliderSimple3d(
        tag, iNode, jNode, *theFrnMdl, kInit, theMaterials, y, x, shearDistI,
        doRayleigh, mass, maxIter, tol, kFactUplift);

    if (theElement == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // then add the flatSliderBearing to the domain
    if (theTclDomain->addElement(theElement) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "flatSliderBearing element: " << tag << endln;
      delete theElement;
      return TCL_ERROR;
    }
  }

  else {
    opserr << "WARNING flatSliderBearing command only works when ndm is 2 or "
              "3, ndm: ";
    opserr << ndm << endln;
    return TCL_ERROR;
  }

  // if get here we have successfully created the flatSliderBearing and added it
  // to the domain
  return TCL_OK;
}


// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 12/13
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the RJWatsonEqsBearing element.

#include <BasicModelBuilder.h>

#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>

#include <RJWatsonEQS2d.h>
#include <RJWatsonEQS3d.h>
#include <FrictionModel.h>
#include <UniaxialMaterial.h>


int
TclBasicBuilder_addRJWatsonEqsBearing(ClientData clientData, Tcl_Interp *interp,
                                      int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  constexpr int eleArgStart = 1;


  Element *theElement = 0;
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();
  int tag;

  if (ndm == 2) {
    // check plane frame problem has 3 dof per node
    if (ndf != 3) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for plane problem need 3 - RJWatsonEqsBearing\n";
      return TCL_ERROR;
    }

    // check the number of arguments is correct
    if ((argc - eleArgStart) < 12) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: RJWatsonEqsBearing eleTag iNode jNode frnMdlTag kInit "
                "-P matTag -Vy matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> "
                "<-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter "
                "tol>\n";
      return TCL_ERROR;
    }

    // get the id and end nodes
    int iNode, jNode, frnMdlTag, matTag, argi, i, j;
    int recvMat = 0;
    double kInit;
    double shearDistI = 1.0;
    int doRayleigh = 0;
    double mass = 0.0;
    int maxIter = 25;
    double tol = 1E-12;
    double kFactUplift = 1E-12;

    if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
      opserr << "WARNING invalid RJWatsonEqsBearing eleTag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[4 + eleArgStart], &frnMdlTag) != TCL_OK) {
      opserr << "WARNING invalid frnMdlTag\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    FrictionModel *theFrnMdl = builder->getTypedObject<FrictionModel>(frnMdlTag);
    if (theFrnMdl == 0) {
      opserr << "WARNING friction model not found\n";
      opserr << "frictionModel: " << frnMdlTag << endln;
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5 + eleArgStart], &kInit) != TCL_OK) {
      opserr << "WARNING invalid kInit\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    UniaxialMaterial *theMaterials[3];
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-P") == 0) {
        theMaterials[0] = 0;
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid axial matTag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[0] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[0] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-Vy") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid shear y matTag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[1] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[1] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-Mz") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid moment z matTag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[2] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[2] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    if (recvMat != 3) {
      opserr << "WARNING wrong number of materials\n";
      opserr << "got " << recvMat << " materials, but want 3 materials\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // check for optional arguments
    Vector x = 0;
    Vector y = 0;
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-orient") == 0) {
        j = i + 1;
        int numOrient = 0;
        while (j < argc && strcmp(argv[j], "-shearDist") != 0 &&
               strcmp(argv[j], "-doRayleigh") != 0 &&
               strcmp(argv[j], "-mass") != 0 && strcmp(argv[j], "-iter") != 0 &&
               strcmp(argv[j], "-kFactUplift") != 0) {
          numOrient++;
          j++;
        }
        if (numOrient == 6) {
          argi = i + 1;
          x.resize(3);
          y.resize(3);
          double value;
          // read the x values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "RJWatsonEqsBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              x(j) = value;
            }
          }
          // read the y values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "RJWatsonEqsBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              y(j) = value;
            }
          }
        } else {
          opserr << "WARNING insufficient arguments after -orient flag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-shearDist") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &shearDistI) != TCL_OK) {
          opserr << "WARNING invalid -shearDist value\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-doRayleigh") == 0)
        doRayleigh = 1;
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-mass") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK) {
          opserr << "WARNING invalid -mass value\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 2 < argc && strcmp(argv[i], "-iter") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &maxIter) != TCL_OK) {
          opserr << "WARNING invalid maxIter\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i + 2], &tol) != TCL_OK) {
          opserr << "WARNING invalid tol\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-kFactUplift") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &kFactUplift) != TCL_OK) {
          opserr << "WARNING invalid kFactUplift\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }

    // now create the RJWatsonEqsBearing
    theElement = new RJWatsonEQS2d(tag, iNode, jNode, *theFrnMdl, kInit,
                                   theMaterials, y, x, shearDistI, doRayleigh,
                                   mass, maxIter, tol, kFactUplift);

    if (theElement == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // then add the RJWatsonEqsBearing to the domain
    if (builder->getDomain()->addElement(theElement) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      delete theElement;
      return TCL_ERROR;
    }
  }

  else if (ndm == 3) {
    // check space frame problem has 6 dof per node
    if (ndf != 6) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for space problem need 6 - RJWatsonEqsBearing \n";
      return TCL_ERROR;
    }

    // check the number of arguments is correct
    if ((argc - eleArgStart) < 18) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: RJWatsonEqsBearing eleTag iNode jNode frnMdlTag kInit "
                "-P matTag -Vy matTag -Vz matTag -T matTag -My matTag -Mz "
                "matTag <-orient <x1 x2 x3> y1 y2 y3> <-shearDist sDratio> "
                "<-doRayleigh> <-mass m> <-iter maxIter tol>\n";
      return TCL_ERROR;
    }

    // get the id and end nodes
    int iNode, jNode, frnMdlTag, matTag, argi, i, j;
    int recvMat = 0;
    double kInit;
    double shearDistI = 1.0;
    int doRayleigh = 0;
    double mass = 0.0;
    int maxIter = 25;
    double tol = 1E-12;
    double kFactUplift = 1E-12;

    if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
      opserr << "WARNING invalid RJWatsonEqsBearing eleTag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[4 + eleArgStart], &frnMdlTag) != TCL_OK) {
      opserr << "WARNING invalid frnMdlTag\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    FrictionModel *theFrnMdl = builder->getTypedObject<FrictionModel>(frnMdlTag);
    if (theFrnMdl == 0) {
      opserr << "WARNING friction model not found\n";
      opserr << "frictionModel: " << frnMdlTag << endln;
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5 + eleArgStart], &kInit) != TCL_OK) {
      opserr << "WARNING invalid kInit\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    UniaxialMaterial *theMaterials[6];
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-P") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid axial matTag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[0] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[0] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-Vy") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid shear y matTag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[1] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[1] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-Vz") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid shear z matTag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[2] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[2] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-T") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid torsional matTag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[3] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[3] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-My") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid moment y matTag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[4] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[4] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-Mz") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid moment z matTag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[5] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[5] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    if (recvMat != 6) {
      opserr << "WARNING wrong number of materials\n";
      opserr << "got " << recvMat << " materials, but want 6 materials\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // check for optional arguments
    Vector x(0);
    Vector y(3);
    y(0) = 0.0;
    y(1) = 1.0;
    y(2) = 0.0;
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-orient") == 0) {
        j = i + 1;
        int numOrient = 0;
        while (j < argc && strcmp(argv[j], "-shearDist") != 0 &&
               strcmp(argv[j], "-doRayleigh") != 0 &&
               strcmp(argv[j], "-mass") != 0 && strcmp(argv[j], "-iter") != 0 &&
               strcmp(argv[j], "-kFactUplift") != 0) {
          numOrient++;
          j++;
        }
        if (numOrient == 3) {
          argi = i + 1;
          double value;
          // read the y values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "RJWatsonEqsBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              y(j) = value;
            }
          }
        } else if (numOrient == 6) {
          argi = i + 1;
          x.resize(3);
          double value;
          // read the x values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "RJWatsonEqsBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              x(j) = value;
            }
          }
          // read the y values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "RJWatsonEqsBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              y(j) = value;
            }
          }
        } else {
          opserr << "WARNING insufficient arguments after -orient flag\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-shearDist") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &shearDistI) != TCL_OK) {
          opserr << "WARNING invalid -shearDist value\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-doRayleigh") == 0)
        doRayleigh = 1;
    }
    for (i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-mass") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK) {
          opserr << "WARNING invalid -mass value\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 2 < argc && strcmp(argv[i], "-iter") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &maxIter) != TCL_OK) {
          opserr << "WARNING invalid maxIter\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i + 2], &tol) != TCL_OK) {
          opserr << "WARNING invalid tol\n";
          opserr << "RJWatsonEqsBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 6 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-kFactUplift") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &kFactUplift) != TCL_OK) {
          opserr << "WARNING invalid kFactUplift\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }

    // now create the RJWatsonEqsBearing
    theElement = new RJWatsonEQS3d(tag, iNode, jNode, *theFrnMdl, kInit,
                                   theMaterials, y, x, shearDistI, doRayleigh,
                                   mass, maxIter, tol, kFactUplift);

    if (theElement == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // then add the RJWatsonEqsBearing to the domain
    if (builder->getDomain()->addElement(theElement) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "RJWatsonEqsBearing element: " << tag << endln;
      delete theElement;
      return TCL_ERROR;
    }
  }

  else {
    opserr << "WARNING RJWatsonEqsBearing command only works when ndm is 2 or "
              "3, ndm: ";
    opserr << ndm << endln;
    return TCL_ERROR;
  }

  // if get here we have successfully created the RJWatsonEqsBearing and added
  // it to the domain
  return TCL_OK;
}
/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Description: This file contains the function to parse the TCL input
// for the singleFPBearing element.
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
#include <assert.h>
#include <BasicModelBuilder.h>

#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>

#include <SingleFPSimple2d.h>
#include <SingleFPSimple3d.h>

#include <FrictionModel.h>
#include <UniaxialMaterial.h>


int
TclCommand_addSingleFPBearing(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  Domain* theTclDomain = builder->getDomain();
  const int eleArgStart = 1;


  Element *theElement = 0;
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();
  int tag;

  if (ndm == 2) {
    // check plane frame problem has 3 dof per node
    if (ndf != 3) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for plane problem need 3 - singleFPBearing\n";
      return TCL_ERROR;
    }

    // check the number of arguments is correct
    if ((argc - eleArgStart) < 11) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: singleFPBearing eleTag iNode jNode frnMdlTag Reff kInit "
                "-P matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-shearDist "
                "sDratio> <-doRayleigh> <-inclVertDisp> <-mass m> <-iter "
                "maxIter tol>\n";
      return TCL_ERROR;
    }

    // get the id and end nodes
    int iNode, jNode, frnMdlTag, matTag, argi, i, j;
    int recvMat = 0;
    double Reff, kInit;
    double shearDistI = 0.0;
    int doRayleigh = 0;
    int inclVertDisp = 0;
    double mass = 0.0;
    int maxIter = 25;
    double tol = 1E-12;
    double kFactUplift = 1E-12;

    if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
      opserr << "WARNING invalid singleFPBearing eleTag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[4 + eleArgStart], &frnMdlTag) != TCL_OK) {
      opserr << "WARNING invalid frnMdlTag\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    FrictionModel *theFrnMdl = builder->getTypedObject<FrictionModel>(frnMdlTag);
    if (theFrnMdl == 0) {
      opserr << "WARNING friction model not found\n";
      opserr << "frictionModel: " << frnMdlTag << endln;
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5 + eleArgStart], &Reff) != TCL_OK) {
      opserr << "WARNING invalid Reff\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6 + eleArgStart], &kInit) != TCL_OK) {
      opserr << "WARNING invalid kInit\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    UniaxialMaterial *theMaterials[2];
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-P") == 0) {
        theMaterials[0] = 0;
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid matTag\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[0] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[0] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-Mz") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid matTag\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[1] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[1] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    if (recvMat != 2) {
      opserr << "WARNING wrong number of materials\n";
      opserr << "got " << recvMat << " materials, but want 2 materials\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // check for optional arguments
    Vector x = 0;
    Vector y = 0;
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-orient") == 0) {
        j = i + 1;
        int numOrient = 0;
        while (j < argc && strcmp(argv[j], "-shearDist") != 0 &&
               strcmp(argv[j], "-doRayleigh") != 0 &&
               strcmp(argv[j], "-inclVertDisp") != 0 &&
               strcmp(argv[j], "-mass") != 0 && strcmp(argv[j], "-iter") != 0 &&
               strcmp(argv[j], "-kFactUplift") != 0) {
          numOrient++;
          j++;
        }
        if (numOrient == 6) {
          argi = i + 1;
          x.resize(3);
          y.resize(3);
          double value;
          // read the x values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "singleFPBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              x(j) = value;
            }
          }
          // read the y values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "singleFPBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              y(j) = value;
            }
          }
        } else {
          opserr << "WARNING insufficient arguments after -orient flag\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-shearDist") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &shearDistI) != TCL_OK) {
          opserr << "WARNING invalid -shearDist value\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 7 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-doRayleigh") == 0)
        doRayleigh = 1;
    }
    for (int i = 7 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-inclVertDisp") == 0)
        inclVertDisp = 1;
    }
    for (int i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-mass") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK) {
          opserr << "WARNING invalid -mass value\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 7 + eleArgStart; i < argc; i++) {
      if (i + 2 < argc && strcmp(argv[i], "-iter") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &maxIter) != TCL_OK) {
          opserr << "WARNING invalid maxIter\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i + 2], &tol) != TCL_OK) {
          opserr << "WARNING invalid tol\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-kFactUplift") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &kFactUplift) != TCL_OK) {
          opserr << "WARNING invalid kFactUplift\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }

    // now create the singleFPBearing
    theElement = new SingleFPSimple2d(
        tag, iNode, jNode, *theFrnMdl, Reff, kInit, theMaterials, y, x,
        shearDistI, doRayleigh, inclVertDisp, mass, maxIter, tol, kFactUplift);

    if (theElement == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // then add the singleFPBearing to the domain
    if (theTclDomain->addElement(theElement) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "singleFPBearing element: " << tag << endln;
      delete theElement;
      return TCL_ERROR;
    }
  }

  else if (ndm == 3) {
    // check space frame problem has 6 dof per node
    if (ndf != 6) {
      opserr << "WARNING invalid ndf: " << ndf;
      opserr << ", for space problem need 6 - singleFPBearing \n";
      return TCL_ERROR;
    }

    // check the number of arguments is correct
    if ((argc - eleArgStart) < 15) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: singleFPBearing eleTag iNode jNode frnMdlTag Reff kInit "
                "-P matTag -T matTag -My matTag -Mz matTag <-orient <x1 x2 x3> "
                "y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-inclVertDsip> "
                "<-mass m> <-iter maxIter tol>\n";
      return TCL_ERROR;
    }

    // get the id and end nodes
    int iNode, jNode, frnMdlTag, matTag, argi, i, j;
    int recvMat = 0;
    double Reff, kInit;
    double shearDistI = 0.0;
    int doRayleigh = 0;
    int inclVertDisp = 0;
    double mass = 0.0;
    int maxIter = 25;
    double tol = 1E-12;
    double kFactUplift = 1E-12;

    if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
      opserr << "WARNING invalid singleFPBearing eleTag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
      opserr << "WARNING invalid iNode\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
      opserr << "WARNING invalid jNode\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[4 + eleArgStart], &frnMdlTag) != TCL_OK) {
      opserr << "WARNING invalid frnMdlTag\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    FrictionModel *theFrnMdl = builder->getTypedObject<FrictionModel>(frnMdlTag);
    if (theFrnMdl == 0) {
      opserr << "WARNING friction model not found\n";
      opserr << "frictionModel: " << frnMdlTag << endln;
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5 + eleArgStart], &Reff) != TCL_OK) {
      opserr << "WARNING invalid Reff\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6 + eleArgStart], &kInit) != TCL_OK) {
      opserr << "WARNING invalid kInit\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }
    UniaxialMaterial *theMaterials[4];
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-P") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid axial matTag\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[0] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[0] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-T") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid torsional matTag\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[1] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[1] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-My") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid moment y matTag\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[2] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[2] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-Mz") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
          opserr << "WARNING invalid moment z matTag\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        theMaterials[3] = builder->getTypedObject<UniaxialMaterial>(matTag);
        if (theMaterials[3] == 0) {
          opserr << "WARNING material model not found\n";
          opserr << "uniaxialMaterial: " << matTag << endln;
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        recvMat++;
      }
    }
    if (recvMat != 4) {
      opserr << "WARNING wrong number of materials\n";
      opserr << "got " << recvMat << " materials, but want 4 materials\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // check for optional arguments
    Vector x(0);
    Vector y(3);
    y(0) = 0.0;
    y(1) = 1.0;
    y(2) = 0.0;
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-orient") == 0) {
        j = i + 1;
        int numOrient = 0;
        while (j < argc && strcmp(argv[j], "-shearDist") != 0 &&
               strcmp(argv[j], "-doRayleigh") != 0 &&
               strcmp(argv[j], "-inclVertDisp") != 0 &&
               strcmp(argv[j], "-mass") != 0 && strcmp(argv[j], "-iter") != 0 &&
               strcmp(argv[j], "-kFactUplift") != 0) {
          numOrient++;
          j++;
        }
        if (numOrient == 3) {
          argi = i + 1;
          double value;
          // read the y values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "singleFPBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              y(j) = value;
            }
          }
        } else if (numOrient == 6) {
          argi = i + 1;
          x.resize(3);
          double value;
          // read the x values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "singleFPBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              x(j) = value;
            }
          }
          // read the y values
          for (j = 0; j < 3; j++) {
            if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
              opserr << "WARNING invalid -orient value\n";
              opserr << "singleFPBearing element: " << tag << endln;
              return TCL_ERROR;
            } else {
              argi++;
              y(j) = value;
            }
          }
        } else {
          opserr << "WARNING insufficient arguments after -orient flag\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-shearDist") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &shearDistI) != TCL_OK) {
          opserr << "WARNING invalid -shearDist value\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-doRayleigh") == 0)
        doRayleigh = 1;
    }
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (strcmp(argv[i], "-inclVertDisp") == 0)
        inclVertDisp = 1;
    }
    for (i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-mass") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK) {
          opserr << "WARNING invalid -mass value\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 7 + eleArgStart; i < argc; i++) {
      if (i + 2 < argc && strcmp(argv[i], "-iter") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &maxIter) != TCL_OK) {
          opserr << "WARNING invalid maxIter\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i + 2], &tol) != TCL_OK) {
          opserr << "WARNING invalid tol\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
    for (int i = 7 + eleArgStart; i < argc; i++) {
      if (i + 1 < argc && strcmp(argv[i], "-kFactUplift") == 0) {
        if (Tcl_GetDouble(interp, argv[i + 1], &kFactUplift) != TCL_OK) {
          opserr << "WARNING invalid kFactUplift\n";
          opserr << "singleFPBearing element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }

    // now create the singleFPBearing
    theElement = new SingleFPSimple3d(
        tag, iNode, jNode, *theFrnMdl, Reff, kInit, theMaterials, y, x,
        shearDistI, doRayleigh, inclVertDisp, mass, maxIter, tol, kFactUplift);

    if (theElement == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "singleFPBearing element: " << tag << endln;
      return TCL_ERROR;
    }

    // then add the singleFPBearing to the domain
    if (builder->getDomain()->addElement(theElement) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "singleFPBearing element: " << tag << endln;
      delete theElement;
      return TCL_ERROR;
    }
  }

  else {
    opserr << "WARNING singleFPBearing command only works when ndm is 2 or 3, "
              "ndm: ";
    opserr << ndm << endln;
    return TCL_ERROR;
  }

  // if get here we have successfully created the singleFPBearing and added it
  // to the domain
  return TCL_OK;
}


// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the function invoked when the user invokes
// the frictionModel command in the interpreter.

#include <FrictionModel.h>
#include <tcl.h>
#include <elementAPI.h>
#include <BasicModelBuilder.h>
extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp,
                                       int cArg, int mArg, TCL_Char ** const argv,
                                       Domain *domain);

extern OPS_Routine OPS_Coulomb;
extern OPS_Routine OPS_VelDependent;
extern OPS_Routine OPS_VelDepMultiLinear;
extern OPS_Routine OPS_VelNormalFrcDep;
extern OPS_Routine OPS_VelPressureDep;

int
TclCommand_addFrictionModel(ClientData clientData, Tcl_Interp *interp,
                            int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  G3_Runtime *rt = G3_getRuntime(interp);

  // make sure there is a minimum number of arguments
  if (argc < 3) {
    opserr << "WARNING insufficient number of friction model arguments\n";
    opserr << "Want: frictionModel type tag <specific friction model args>\n";
    return TCL_ERROR;
  }

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, nullptr);

  // pointer to a friction model that will be added to the model builder
  FrictionModel *theFrnMdl = nullptr;

  // ----------------------------------------------------------------------------
  if (strcmp(argv[1], "Coulomb") == 0 || strcmp(argv[1], "Constant") == 0) {
    void *theFrn = OPS_Coulomb(rt, argc, argv);
    if (theFrn != 0)
      theFrnMdl = (FrictionModel *)theFrn;
    else
      return TCL_ERROR;
  }

  // ----------------------------------------------------------------------------
  if (strcmp(argv[1], "VelDependent") == 0 ||
      strcmp(argv[1], "VDependent") == 0) {
    void *theFrn = OPS_VelDependent(rt, argc, argv);
    if (theFrn != 0)
      theFrnMdl = (FrictionModel *)theFrn;
    else
      return TCL_ERROR;
  }

  // ----------------------------------------------------------------------------
  if (strcmp(argv[1], "VelDepMultiLinear") == 0 ||
      strcmp(argv[1], "VDependentMultiLinear") == 0) {
    void *theFrn = OPS_VelDepMultiLinear(rt, argc, argv);
    if (theFrn != 0)
      theFrnMdl = (FrictionModel *)theFrn;
    else
      return TCL_ERROR;
  }

  // ----------------------------------------------------------------------------
  if (strcmp(argv[1], "VelNormalFrcDep") == 0 ||
      strcmp(argv[1], "VNDependent") == 0) {
    void *theFrn = OPS_VelNormalFrcDep(rt, argc, argv);
    if (theFrn != 0)
      theFrnMdl = (FrictionModel *)theFrn;
    else
      return TCL_ERROR;
  }

  // ----------------------------------------------------------------------------
  if (strcmp(argv[1], "VelPressureDep") == 0 ||
      strcmp(argv[1], "VPDependent") == 0) {
    void *theFrn = OPS_VelPressureDep(rt, argc, argv);
    if (theFrn != 0)
      theFrnMdl = (FrictionModel *)theFrn;
    else
      return TCL_ERROR;
  }

  // ----------------------------------------------------------------------------
  if (theFrnMdl == nullptr) {
    opserr << "WARNING could not create friction model " << argv[1] << endln;
    return TCL_ERROR;
  }

  // now add the friction model to the modelBuilder
  if (builder->addTypedObject<FrictionModel>(theFrnMdl->getTag(), theFrnMdl) < 0) {
    opserr << "WARNING could not add friction model to the domain\n";
    opserr << *theFrnMdl << endln;
    delete theFrnMdl; // invoke the destructor, otherwise mem leak
    return TCL_ERROR;
  }

  return TCL_OK;
}

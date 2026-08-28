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
// Written: fmk
// Created: 03/01
//
#include <Parsing.h>
#include <Logging.h>
#include <string.h>
#include <Domain.h>
#include <ModelRegistry.h>

#include <Brick.h>
#include <Brick02.h>
// #define XARA_HAVE_H8E12
#ifdef XARA_HAVE_H8E12
# include <H8E12.h>
#endif
#include <BbarBrick.h>
#include <BbarBrickWithSensitivity.h>

#ifdef _FLBrick
#  include <FLBrick.h>
#endif

int
XaraElemCmd_H8(ClientData clientData, 
                         Tcl_Interp *interp,
                         ArgSize argc,
                         TCL_Char **const argv)
{
  const int eleArgStart = 1;

  ModelRegistry* builder = (ModelRegistry*)clientData;
  Domain* theTclDomain = builder->getDomain();

  // check the number of arguments is correct
  if ((argc - eleArgStart) < 11) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element Brick eleTag? Node1? Node2? Node3? Node4? Node5? "
              "Node6? Node7? Node 8? matTag?\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int BrickId;
  int matID;
  std::array<int, 8> node_tags;

  if (Tcl_GetInt(interp, argv[1 + eleArgStart], &BrickId) != TCL_OK) {
    opserr << "WARNING invalid Brick eleTag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2 + eleArgStart], &node_tags[0]) != TCL_OK) {
    opserr << "WARNING invalid Node1\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3 + eleArgStart], &node_tags[1]) != TCL_OK) {
    opserr << "WARNING invalid Node2\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4 + eleArgStart], &node_tags[2]) != TCL_OK) {
    opserr << "WARNING invalid Node3\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5 + eleArgStart], &node_tags[3]) != TCL_OK) {
    opserr << "WARNING invalid Node4\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6 + eleArgStart], &node_tags[4]) != TCL_OK) {
    opserr << "WARNING invalid Node5\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[7 + eleArgStart], &node_tags[5]) != TCL_OK) {
    opserr << "WARNING invalid Node6\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[8 + eleArgStart], &node_tags[6]) != TCL_OK) {
    opserr << "WARNING invalid Node7\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[9 + eleArgStart], &node_tags[7]) != TCL_OK) {
    opserr << "WARNING invalid Node8\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[10 + eleArgStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matTag\n";
    return TCL_ERROR;
  }


  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;
  if ((argc - eleArgStart) > 11) {
    if (Tcl_GetDouble(interp, argv[11 + eleArgStart], &b1) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      opserr << "Brick element: " << BrickId << "\n";
      return TCL_ERROR;
    }
  }

  if ((argc - eleArgStart) > 12) {
    if (Tcl_GetDouble(interp, argv[12 + eleArgStart], &b2) != TCL_OK) {
      opserr << "WARNING invalid b2\n";
      return TCL_ERROR;
    }
  }

  if ((argc - eleArgStart) > 13) {
    if (Tcl_GetDouble(interp, argv[13 + eleArgStart], &b3) != TCL_OK) {
      opserr << "WARNING invalid b3\n";
      return TCL_ERROR;
    }
  }

  // now create the Brick and add it to the Domain
  Element *theBrick = nullptr;
  if (strcmp(argv[1], "stdBrick") == 0) {
    if (getenv("XARA_BRICK02") != nullptr) {
      theBrick = new Brick02(BrickId, node_tags, *theMaterial, b1, b2, b3);
    } else {
      theBrick = new Brick(BrickId, node_tags, *theMaterial, b1, b2, b3);
    }
  }
  else if (strcmp(argv[1], "bbarBrickWithSensitivity") == 0) {
    theBrick = new BbarBrickWithSensitivity(BrickId, 
                                            node_tags[0], node_tags[1], node_tags[2], node_tags[3],
                                            node_tags[4], node_tags[5], node_tags[6], node_tags[7],
                                            *theMaterial, b1, b2, b3);
  } else if (strcmp(argv[1], "bbarBrick") == 0) {
    theBrick = new BbarBrick(BrickId, node_tags, *theMaterial, b1, b2, b3);
  }
#ifdef XARA_HAVE_H8E12
  else if (strcmp(argv[1], "H8E12") == 0) {
    theBrick = new H8E12(BrickId, node_tags, *theMaterial, b1, b2, b3);
  }
#endif
#ifdef _FLBrick
  else if (strcmp(argv[1], "flBrick") == 0) {
    theBrick = new FLBrick(BrickId, node_tags[0], node_tags[1], node_tags[2], node_tags[3], node_tags[4], node_tags[5],
                           node_tags[6], node_tags[7], *theMaterial, b1, b2, b3);
  }
#endif

  else {
    opserr << "WARNING brick element " << argv[1] << " type not recognized\n";
    return TCL_ERROR;
  }


  if (theTclDomain->addElement(theBrick) == false) {
    opserr << "WARNING could not add element to the domain\n";
    delete theBrick;
    return TCL_ERROR;
  }

  return TCL_OK;
}

//
// Description: This file contains the implementation of
//    XaraElemCmd_H8UP() ,
//    TclBasicBuilder_addTwentyEightNodeBrickUP(),
//    TclBasicBuilder_addBBarBrickUP()
//
//
// Jinchi Lu and Zhaohui Yang (May 2004)
//
#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <BrickUP.h>
#include <BBarBrickUP.h>
#include <Twenty_Node_Brick.h>
#include <Twenty_Eight_Node_BrickUP.h>

int
XaraElemCmd_H20(ClientData clientData, 
                Tcl_Interp *interp,
                ArgSize argc,
                TCL_Char ** const argv)
{
  ModelRegistry *builder = (ModelRegistry*)clientData;
  Domain* theTclDomain = builder->getDomain();

  if (builder == 0 || clientData == 0) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }

  if (builder->getNDM() != 3) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with 20NodeBrick element\n";
    return TCL_ERROR;
  }
  // check the number of arguments is correct
  int argStart = 2;
  if ((argc - argStart) < 22) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element 20NodeBrick eleTag? N1? N2? N3? N4? N5? N6? N7? "
              "N8? N9? N10? N11? N12? N13? N14? N15? N16? N17? N18? N19? N20? "
              "matTag? <b1? b2? b3?>\n";
    return TCL_ERROR;
  }
  // get the id and end nodes
  int brickId;
  std::array<int, 20> Nod{};
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;
  if (Tcl_GetInt(interp, argv[argStart], &brickId) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "invalid element tag " << argv[argStart] << "\n";
    return TCL_ERROR;
  }
  for (int i = 0; i < 20; i++)
    if (Tcl_GetInt(interp, argv[1 + argStart + i], &(Nod[i])) != TCL_OK) {
      opserr << "WARNING invalid Node number\n";
      return TCL_ERROR;
    }

  int matID;
  if (Tcl_GetInt(interp, argv[21 + argStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matID\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) >= 23) {
    if (Tcl_GetDouble(interp, argv[22 + argStart], &b1) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 24) {
    if (Tcl_GetDouble(interp, argv[23 + argStart], &b2) != TCL_OK) {
      opserr << "WARNING invalid b2\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 25) {
    if (Tcl_GetDouble(interp, argv[24 + argStart], &b3) != TCL_OK) {
      opserr << "WARNING invalid b3\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }

  // now create the brick and add it to the Domain
  Twenty_Node_Brick *theTwentyNodeBrick =
      new Twenty_Node_Brick(brickId, 
                            Nod,
                            *theMaterial, b1, b2, b3);

  if (theTclDomain->addElement(theTwentyNodeBrick) == false) {
    opserr << "WARNING could not add element to the domain\n";
    delete theTwentyNodeBrick;
    return TCL_ERROR;
  }
  // if get here we have successfully created the element and added it to the
  // domain
  return TCL_OK;
}


int
XaraElemCmd_H8UP(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
                           TCL_Char ** const argv)
{
  ModelRegistry *builder = (ModelRegistry*)clientData;
  Domain* theTclDomain = builder->getDomain();

  if (builder == 0 || clientData == 0) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }

  if (builder->getNDM() != 3 || builder->getNDF() != 4) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with QuadUP element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 15) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element brickUP eleTag? N1? N2? N3? N4? N5? N6? N7? N8? "
              "matTag? bulk? rhof? perm_x? perm_y? perm_z? <b1? b2? b3?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int brickUPId, Nod[8], matID;
  double bk, r, perm1, perm2, perm3;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &brickUPId) != TCL_OK) {
    opserr << "WARNING invalid brickUP eleTag" << "\n";
    return TCL_ERROR;
  }

  for (int i = 0; i < 8; i++)
    if (Tcl_GetInt(interp, argv[1 + argStart + i], &(Nod[i])) != TCL_OK) {
      opserr << "WARNING invalid Node number\n";
      return TCL_ERROR;
    }

  if (Tcl_GetInt(interp, argv[9 + argStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matID\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[10 + argStart], &bk) != TCL_OK) {
    opserr << "WARNING invalid fluid bulk modulus\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[11 + argStart], &r) != TCL_OK) {
    opserr << "WARNING invalid fluid mass density\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[12 + argStart], &perm1) != TCL_OK) {
    opserr << "WARNING invalid permeability_x\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[13 + argStart], &perm2) != TCL_OK) {
    opserr << "WARNING invalid permeability_y\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[14 + argStart], &perm3) != TCL_OK) {
    opserr << "WARNING invalid permeability_z\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) >= 16) {
    if (Tcl_GetDouble(interp, argv[15 + argStart], &b1) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 17) {
    if (Tcl_GetDouble(interp, argv[16 + argStart], &b2) != TCL_OK) {
      opserr << "WARNING invalid b2\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 18) {
    if (Tcl_GetDouble(interp, argv[17 + argStart], &b3) != TCL_OK) {
      opserr << "WARNING invalid b3\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }

  // now create the brickUP and add it to the Domain
  BrickUP *theBrickUP = new BrickUP(
      brickUPId, Nod[0], Nod[1], Nod[2], Nod[3], Nod[4], Nod[5], Nod[6], Nod[7],
      *theMaterial, bk, r, perm1, perm2, perm3, b1, b2, b3);


  if (theTclDomain->addElement(theBrickUP) == false) {
    opserr << "WARNING could not add element to the domain\n";
    delete theBrickUP;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclBasicBuilder_addTwentyEightNodeBrickUP(ClientData clientData, Tcl_Interp *interp,
                                          int argc, TCL_Char ** const argv)

{
  ModelRegistry *builder = (ModelRegistry*)clientData;
  Domain* theTclDomain = builder->getDomain();


  if (builder->getNDM() != 3) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with 20_8_BrickUP element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 27) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element 20_8_BrickUP eleTag? N1? N2? N3? N4? N5? N6? N7? "
              "N8? N9? N10? N11? N12? N13? N14? N15? N16? N17? N18? N19? N20? "
              "matTag? bulk? rhof? perm_x? perm_y? perm_z? <b1? b2? b3?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int brickUPId, Nod[20], matID;
  double bk, r, perm1, perm2, perm3;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &brickUPId) != TCL_OK) {
    opserr << "WARNING invalid 20_8_BrickUP eleTag" << "\n";
    return TCL_ERROR;
  }

  for (int i = 0; i < 20; i++)
    if (Tcl_GetInt(interp, argv[1 + argStart + i], &(Nod[i])) != TCL_OK) {
      opserr << "WARNING invalid Node number\n";
      return TCL_ERROR;
    }

  if (Tcl_GetInt(interp, argv[21 + argStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matID\n";
    opserr << "20_8_BrickUP element: " << brickUPId << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[22 + argStart], &bk) != TCL_OK) {
    opserr << "WARNING invalid fluid bulk modulus\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[23 + argStart], &r) != TCL_OK) {
    opserr << "WARNING invalid fluid mass density\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[24 + argStart], &perm1) != TCL_OK) {
    opserr << "WARNING invalid permeability_x\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[25 + argStart], &perm2) != TCL_OK) {
    opserr << "WARNING invalid permeability_y\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[26 + argStart], &perm3) != TCL_OK) {
    opserr << "WARNING invalid permeability_z\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) >= 28) {
    if (Tcl_GetDouble(interp, argv[27 + argStart], &b1) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 29) {
    if (Tcl_GetDouble(interp, argv[28 + argStart], &b2) != TCL_OK) {
      opserr << "WARNING invalid b2\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 30) {
    if (Tcl_GetDouble(interp, argv[29 + argStart], &b3) != TCL_OK) {
      opserr << "WARNING invalid b3\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);

  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }

  // now create the brickUP and add it to the Domain
  TwentyEightNodeBrickUP *theTwentyEightNodeBrickUP =
      new TwentyEightNodeBrickUP(
          brickUPId, Nod[0], Nod[1], Nod[2], Nod[3], Nod[4], Nod[5], Nod[6],
          Nod[7], Nod[8], Nod[9], Nod[10], Nod[11], Nod[12], Nod[13], Nod[14],
          Nod[15], Nod[16], Nod[17], Nod[18], Nod[19], *theMaterial, bk, r,
          perm1, perm2, perm3, b1, b2, b3);

  if (theTclDomain->addElement(theTwentyEightNodeBrickUP) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "20_8_BrickUP element: " << brickUPId << "\n";
    delete theTwentyEightNodeBrickUP;
    return TCL_ERROR;
  }

  // if get here we have successfully created the element and added it to the
  // domain
  return TCL_OK;
}



/*  *****************************************************************************
    BBAR  BRICK  U_P
    *****************************************************************************
 */

int
TclBasicBuilder_addBBarBrickUP(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
                               TCL_Char ** const argv)
{
  ModelRegistry *builder = (ModelRegistry*)clientData;
  Domain* theTclDomain = builder->getDomain();


  if (builder == nullptr || clientData == nullptr) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }

  if (builder->getNDM() != 3 || builder->getNDF() != 4) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with QuadUP element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 15) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element BBarBrickUP eleTag? N1? N2? N3? N4? N5? N6? N7? "
              "N8? matTag? bulk? rhof? perm_x? perm_y? perm_z? <b1? b2? b3?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int BBarBrickUPId, Nod[8], matID;
  double bk, r, perm1, perm2, perm3;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &BBarBrickUPId) != TCL_OK) {
    opserr << "WARNING invalid BBarBrickUP eleTag" << "\n";
    return TCL_ERROR;
  }

  for (int i = 0; i < 8; i++)
    if (Tcl_GetInt(interp, argv[1 + argStart + i], &(Nod[i])) != TCL_OK) {
      opserr << "WARNING invalid Node number\n";
      opserr << "BBarBrickUP element: " << BBarBrickUPId << "\n";
      return TCL_ERROR;
    }

  if (Tcl_GetInt(interp, argv[9 + argStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matID\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[10 + argStart], &bk) != TCL_OK) {
    opserr << "WARNING invalid fluid bulk modulus\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[11 + argStart], &r) != TCL_OK) {
    opserr << "WARNING invalid fluid mass density\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[12 + argStart], &perm1) != TCL_OK) {
    opserr << "WARNING invalid permeability_x\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[13 + argStart], &perm2) != TCL_OK) {
    opserr << "WARNING invalid permeability_y\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[14 + argStart], &perm3) != TCL_OK) {
    opserr << "WARNING invalid permeability_z\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) >= 16) {
    if (Tcl_GetDouble(interp, argv[15 + argStart], &b1) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      opserr << "BBarBrickUP element: " << BBarBrickUPId << "\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 17) {
    if (Tcl_GetDouble(interp, argv[16 + argStart], &b2) != TCL_OK) {
      opserr << "WARNING invalid b2\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 18) {
    if (Tcl_GetDouble(interp, argv[17 + argStart], &b3) != TCL_OK) {
      opserr << "WARNING invalid b3\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }

  // now create the BBarBrickUP and add it to the Domain
  BBarBrickUP *theBBarBrickUP = new BBarBrickUP(
      BBarBrickUPId, Nod[0], Nod[1], Nod[2], Nod[3], Nod[4], Nod[5], Nod[6],
      Nod[7], *theMaterial, bk, r, perm1, perm2, perm3, b1, b2, b3);


  if (theTclDomain->addElement(theBBarBrickUP) == false) {
    opserr << "WARNING could not add element to the domain\n";
    delete theBBarBrickUP;
    return TCL_ERROR;
  }

  // if get here we have successfully created the element and added it to the
  // domain
  return TCL_OK;
}
